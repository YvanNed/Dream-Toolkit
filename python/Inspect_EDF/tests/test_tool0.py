"""Tool 0 (one participant, the whole pipeline) must give the SAME results as the batch tools on the same
participant: it shares their code (tools/*_lib.py, see tools/plan_tool0_redesign.md). Each step of tool 0 is
compared with the chain outputs of the batch tool it mirrors, on a copy of the study inputs (tool 0 writes
its report and annotations into the study folder, never into the chain's).

Phase 1: step 1 (loading: channels, hypnogram, events, context channels) against tool 7, step 2 (header)
against tool 1's known answer, and the explorer's writes (rescoring, annotations, inspection report)."""
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import minidb
import pipeline
from nbdriver import NotebookSession

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import edf_io_lib as EIO                    # noqa: E402
import qc_rejected_epochs_lib as L          # noqa: E402

pytestmark = pytest.mark.tool0

T0 = 'tools/0_live_explore_1file_voila.ipynb'
INPUT_SUFFIXES = ('.edf', '.txt', '.csv', '.xml', '.tsv')


def study_copy(chain, dest, extra=()):
    """The study INPUTS of the chain (recordings, companions, config_param, participant table) copied to
    `dest`, plus the chain output folders listed in `extra` (relative paths)."""
    data = Path(chain['data'])
    dest.mkdir(parents=True, exist_ok=True)
    for f in data.iterdir():
        if f.is_file() and f.suffix.lower() in INPUT_SUFFIXES:
            shutil.copy2(f, dest / f.name)
    shutil.copytree(data / 'config_param', dest / 'config_param')
    for rel in extra:
        src = data / rel
        if src.is_dir():
            shutil.copytree(src, dest / rel)
        else:
            (dest / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest / rel)
    return dest


def load_participant(nb, study, fid):
    nb.pick('fc_data', study)
    nb.set('dd_part', fid)
    nb.click('btn_run1')
    assert nb.value("STATUS['1']") == 'up_to_date', nb.widget_text()[-2000:]


def tool7_fif(chain, fid):
    return next(Path(chain['data'], 'derivatives', 'raw_epo').rglob(f'{fid}_all-epo.fif'))


def test_step1_reads_the_recording_like_tool7(chain, tmp_path):
    """Same epochs, same stages, same EEG channels, same mapped events as tool 7 for every participant."""
    import mne
    study = study_copy(chain, tmp_path / 'study')
    with NotebookSession(T0) as nb:
        for fid in minidb.PARTICIPANTS[chain['dataset']]:
            load_participant(nb, study, fid)
            ep = mne.read_epochs(str(tool7_fif(chain, fid)), preload=False, verbose=False)
            meta = ep.metadata.reset_index(drop=True)
            assert nb.value("STATE['n']") == len(meta), fid
            assert nb.value("[str(s) for s in STATE['hypno_work']]") == meta['stage'].astype(str).tolist(), fid
            # tool 7 re-references to M2 then drops it; tool 0 shows the EEG as recorded, M2 included
            assert set(nb.value("STATE['P']['ch_names']")) - {'M2'} == set(ep.ch_names), fid
            on7 = pd.read_csv(tool7_fif(chain, fid).with_name(f'{fid}_event_onsets.tsv'), sep='\t')
            on0 = pd.DataFrame(nb.value("STATE['onsets'][STATE['onsets']['mapped']][['type', 'onset_s']]"
                                        ".to_dict('records')"))
            key = lambda df: sorted(zip(df['type'].astype(str), df['onset_s'].astype(float).round(3)))
            assert key(on0) == key(on7), fid


def test_step2_flags_the_narrow_dynamic_range(chain, tmp_path):
    """Tool 1's known answer (test_tools_1_3_4): the clipping recording exports Fp1 at ±50 µV."""
    study = study_copy(chain, tmp_path / 'study')
    fid = minidb.PARTICIPANTS[chain['dataset']][0]
    with NotebookSession(T0) as nb:
        nb.pick('fc_data', study)
        nb.set('dd_part', fid)
        nb.click('btn_run2')
        assert nb.value("STATUS['2']") == 'up_to_date'
        content = nb.value("REPORT_BLOCKS['2'][0]")
    if chain['dataset'] == 'synthetic':
        assert 'EEG <b>Fp1</b>: dynamic range 100 µV' in content
        assert 'header anonymized' in content


def test_explorer_writes_only_its_own_files(chain, tmp_path):
    """Rescoring and annotations are written beside the recording in the documented files (the annotation in
    tool 8's format, readable by tool 7's event chain); the inspection log + report go to
    reports_live_explore/; nothing else in the study is created or modified."""
    study = study_copy(chain, tmp_path / 'study')
    fid = minidb.PARTICIPANTS[chain['dataset']][0]
    before = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in study.rglob('*') if p.is_file()}
    with NotebookSession(T0) as nb:
        load_participant(nb, study, fid)
        nb.run('goto_epoch(10)')
        old = nb.value("str(STATE['hypno_work'][10])")
        new = 'N3' if old != 'N3' else 'N2'
        nb.run(f"stage_btns[{new!r}].click()")
        assert nb.value("STATUS['1']") == 'up_to_date'
        nb.click('btn_run1')                             # a reload would lose the unsaved rescoring:
        assert 'not saved' in nb.value('html_guard.value')            # refused once, with a warning
        assert nb.value("str(STATE['hypno_work'][10])") == new
        nb.click('btn_save_hypno')
        assert nb.value('html_guard.value') == '' and nb.value('unsaved_rescoring()') == 0
        nb.run('goto_epoch(20)')
        nb.run("dd_annot_label.value = 'arousal_spontaneous'")
        nb.set('ft_annot_offset', 12.0)
        nb.click('btn_annot_add')
        nb.set('dd_signal', 'bp')                        # a display filter: nothing written
        assert nb.value("STATE['display_mode']") == 'bp'
    after = {p for p in study.rglob('*') if p.is_file()}
    hyp = study / f'{fid}_Hypnogram_remapped.txt'
    expected = {hyp.with_name(f'{hyp.stem}_rescored.txt'), hyp.with_name(f'{hyp.stem}_rescore_log.tsv'),
                study / f'{fid}_manual_events.tsv',
                study / 'reports_live_explore' / f'{fid}_inspection_log.tsv',
                study / 'reports_live_explore' / f'{fid}_inspection_report.html'}
    assert after - set(before) == expected
    assert [p for p, s in before.items() if (p.stat().st_size, p.stat().st_mtime_ns) != s] == []
    # the rescoring
    rescored = np.loadtxt(hyp.with_name(f'{hyp.stem}_rescored.txt'), dtype=str)
    original = np.loadtxt(hyp, dtype=str)
    assert rescored[10] == new and (rescored != original).sum() == 1
    log = pd.read_csv(hyp.with_name(f'{hyp.stem}_rescore_log.tsv'), sep='\t')
    assert log[['epoch_idx', 'original_label', 'new_label']].values.tolist() == [[10, old, new]]
    # the annotation: tool 8's table, read back by the shared loader and by tool 7's event chain
    man = L.load_manual_events(study, fid)
    assert man[['type', 'onset_s', 'epoch_idx', 'source']].values.tolist() == \
        [['arousal_spontaneous', 612.0, 20, 'manual_tool0']]
    ev, src = EIO.load_events(study / f'{fid}.edf', include_manual=True)
    assert src.endswith('+manual') and (ev['Source'] == 'manual').sum() == 1
    # the inspection log: one row per action
    steps = pd.read_csv(study / 'reports_live_explore' / f'{fid}_inspection_log.tsv', sep='\t')
    assert steps['action'].tolist() == ['recording loaded', 'rescored hypnogram saved', 'annotation added']


def test_step3_lists_every_scored_event(chain, tmp_path):
    """Step 3 counts every scored event of the companion under its tool-4 label, and per sleep stage."""
    study = study_copy(chain, tmp_path / 'study')
    fid = minidb.PARTICIPANTS[chain['dataset']][0]
    remap = json.loads((study / 'config_param' / 'event_remap.json').read_text(encoding='utf-8'))
    ev, _src = EIO.load_events(study / f'{fid}.edf')
    expected = ev['Name'].map(lambda n: remap.get(n) or n).value_counts().to_dict()
    with NotebookSession(T0) as nb:
        load_participant(nb, study, fid)
        nb.click('btn_run3')
        assert nb.value("STATUS['3']") == 'up_to_date'
        counts = nb.value("STATE['onsets']['type'].value_counts().to_dict()")
        html = nb.value("''.join(REPORT_BLOCKS['3'])")
    assert counts == expected and sum(counts.values()) == len(ev)
    for typ in expected:
        assert typ in html


@pytest.mark.tool3bis
def test_step3_gssc_equals_tool3bis(chain, tmp_path):
    """Tool 0's automatic scoring is tool 3bis's: the scoring read from tool 3bis's outputs, and the scoring
    redone in tool 0 with the same settings, give tool 3bis's stages, probabilities and kappa."""
    pytest.importorskip('gssc')
    if chain['dataset'] != 'synthetic':
        pytest.skip('one dataset is enough: GSSC is deterministic and the code path is the same')
    fid = minidb.PARTICIPANTS['synthetic'][0]
    study = study_copy(chain, tmp_path / 'study')
    with NotebookSession('tools/3bis_autoscoring_gssc_voila.ipynb') as nb:
        nb.pick('fc_data', study)
        nb.click('btn_scan')
        pipeline.select_subset(nb, [fid], skip=False)
        nb.click('btn_run')
    dens = pd.read_csv(next(study.rglob(f'{fid}_hypnodensity.tsv')), sep='\t')
    kappa_3bis = float(pd.read_csv(next(study.rglob(f'{fid}_autoscoring_comparison.tsv')), sep='\t')['kappa'].iloc[0])
    probs_3bis = dens[[f'p_{s}' for s in ('W', 'N1', 'N2', 'N3', 'R')]].values
    with NotebookSession(T0) as nb:
        load_participant(nb, study, fid)
        nb.click('btn_run3')
        for button in ('btn_gssc_read', 'btn_gssc_run'):          # read tool 3bis's, then score here
            nb.click(button)
            assert nb.value("list(STATE['gssc']['stages'])") == dens['stage_gssc'].tolist(), button
            probs = np.array(nb.value("STATE['gssc']['probs'].tolist()"))
            np.testing.assert_allclose(probs, probs_3bis, atol=6e-5)   # the TSV keeps 4 decimals
            assert abs(nb.value("STATE['gssc']['result']['kappa']") - kappa_3bis) < 6e-5, button
        assert nb.value("box_gssc_nav.layout.display") == ''         # the explorer offers the GSSC layer


@pytest.mark.tool5
def test_step4_equals_tool5(chain, tmp_path):
    """Step 4 is tool 5's computation: every registered metric and the lights source of every participant
    equal tool 5's {file_id}_sleep_metrics.tsv (4 decimals, as tool 5 writes them)."""
    study = study_copy(chain, tmp_path / 'study')
    with NotebookSession(T0) as nb:
        for fid in minidb.PARTICIPANTS[chain['dataset']]:
            load_participant(nb, study, fid)
            nb.click('btn_run4')
            assert nb.value("STATUS['4']") == 'up_to_date', nb.widget_text()[-1500:]
            mine = nb.value("{k: (None if v is None or not np.isfinite(float(v)) else round(float(v), 4))"
                            " for k, v in STATE['macro']['metrics'].items()}")
            lights = nb.value("STATE['macro']['lights']['source']")
            saved = pd.read_csv(next(Path(chain['data'], 'derivatives', 'features_macrostructure')
                                     .rglob(f'{fid}_sleep_metrics.tsv')), sep='\t')
            ref = dict(zip(saved['metric'], saved['value']))
            assert lights == ref['lights_source'], fid
            assert mine, fid
            for key, value in mine.items():
                expected = float(ref[key])
                if value is None:
                    assert np.isnan(expected), (fid, key)
                else:
                    assert abs(value - expected) < 1e-9, (fid, key, value, expected)


def test_context_channels_match_tool7(chain, tmp_path):
    """The EOG / EMG context traces of the explorer are tool 7's _context-epo.fif companion, sample for
    sample. The test montage has no EOG / EMG, so two EEG channels are declared as context channels."""
    if chain['dataset'] != 'synthetic':
        pytest.skip('one dataset is enough: the code path does not depend on the signal')
    fid = minidb.PARTICIPANTS['synthetic'][0]
    study = study_copy(chain, tmp_path / 'study', extra=['reports_quality_overview/quality_summary.tsv'])
    cfg_path = study / 'config_param' / 'remap_reref_persubject.json'
    cfg = json.loads(cfg_path.read_text(encoding='utf-8'))
    cfg[fid]['context_channels'] = {'eog_left': 'Fp1', 'emg': 'O1'}
    cfg_path.write_text(json.dumps(cfg, indent=2), encoding='utf-8')
    pipeline.run_tool7(study, events=False, subset=[fid], skip=False)
    ctx7 = L.load_context_epochs(next(study.rglob(f'{fid}_all-epo.fif')).parent, fid)
    assert ctx7 is not None and ctx7['labels'] == ['EOG-L', 'EMG']
    with NotebookSession(T0) as nb:
        load_participant(nb, study, fid)
        assert nb.value("STATE['ctx']['labels']") == ['EOG-L', 'EMG']
        nb.run(f"np.save(r'{tmp_path / 'ctx0.npy'}', STATE['ctx']['data_uV'])")
    ctx0 = np.load(tmp_path / 'ctx0.npy')
    assert ctx0.shape == ctx7['data_uV'].shape
    # the .fif stores float32 (~7 significant digits): 1e-6 relative + 1 nV absolute is that precision
    np.testing.assert_allclose(ctx0, ctx7['data_uV'], rtol=1e-6, atol=1e-3)
