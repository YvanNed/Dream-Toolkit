"""Small runs of the preparation tools 1, 3 and 4 on a fresh synthetic database, with known answers.

    tool 1   header inspection: 4 files x 4 channels at 256 Hz, the clipped Fp1 of sim01 (physical range
             +-50 uV) reported as a bad dynamic range, every header anonymous
    tool 3   hypnogram remapping: a raw Compumedics-style export (0-5 / W / R / ?) is written from the
             generated hypnogram; the remapped file must give back the AASM stages epoch by epoch, and the
             mid-night '?' must be listed for review
    tool 4   event harmonization: the raw labels are listed and event_remap.json gets the canonical names,
             including a French label (suggested) and a label the user ticks 'ignore' (saved as null)

Tools 0 and 2 stay at the quick level (they execute): too interactive for a useful run.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import minidb
from nbdriver import NotebookSession

T1 = 'tools/1_inspect_edf_voila.ipynb'
T3 = 'tools/3_remap_hypno_voila.ipynb'
T4 = 'tools/4_remap_events_edf_voila.ipynb'
IDS = minidb.PARTICIPANTS['synthetic']


@pytest.fixture(scope='module')
def fresh_db(tmp_path_factory):
    """The synthetic recordings, built once for the three tools (~20 s); each tool writes to its own place."""
    return minidb.build(tmp_path_factory.mktemp('tools134'), 'synthetic')


# ----------------------------------------------------------------------------------------------- tool 1
@pytest.mark.tool1
def test_tool1_inspects_the_headers(fresh_db):
    with NotebookSession(T1) as nb:
        nb.pick('chooser', fresh_db)
        nb.click('run_button')
    out = fresh_db / 'summary_inspection'
    full = pd.read_csv(out / 'FULL_summary_table_edf.tsv', sep='\t', dtype=str)
    assert sorted(full['subject'].unique()) == IDS
    channels = {fid: sorted(chs) for fid, chs in full.groupby('subject')['channel']}
    assert channels == {fid: ['A2', 'C3', 'Fp1', 'O1'] for fid in IDS}
    assert set(full['sampling_frequency']) == {'256'}

    bad_range = pd.read_csv(out / 'EEG_bad_dynamic_range_edf.tsv', sep='\t', dtype=str)
    assert bad_range[['subject', 'channel']].values.tolist() == [['sim01_clipping', 'Fp1']]

    anon = pd.read_csv(out / 'anonymization_check_edf.tsv', sep='\t')
    assert len(anon) == 4 and anon['header_anonymized'].all()


# ----------------------------------------------------------------------------------------------- tool 3
RAW_CODE = {'W': 'W', 'N1': '1', 'N2': '2', 'R': 'R'}      # the Compumedics export codes


@pytest.mark.tool3
def test_tool3_remaps_a_raw_export_to_aasm(fresh_db):
    """The raw export codes N3 as '3' or '4' (old R&K stage 4), starts with an unscored '?' and, in sim01,
    holds one '?' in the middle of the night. Expected after remapping: the generated stages, with both
    '?' as W (the default rule); the mid-night one is listed for review in mid_uncertain_epochs_to_verify.tsv."""
    truth, mid = {}, 200
    for fid in IDS:
        stages = np.loadtxt(fresh_db / f'{fid}_Hypnogram_remapped.txt', dtype=str)
        raw = [RAW_CODE.get(s, '3' if i % 2 else '4') for i, s in enumerate(stages)]
        raw[0] = '?'
        expected = stages.copy()
        expected[0] = 'W'
        if fid == IDS[0]:
            raw[mid] = '?'
            expected[mid] = 'W'
        (fresh_db / f'{fid}_Hypnogram_Export.txt').write_text('\n'.join(raw) + '\n', encoding='utf-8')
        truth[fid] = expected

    with NotebookSession(T3) as nb:
        nb.pick('chooser', fresh_db)
        nb.set('suffix_box', '_Hypnogram_Export.txt')
        nb.set('output_suffix_box', '_Hypnogram_tool3.txt')     # keep the generated truth untouched
        nb.click('btn_scan')
        nb.click('btn_run_remap')
        nb.click_labelled('Save remapping')
        nb.click('btn_run_save')
        nb.click_labelled('Save files')
        nb.click('btn_run_verify')
        text = nb.widget_text()
    assert 'Remapping successful' in text, text[-2000:]

    for fid in IDS:
        remapped = np.loadtxt(fresh_db / f'{fid}_Hypnogram_tool3.txt', dtype=str)
        assert remapped.tolist() == truth[fid].tolist(), fid

    review = pd.read_csv(fresh_db / 'mid_uncertain_epochs_to_verify.tsv', sep='\t', dtype={'participant_id': str})
    assert review[['participant_id', 'epoch_index']].values.tolist() == [[IDS[0], mid]]


# ----------------------------------------------------------------------------------------------- tool 4
@pytest.mark.tool4
def test_tool4_harmonizes_the_event_labels(fresh_db):
    """sim02's event CSV gains a French label ('Hypopnée', suggested from the French rules) and a 'Lights Off'
    marker the user ticks 'ignore'. A mapping already in event_remap.json and absent from the files is kept."""
    csv = fresh_db / f'{IDS[1]}_event_xml.csv'
    lines = csv.read_text(encoding='utf-8').splitlines()
    n = len(lines) - 1
    lines += [f'{n},Hypopnée,5000.0,12.0', f'{n + 1},Lights Off,0.5,0.0']
    csv.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    remap_json = fresh_db / 'config_param' / 'event_remap.json'
    remap_json.write_text(json.dumps({'Central Apnea': 'apnea_central'}), encoding='utf-8')

    with NotebookSession(T4) as nb:
        nb.pick('chooser', fresh_db)
        nb.click('run_scan_button')
        listed = nb.value('{k: sorted(v) for k, v in STATE["label_files"].items()}')
        nb.run('STATE["rows"]["Lights Off"][1].value = True')          # the user ticks 'ignore'
        nb.click('validate_button')
        nb.click('preview_save_button')
        nb.click('verify_button')
        text = nb.widget_text()

    synthetic_labels = ['Arousal (ARO RES)', 'Arousal (ARO SPONT)', 'Hypopnea', 'Obstructive Apnea',
                        'SpO2 desaturation']
    assert {k: v for k, v in listed.items() if k in synthetic_labels} == {k: IDS for k in synthetic_labels}
    assert listed['Hypopnée'] == [IDS[1]] and listed['Lights Off'] == [IDS[1]]

    expected = {k: minidb.EVENT_REMAP[k] for k in synthetic_labels}
    expected.update({'Hypopnée': 'hypopnea', 'Lights Off': None, 'Central Apnea': 'apnea_central'})
    assert json.loads(remap_json.read_text(encoding='utf-8')) == expected
    assert 'All raw labels are mapped' in text
