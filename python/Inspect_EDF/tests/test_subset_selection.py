"""The shared participant selector inside the tools (invariant 6 of the plan): in 'Subset only' mode
only the chosen participants are touched; 'Skip already processed' still applies to a subset, with a
visible warning naming the participants it leaves out.

Each test works on a private copy of the chain's data folder, so the golden comparison is unaffected.
"""
import shutil

import pytest

import minidb
import pipeline
from nbdriver import NotebookSession
from pipeline import select_subset


@pytest.fixture
def data_copy(chain, tmp_path):
    dst = tmp_path / 'data'
    shutil.copytree(chain['data'], dst)
    return dst


def mtimes(folder, pattern):
    return {f.name: f.stat().st_mtime_ns for f in folder.rglob(pattern)}


@pytest.mark.tool5
def test_tool5_subset(chain, data_copy):
    p = minidb.PARTICIPANTS[chain['dataset']]
    per_file = data_copy / 'derivatives' / 'features_macrostructure'
    before = mtimes(per_file, '*_sleep_metrics.tsv')
    assert len(before) == 4

    with NotebookSession(pipeline.T5) as nb:
        nb.pick('fc_data', data_copy)
        nb.pick('fc_subj', data_copy / 'participants.tsv')
        nb.set('dd_join_col', 'participant_id')
        nb.click('btn_scan')
        # every participant is already processed: in 'All' mode nothing would run
        assert '0</b> will run' in nb.value('participant_selector.lbl_summary.value')

        select_subset(nb, [p[0]], skip=True)                     # skip ticked -> warned, nothing run
        warn = nb.value('participant_selector.html_warn.value')
        assert 'will NOT be processed' in warn and p[0] in warn
        nb.click('btn_run')
        assert 'will NOT be processed' in nb.widget_text()     # repeated at the top of the run log
        assert mtimes(per_file, '*_sleep_metrics.tsv') == before

        nb.set('participant_selector.cb_skip', False)          # untick -> p[0] reprocessed, only it
        assert nb.value('participant_selector.html_warn.value') == ''
        nb.click('btn_run')
    after = mtimes(per_file, '*_sleep_metrics.tsv')
    changed = sorted(k for k in after if after[k] != before[k])
    assert changed == [f'{p[0]}_sleep_metrics.tsv']
    # the database table is rebuilt from every per-recording file: still the 4 recordings
    gtable = (data_copy / 'reports_features_macrostructure' / 'global_sleep_metrics.tsv').read_text()
    assert all(pid in gtable for pid in p)


@pytest.mark.tool6
def test_tool6_subset(chain, data_copy):
    p = minidb.PARTICIPANTS[chain['dataset']]
    reports = data_copy / 'reports_quality_overview'
    before = mtimes(reports, '*_quality_metrics.tsv')
    assert len(before) == 4
    log = pipeline.run_tool6(data_copy, subset=[p[2]], skip=True)          # already processed: warned
    assert 'will NOT be processed' in log
    assert mtimes(reports, '*_quality_metrics.tsv') == before
    pipeline.run_tool6(data_copy, subset=[p[2]], skip=False)               # reprocessed, only it
    after = mtimes(reports, '*_quality_metrics.tsv')
    assert sorted(k for k in after if after[k] != before[k]) == [f'{p[2]}_quality_metrics.tsv']
    summary = (reports / 'quality_summary.tsv').read_text(encoding='utf-8')
    assert all(pid in summary for pid in p)                                # still the whole database
