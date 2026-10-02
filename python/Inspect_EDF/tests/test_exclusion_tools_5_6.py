"""Step 5 of the exclusion refactor: tools 5 and 6 flag the excluded participants in their database tables
and leave them out of the database figures / statistics (the participant 5bis excluded in the chain)."""
import re
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

import minidb
import pipeline
from nbdriver import NotebookSession

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import participant_selection_lib as P     # noqa: E402


def tsv(path):
    return pd.read_csv(path, sep='\t', dtype={'file_id': str})


@pytest.fixture
def data_copy(chain, tmp_path):
    dst = tmp_path / 'data'
    shutil.copytree(chain['data'], dst)
    return dst


@pytest.mark.tool5
def test_tool5_knows_nothing_of_exclusions(chain):
    """Tool 5 extracts the raw macrostructure of every recording; the decisions belong to 5bis (its
    selection files + the registry). No exclusion column, every participant in the report's figures."""
    p = minidb.PARTICIPANTS[chain['dataset']]
    rep = chain['data'] / 'reports_features_macrostructure'
    g = tsv(rep / 'global_sleep_metrics.tsv')
    assert 'participant_excluded' not in g.columns and 'participant_exclude_reason' not in g.columns
    assert sorted(g['file_id']) == sorted(p)
    html = (rep / 'sleep_metrics_database_report.html').read_text(encoding='utf-8')
    assert 'not in the figures' not in html and 'EXCLUDED' not in html


@pytest.mark.tool7
def test_tool7_rerun_with_nothing_to_process_refreshes_its_tables(chain, data_copy):
    """A decision taken after tool 7 ran (here a tool-8 exclusion) reaches tool 7's database tables when
    tool 7 is run again, even with every participant skipped."""
    p = minidb.PARTICIPANTS[chain['dataset']]
    P.write_registry_rows(data_copy, '8_manual', [{'file_id': p[3], 'reason': 'bad night'}])
    stale = {}

    def check_stale_box(nb):                                       # shown at Load, before the run
        stale['before'] = nb.value('load_info.value')

    log = pipeline.run_tool7(data_copy, before_run=check_stale_box)   # everyone already processed
    assert 'registry changed after' in stale['before']
    assert 'No participant to process' in log
    g = tsv(data_copy / 'reports_preprocessing' / 'global_epoch_rejection.tsv')
    assert g.groupby('file_id')['participant_excluded'].first().to_dict() == {
        p[0]: False, p[1]: True, p[2]: False, p[3]: True}
    stage = tsv(data_copy / 'reports_preprocessing' / 'global_rejection_by_stage.tsv')
    assert (stage['n_participants'] == 2).all() and (stage['n_participants_excluded'] == 2).all()


@pytest.mark.tool6
def test_tool6_tables_flag_and_overview_leaves_out_the_excluded(chain):
    p = minidb.PARTICIPANTS[chain['dataset']]
    rep = chain['data'] / 'reports_quality_overview'
    for name in ('quality_summary.tsv', 'quality_summary_by_stage.tsv'):
        g = tsv(rep / name)
        assert g.groupby('file_id')['participant_excluded'].first().to_dict() == {
            p[0]: False, p[1]: True, p[2]: False, p[3]: False}
    html = (rep / 'dataset_overview.html').read_text(encoding='utf-8')
    assert re.search(r'<b>Participants:</b> 3\b', html)                   # 4 in the table, 3 in the stats
    assert '1 excluded participant(s) left out' in html and p[1] in html
    report = next(rep.rglob(f'{p[1]}_quality_overview.html')).read_text(encoding='utf-8')
    assert 'EXCLUDED participant' in report


@pytest.mark.tool9
def test_tool9_scan_warns_when_the_registry_is_newer_than_its_tables(chain, data_copy):
    p = minidb.PARTICIPANTS[chain['dataset']]
    raw, dec = data_copy / 'derivatives' / 'raw_epo', data_copy / 'derivatives' / 'rejection_auto'

    def scan_text():
        with NotebookSession(pipeline.T9) as nb:
            nb.pick('fc_data', data_copy)
            nb.pick('fc_raw', raw)
            nb.pick('fc_decision', dec)
            nb.click('btn_scan')
            return nb.widget_text()

    assert 'registry changed after' not in scan_text()             # tables built after the registry
    P.write_registry_rows(data_copy, '8_manual', [{'file_id': p[0], 'reason': 'bad night'}])
    assert 'registry changed after' in scan_text()
    pipeline.run_tool9(data_copy, raw, dec, subset=p, reaggregate=True)
    assert 'registry changed after' not in scan_text()             # refreshed by Re-aggregate only
    g = tsv(data_copy / 'reports_features_spectral' / 'global_spectral_stage.tsv')
    assert g.groupby('file_id')['participant_excluded'].first()[p[0]]
