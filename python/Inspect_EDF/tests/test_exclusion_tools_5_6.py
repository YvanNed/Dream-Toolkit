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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import participant_selection_lib as P     # noqa: E402


def tsv(path):
    return pd.read_csv(path, sep='\t', dtype={'file_id': str})


@pytest.fixture
def data_copy(chain, tmp_path):
    dst = tmp_path / 'data'
    shutil.copytree(chain['data'], dst)
    return dst


def test_tool5_picks_up_the_registry_at_its_next_run(chain, data_copy):
    """In the chain tool 5 runs BEFORE 5bis decides (5bis needs tool 5's table), so its table carries
    no exclusion yet. Running tool 5 again (every recording skipped: nothing is recomputed) rebuilds
    the database table with the current registry."""
    p = minidb.PARTICIPANTS[chain['dataset']]
    rep = data_copy / 'reports_features_macrostructure'
    assert not tsv(rep / 'global_sleep_metrics.tsv')['participant_excluded'].any()
    before = {f.name: f.stat().st_mtime_ns for f in (data_copy / 'derivatives').rglob('*_sleep_metrics.tsv')}
    pipeline.run_tool5(data_copy)
    after = {f.name: f.stat().st_mtime_ns for f in (data_copy / 'derivatives').rglob('*_sleep_metrics.tsv')}
    assert before == after                                         # nothing recomputed
    g = tsv(rep / 'global_sleep_metrics.tsv')
    assert g.set_index('file_id')['participant_excluded'].to_dict() == {
        p[0]: False, p[1]: True, p[2]: False, p[3]: False}
    assert g.set_index('file_id').loc[p[1], 'participant_exclude_reason'] == '5bis_criteria: age = 70 > 60'
    html = (rep / 'sleep_metrics_database_report.html').read_text(encoding='utf-8')
    assert p[1] in html and 'not in the figures' in html


def test_tool7_rerun_with_nothing_to_process_refreshes_its_tables(chain, data_copy):
    """A decision taken after tool 7 ran (here a tool-8 exclusion) reaches tool 7's database tables when
    tool 7 is run again, even with every participant skipped."""
    p = minidb.PARTICIPANTS[chain['dataset']]
    P.write_registry_rows(data_copy, '8_manual', [{'file_id': p[3], 'reason': 'bad night'}])
    log = pipeline.run_tool7(data_copy)                            # every participant already processed
    assert 'No participant to process' in log
    g = tsv(data_copy / 'reports_preprocessing' / 'global_epoch_rejection.tsv')
    assert g.groupby('file_id')['participant_excluded'].first().to_dict() == {
        p[0]: False, p[1]: True, p[2]: False, p[3]: True}
    stage = tsv(data_copy / 'reports_preprocessing' / 'global_rejection_by_stage.tsv')
    assert (stage['n_participants'] == 2).all() and (stage['n_participants_excluded'] == 2).all()


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
