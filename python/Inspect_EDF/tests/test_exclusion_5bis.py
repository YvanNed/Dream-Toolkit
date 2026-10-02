"""Tool 5bis and the exclusion registry (step 2 of the exclusion refactor), driven through its widgets.

- Validation writes the participants the criteria exclude to config_param/participant_exclusions.tsv.
- Section 3 (manual decisions): exclusion, forced inclusion (comment required), deletion; a new
  validation replaces the criteria rows only.
- The exclusion columns other tools will add to global_sleep_metrics.tsv are never taken for metrics.
"""
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


def _rows(data):
    reg, warn = P.load_registry(data)
    assert warn is None
    return sorted(zip(reg['file_id'], reg['source'], reg['decision']))


@pytest.fixture
def macro_copy(chain, tmp_path):
    """A private copy of the chain's tool-5 outputs (+ config_param): these tests write decisions and
    must not touch the chain the other tests compare to the golden snapshot."""
    data = tmp_path / 'data'
    for sub in ('reports_features_macrostructure', 'config_param'):
        shutil.copytree(chain['data'] / sub, data / sub)
    return data


def test_validation_writes_the_criteria_exclusions_to_the_registry(chain):
    p = minidb.PARTICIPANTS[chain['dataset']]
    reg, warn = P.load_registry(chain['data'])
    assert warn is None
    assert list(zip(reg['file_id'], reg['source'], reg['decision'], reg['reason'])) == \
        [(p[1], '5bis_criteria', 'exclude', 'age = 70 > 60')]
    assert reg['comment'].tolist() == ['test selection']


def test_manual_decisions_combine_with_the_criteria(chain, macro_copy):
    p = minidb.PARTICIPANTS[chain['dataset']]
    with NotebookSession(pipeline.T5BIS) as nb:
        nb.pick('fc_data', macro_copy)
        nb.click('btn_load')
        assert '<b>3</b> kept / 4' in nb.value('lbl_summary.value')      # age > 60 pre-filled -> p1 out

        # a forced inclusion without a comment is refused, nothing written
        nb.set('dd_man_part', p[1])
        nb.run('tb_man_decision.value = PSL.INCLUDE')
        nb.click('btn_man_save')
        assert 'needs a comment' in nb.value('lbl_man.value')
        assert _rows(macro_copy) == [(p[1], '5bis_criteria', 'exclude')]

        nb.set('txt_man_comment', 'age checked: within the protocol limits')
        nb.click('btn_man_save')
        nb.set('dd_man_part', p[2])
        nb.run('tb_man_decision.value = PSL.EXCLUDE')
        nb.set('txt_man_reason', 'electrode off from 2 a.m.')
        nb.click('btn_man_save')
        summary = nb.value('lbl_summary.value')
        assert '<b>3</b> kept / 4' in summary and 'forced inclusions: 1' in summary   # p1 back in, p2 out
        assert p[2] in nb.value('html_excluded.value')
        n_lines = nb.value('len(box_man_rows.children)')
        assert n_lines == 3                                       # criteria row + 2 manual rows listed

        nb.click('btn_validate')                                  # criteria rows replaced, manual kept
        assert 'Exclusion registry updated' in nb.value('lbl_validate.value')

    assert _rows(macro_copy) == sorted([(p[1], '5bis_criteria', 'exclude'), (p[1], '5bis_manual', 'include'),
                                        (p[2], '5bis_manual', 'exclude')])
    st = P.effective_status(P.load_registry(macro_copy)[0])
    assert P.status_of(st, p[1]) == (False, '')
    assert P.status_of(st, p[2]) == (True, '5bis_manual: electrode off from 2 a.m.')
    # the selection table itself still records the criteria decision only
    sel = pd.read_csv(macro_copy / 'reports_features_macrostructure' / 'participant_selection.tsv', sep='\t')
    assert sel.loc[sel['file_id'] == p[1], 'excluded'].item()

    with NotebookSession(pipeline.T5BIS) as nb:                  # delete a manual decision
        nb.pick('fc_data', macro_copy)
        nb.click('btn_load')
        nb.run(f'on_man_delete({p[2]!r})')
        assert 'deleted' in nb.value('lbl_man.value')
        nb.run('for r in ROWS.values(): r["cb"].value = False')  # no criterion any more
        nb.click('btn_validate')
    assert _rows(macro_copy) == [(p[1], '5bis_manual', 'include')]


def test_exclusion_columns_are_not_metrics(chain, macro_copy):
    table = macro_copy / 'reports_features_macrostructure' / 'global_sleep_metrics.tsv'
    df = pd.read_csv(table, sep='\t', dtype={'file_id': str})
    df[P.COL_EXCLUDED] = [False, True, False, False]
    df[P.COL_REASON] = ['', 'x', '', '']
    df.to_csv(table, sep='\t', index=False)
    with NotebookSession(pipeline.T5BIS) as nb:
        nb.pick('fc_data', macro_copy)
        nb.click('btn_load')
        keys = nb.value("[m[0] for m in S['metrics']]")
    assert P.COL_EXCLUDED not in keys and P.COL_REASON not in keys
    assert 'age' in keys
