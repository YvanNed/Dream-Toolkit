"""Unit tests of tools/participant_selection_lib.py: the exclusion registry rules and the selector.
Fast (no notebook, no EDF): they run in a few seconds."""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import participant_selection_lib as P     # noqa: E402

pytestmark = pytest.mark.quick


@pytest.fixture
def data(tmp_path):
    (tmp_path / 'config_param').mkdir()
    return tmp_path


def _status(data):
    reg, warn = P.load_registry(data)
    assert warn is None
    return P.effective_status(reg)


# ---------------------------------------------------------------- registry ----
def test_absent_registry_means_everyone_included(data):
    reg, warn = P.load_registry(data)
    assert len(reg) == 0 and warn is None
    assert P.status_of(P.effective_status(reg), 'any') == (False, '')


def test_each_source_only_replaces_its_own_rows(data):
    P.write_registry_rows(data, '5bis_criteria', [{'file_id': 'a', 'reason': 'ahi = 30 > 15'},
                                                  {'file_id': 'b', 'reason': 'tst_min = 200 < 240'}])
    P.write_registry_rows(data, '8bis_auto', [{'file_id': 'b', 'reason': 'all in-scope epochs rejected'}],
                          replace_ids=['b'])
    # a new 5bis validation replaces ALL 5bis_criteria rows, and leaves the 8bis row alone
    P.write_registry_rows(data, '5bis_criteria', [{'file_id': 'c', 'reason': 'se_pct = 60 < 85'}])
    reg, _ = P.load_registry(data)
    assert sorted(zip(reg['file_id'], reg['source'])) == [('b', '8bis_auto'), ('c', '5bis_criteria')]
    st = P.effective_status(reg)
    assert P.status_of(st, 'a') == (False, '')
    assert P.status_of(st, 'b') == (True, '8bis_auto: all in-scope epochs rejected')
    assert P.excluded_ids(st) == ['b', 'c']


def test_replace_ids_limits_the_rewrite_to_the_processed_participants(data):
    P.write_registry_rows(data, '8bis_auto', [{'file_id': 'a', 'reason': 'x'}, {'file_id': 'b', 'reason': 'y'}])
    P.write_registry_rows(data, '8bis_auto', [], replace_ids=['a'])        # 'a' reprocessed, no longer excluded
    assert P.excluded_ids(_status(data)) == ['b']


def test_reasons_of_several_sources_are_joined(data):
    P.write_registry_rows(data, '5bis_criteria', [{'file_id': 'a', 'reason': 'ahi = 30 > 15'}])
    P.write_registry_rows(data, '7_manual', [{'file_id': 'a', 'reason': 'noisy night'}], replace_ids=['a'])
    excluded, reason = P.status_of(_status(data), 'a')
    assert excluded
    assert reason == '5bis_criteria: ahi = 30 > 15 | 7_manual: noisy night'


def test_forced_inclusion_overrides_every_exclusion(data):
    P.write_registry_rows(data, '5bis_criteria', [{'file_id': 'a', 'reason': 'ahi = 30 > 15'}])
    P.write_registry_rows(data, '8bis_auto', [{'file_id': 'a', 'reason': '312 arousal events > 200'}])
    P.write_registry_rows(data, '5bis_manual', [{'file_id': 'a', 'decision': 'include',
                                                 'comment': 'arousals checked by hand'}], replace_ids=['a'])
    st = _status(data)
    assert P.status_of(st, 'a') == (False, '')
    assert st[P.id_key('a')]['forced_include']
    assert '8bis_auto' in st[P.id_key('a')]['overridden_reason']


def test_only_the_manual_5bis_editor_may_force_an_inclusion(data):
    with pytest.raises(ValueError):
        P.write_registry_rows(data, '8bis_auto', [{'file_id': 'a', 'decision': 'include'}])
    with pytest.raises(ValueError):
        P.write_registry_rows(data, 'tool_x', [{'file_id': 'a'}])


def test_ids_match_case_insensitively_on_windows_like_paths(data):
    P.write_registry_rows(data, '7_manual', [{'file_id': 'Sub01', 'reason': 'r'}])
    excluded, _ = P.status_of(_status(data), 'Sub01')
    assert excluded


def test_unreadable_registry_warns_and_refuses_to_be_overwritten(data):
    P.registry_path(data).write_text('this is\tnot\na registry\n', encoding='utf-8')
    reg, warn = P.load_registry(data)
    assert len(reg) == 0 and 'everyone treated as included' in warn
    with pytest.raises(RuntimeError):
        P.write_registry_rows(data, '7_manual', [{'file_id': 'a', 'reason': 'r'}])
    assert 'not\na registry' in P.registry_path(data).read_text(encoding='utf-8')   # untouched


def test_no_temporary_file_left_behind(data):
    P.write_registry_rows(data, '7_manual', [{'file_id': 'a', 'reason': 'r'}])
    assert [f.name for f in (data / 'config_param').iterdir()] == [P.REGISTRY_NAME]


def test_add_exclusion_columns_reflects_the_current_registry(data):
    P.write_registry_rows(data, '7_manual', [{'file_id': 'b', 'reason': 'r'}])
    df = pd.DataFrame({'file_id': ['a', 'b'], 'x': [1, 2],
                       P.COL_EXCLUDED: [True, False], P.COL_REASON: ['stale', '']})   # stale columns
    out = P.add_exclusion_columns(df, _status(data))
    assert list(out.columns) == ['file_id', 'x', P.COL_EXCLUDED, P.COL_REASON]
    assert out[P.COL_EXCLUDED].tolist() == [False, True]
    assert out[P.COL_REASON].tolist() == ['', '7_manual: r']


def test_find_data_folder(data):
    deep = data / 'derivatives' / 'raw_epo' / 'group1'
    deep.mkdir(parents=True)
    assert P.find_data_folder(deep) == data                               # via config_param/
    other = data.parent / (data.name + '_noconfig') / 'reports_preprocessing' / 'g'
    other.mkdir(parents=True)
    assert P.find_data_folder(other) == other.parent.parent               # via reports_*


# ---------------------------------------------------------------- selector ----
def make_selector():
    sel = P.ParticipantSelector()
    sel.set_participants(['p1', 'p2', 'p3', 'p4'], done={'p1', 'p2'}, excluded={'p3': '5bis: ahi'})
    return sel


def test_all_mode_honours_skip():
    sel = make_selector()
    assert sel.selected_ids() == ['p3', 'p4']
    sel.cb_skip.value = False
    assert sel.selected_ids() == ['p1', 'p2', 'p3', 'p4']
    assert sel.warning_html() == ''
    assert 'excluded' in sel.lbl_summary.value


def test_subset_with_skip_warns_about_the_processed_ones():
    sel = make_selector()
    sel.dd_mode.value = 'Subset only'
    sel.tags.value = ['p2', 'p4']
    assert sel.selected_ids() == ['p4']
    assert sel.skipped_ids() == ['p2']
    assert 'will NOT be processed' in sel.html_warn.value and 'p2' in sel.html_warn.value
    sel.cb_skip.value = False                                  # untick -> reprocessed, warning gone
    assert sel.selected_ids() == ['p2', 'p4']
    assert sel.html_warn.value == ''


def test_paste_and_list_pick_switch_to_subset_mode():
    sel = make_selector()
    sel.txt_paste.value = 'p4, p3; nope'
    sel._on_paste()
    assert sel.subset_mode and sel.tags.value == ['p4', 'p3']
    assert 'nope' in sel.lbl_paste.value
    sel.lst.value = 'p1'                                         # click in the participant list
    assert sel.tags.value == ['p4', 'p3', 'p1'] and sel.lst.value is None
    assert sel.chosen_ids() == ['p1', 'p3', 'p4']                # display order, not click order


def test_rescan_drops_vanished_ids_from_the_subset():
    sel = make_selector()
    sel.dd_mode.value = 'Subset only'
    sel.tags.value = ['p1', 'p4']
    sel.set_participants(['p1', 'p2'], done=set())
    assert sel.tags.value == ['p1']


def test_add_box_offers_only_the_participants_not_selected_yet():
    sel = make_selector()
    assert list(sel.cmb_add.options) == ['p1', 'p2', 'p3', 'p4']
    sel.cmb_add.value = 'p2'                                     # picked in the autocompletion + Enter
    assert sel.subset_mode and sel.tags.value == ['p2'] and sel.cmb_add.value == ''
    assert list(sel.cmb_add.options) == ['p1', 'p3', 'p4']
    sel.tags.value = ['p2', 'zzz', 'p2']                         # typed straight into the tag box
    assert sel.tags.value == ['p2'] and 'zzz' in sel.lbl_paste.value
    sel.tags.value = []                                          # 'x' on the tag: offered again
    assert list(sel.cmb_add.options) == ['p1', 'p2', 'p3', 'p4']


def test_empty_registry_header_only_means_everyone_included(data):
    """A registry whose rows were all removed (e.g. a 5bis validation without any criterion) keeps its
    header only: it must read as 'nobody excluded', and the exclusion columns still be added."""
    P.write_registry_rows(data, '5bis_criteria', [{'file_id': 'a', 'reason': 'r'}])
    P.write_registry_rows(data, '5bis_criteria', [])                # new validation, no criterion
    assert P.registry_path(data).read_text(encoding='utf-8').strip() == '\t'.join(P.REGISTRY_COLUMNS)
    reg, warn = P.load_registry(data)
    assert warn is None and len(reg) == 0
    out = P.add_exclusion_columns(pd.DataFrame({'file_id': ['a', 'b']}), P.effective_status(reg))
    assert out[P.COL_EXCLUDED].tolist() == [False, False] and out[P.COL_REASON].tolist() == ['', '']
    P.registry_path(data).write_text('', encoding='utf-8')         # an empty FILE (0 byte) too:
    reg, warn = P.load_registry(data)                              # no decision, no warning, and
    assert len(reg) == 0 and warn is None                          # a tool can still write into it
    P.write_registry_rows(data, '7_manual', [{'file_id': 'a', 'reason': 'r'}])
    assert P.excluded_ids(P.effective_status(P.load_registry(data)[0])) == ['a']


def test_registry_newer_than_its_tables(data, tmp_path):
    import os, time
    table = tmp_path / 'global.tsv'
    assert not P.registry_newer_than(data, table)                  # nothing to compare yet
    table.write_text('x\n', encoding='utf-8')
    assert not P.registry_newer_than(data, table)                  # no registry
    P.write_registry_rows(data, '7_manual', [{'file_id': 'a', 'reason': 'r'}])
    later = time.time() + 5
    os.utime(P.registry_path(data), (later, later))
    assert P.registry_newer_than(data, table)                      # decision taken after the table
    os.utime(table, (later + 5, later + 5))
    assert not P.registry_newer_than(data, table)                  # table rebuilt since
    assert 'registry changed after' in P.stale_tables_html(P.STALE_RERUN)


def test_rewriting_the_same_decisions_leaves_the_registry_untouched(data):
    """Tool 8bis re-run with the same thresholds: no new decision, so no new file date (the stale-table
    warning of the other tools must not fire) and the original decided_at is kept."""
    import os
    P.write_registry_rows(data, '8bis_auto', [{'file_id': 'a', 'reason': 'all channels dropped'}],
                          replace_ids=['a'])
    P.write_registry_rows(data, '8bis_auto', [], replace_ids=['b'])      # nothing for b, before or after
    path = P.registry_path(data)
    os.utime(path, (1_000_000, 1_000_000))
    before = path.read_text(encoding='utf-8')
    P.write_registry_rows(data, '8bis_auto', [{'file_id': 'a', 'reason': 'all channels dropped'}],
                          replace_ids=['a'])
    P.write_registry_rows(data, '8bis_auto', [], replace_ids=['b'])
    assert path.stat().st_mtime == 1_000_000 and path.read_text(encoding='utf-8') == before
    P.write_registry_rows(data, '8bis_auto', [{'file_id': 'a', 'reason': 'all in-scope epochs rejected'}],
                          replace_ids=['a'])                               # a real change: written
    assert path.stat().st_mtime != 1_000_000
