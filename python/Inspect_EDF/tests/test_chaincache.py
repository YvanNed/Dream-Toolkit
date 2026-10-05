"""Unit tests of the chain reuse (chaincache.py): which tools rerun, and what is deleted before they do."""
import sys
from pathlib import Path

import pytest

import chaincache as C

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import participant_selection_lib as P     # noqa: E402

pytestmark = pytest.mark.quick
ORDER = ['5', '5bis', '6', '7', '8bis', '8', '9']


def test_first_tool_to_rerun():
    base, fp = 'B', {t: t + '-v1' for t in ORDER}
    state = {'base': 'B', 'tools': dict(fp)}
    assert C.first_tool_to_rerun(ORDER, state, base, fp) == len(ORDER)          # nothing changed
    assert C.first_tool_to_rerun(ORDER, state, base, dict(fp, **{'9': '9-v2'})) == 6   # tool 9 only
    assert C.first_tool_to_rerun(ORDER, state, base, dict(fp, **{'7': '7-v2'})) == 3   # 7 and after
    assert C.first_tool_to_rerun(ORDER, state, 'B2', fp) is None               # base changed: scratch
    assert C.first_tool_to_rerun(ORDER, None, base, fp) is None                # no state: scratch
    partial = {'base': 'B', 'tools': {t: fp[t] for t in ORDER[:4]}}            # interrupted after 7
    assert C.first_tool_to_rerun(ORDER, partial, base, fp) == 4


def test_a_shared_module_changes_every_tool_using_it():
    base, fp = C.fingerprints('synthetic')
    assert C.LIB_QC in C.TOOL_DEPS['8'] and C.LIB_QC in C.TOOL_DEPS['8bis']
    assert all(C.LIB_PSL in deps for deps in C.TOOL_DEPS.values())
    assert set(fp) == set(ORDER)
    assert C.fingerprints('real')[0] != base                                   # one state per dataset


def test_rewind_deletes_only_the_rerun_tools(tmp_path):
    data = tmp_path / 'data'
    for tool, rels in C.TOOL_OUTPUTS.items():
        for rel in rels:
            p = data / rel
            if '.' in p.name:
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(tool, encoding='utf-8')
            else:
                (p / 'sub').mkdir(parents=True, exist_ok=True)
                (p / 'sub' / 'f.tsv').write_text(tool, encoding='utf-8')
    (data / 'config_param').mkdir()
    P.write_registry_rows(data, '5bis_criteria', [{'file_id': 'a', 'reason': 'r'}])
    P.write_registry_rows(data, '8bis_auto', [{'file_id': 'b', 'reason': 'r'}])

    C.rewind(data, ['8bis', '8', '9'])
    for tool, rels in C.TOOL_OUTPUTS.items():
        gone = tool in ('8bis', '8', '9')
        for rel in rels:
            assert (data / rel).exists() != gone, (tool, rel)
    reg, _ = P.load_registry(data)
    assert reg['source'].tolist() == ['5bis_criteria']                         # 8bis rows only removed
