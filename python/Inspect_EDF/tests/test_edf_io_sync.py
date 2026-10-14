"""The copies of edf_io_lib's reading functions that some tools still carry must stay identical to the library's.

The Curry generators of tools 2, 4, 6 and 7 rewrite these blocks, so they were not replaced by an import (tool-0
redesign, phase 8; moving them is a TODO, SPEC -> Planned modules), and quality_lib keeps its custom-stage helpers
because the tool-6 Curry generator copies that library. This test turns "keep in sync" into a check: it fails as
soon as one copy differs from tools/edf_io_lib.py (the code is compared, docstrings and comments are not) and names
the copy to update.

Copies that differ on purpose, so not checked: tool 1bis's read_edf_header_custom (it reads only the identity
fields it overwrites), tool 7's _events_df_from_xml (no Desaturation column: tool 5's, for the ODI) and tool 4's
event loaders (the event names only: another interface)."""
import ast
import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.quick

ROOT = Path(__file__).resolve().parent.parent
LIB = ROOT / 'tools' / 'edf_io_lib.py'

# tool file -> the names it still defines, each a copy of edf_io_lib's
COPIES = {
    'tools/1bis_anonymize_edf_voila.ipynb': ['detect_encoding'],
    'tools/2_select&remap_channels_edf_voila.ipynb': ['detect_encoding', 'read_edf_header_custom',
                                                      'extract_filter_value'],
    'tools/4_remap_events_edf_voila.ipynb': ['_TXT_DUR_RE', 'read_edf_start_datetime', 'event_companion_paths'],
    'tools/6_quality_overview_voila.ipynb': ['drop_suffix_duplicates', 'adapt_remap_dict_to_suffixes',
                                             'get_phys_bounds_uV'],
    'tools/7_preprocessing_voila.ipynb': ['drop_suffix_duplicates', 'adapt_remap_dict_to_suffixes', '_TXT_DUR_RE',
                                          'read_edf_start_datetime', 'event_companion_paths', '_events_df_from_txt',
                                          '_events_df_from_csv', '_events_df_from_manual', 'augment_remap_for_manual',
                                          'load_events', 'load_event_remap'],
    'tools/quality_lib.py': ['BASE_STAGE_COLORS', 'CUSTOM_STAGE_PALETTE', 'load_custom_stages', 'parse_custom_field',
                             'custom_stage_style'],
}


def definitions(path):
    """The top-level functions and assignments of a .py file, or of the code cells of a notebook: name -> AST node."""
    if path.suffix == '.ipynb':
        cells = json.loads(path.read_text(encoding='utf-8'))['cells']
        sources = [c['source'] if isinstance(c['source'], str) else ''.join(c['source'])
                   for c in cells if c['cell_type'] == 'code']
    else:
        sources = [path.read_text(encoding='utf-8')]
    found = {}
    for src in sources:
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                found[node.name] = node
            elif isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
                found[node.targets[0].id] = node
    return found


def code(node):
    """The code of a definition, its docstring left out (comments are not part of the syntax tree)."""
    if isinstance(node, ast.FunctionDef):
        body = node.body
        if (body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)):
            body = body[1:]
        return ast.dump(ast.Module(body=body, type_ignores=[])) + ast.dump(node.args)
    return ast.dump(node.value)


@pytest.mark.parametrize('tool', list(COPIES))
def test_copies_match_edf_io_lib(tool):
    lib = definitions(LIB)
    mine = definitions(ROOT / tool)
    for name in COPIES[tool]:
        assert name in mine, f'{tool} no longer defines {name}: remove it from COPIES in this test'
        assert code(mine[name]) == code(lib[name]), (
            f'{tool}: {name} differs from tools/edf_io_lib.py. Make both the same (then re-run the tool\'s Curry '
            'generator, tools_curry/_make_tool*_curry.py, when it has one).')
