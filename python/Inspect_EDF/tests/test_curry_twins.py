"""The Curry twins are generated from the EDF notebooks (tools_curry/_make_tool*_curry.py). They are not
run here (no Curry test data in the suite), but a twin left behind after an EDF change is caught: every
code cell must parse, and the features of the EDF original must be present in its twin."""
import ast
import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.quick

ROOT = Path(__file__).resolve().parent.parent

# twin -> markers that must appear in it (the EDF original carries them too)
TWINS = {
    'tools_curry/6_quality_overview_curry_voila.ipynb': [
        'participant_selector', 'add_exclusion_columns', 'excluded_banner_html', 'read_curry_header'],
    'tools_curry/7_preprocessing_curry_voila.ipynb': [
        'participant_selector', "raw.info['bads'] = bad_present", "'bad_channels'",
        "'7_manual'", 'n_participants_excluded', 'read_raw_curry'],
    'tools_curry/8_reject_manually_curry_voila.ipynb': [
        'participant_selection_lib', "'8_manual'", '_epoch_decision.tsv', '_channel_decision.tsv'],
    'tools_curry/8bis_reject_automatically_curry_voila.ipynb': [
        'participant_selector', "'8bis_auto'", '_epoch_decision.tsv', 'tool7_deselected'],
}
ORIGINALS = {
    'tools_curry/6_quality_overview_curry_voila.ipynb': 'tools/6_quality_overview_voila.ipynb',
    'tools_curry/7_preprocessing_curry_voila.ipynb': 'tools/7_preprocessing_voila.ipynb',
    'tools_curry/8_reject_manually_curry_voila.ipynb': 'tools/8_reject_manually_voila.ipynb',
    'tools_curry/8bis_reject_automatically_curry_voila.ipynb': 'tools/8bis_reject_automatically_voila.ipynb',
}


def code_of(path):
    nb = json.loads((ROOT / path).read_text(encoding='utf-8'))
    return [''.join(c['source']) if isinstance(c['source'], list) else c['source']
            for c in nb['cells'] if c['cell_type'] == 'code']


@pytest.mark.parametrize('twin', list(TWINS))
def test_twin_parses_and_follows_its_edf_original(twin):
    cells = code_of(twin)
    for i, src in enumerate(cells):
        ast.parse(src, filename=f'{twin} cell {i}')
    text = '\n'.join(cells)
    original = '\n'.join(code_of(ORIGINALS[twin]))
    for marker in TWINS[twin]:
        if marker not in ('read_raw_curry', 'read_curry_header'):
            assert marker in original, f'marker {marker!r} not in the EDF original: update this test'
        assert marker in text, f'{twin} lacks {marker!r}: re-run its tools_curry/_make_tool*_curry.py'
