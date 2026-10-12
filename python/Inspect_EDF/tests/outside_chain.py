r"""Tools that are not part of the test chain (3bis today), run on a copy of the INPUTS of a chain dataset,
for the before / after check of the tool-0 redesign (code moved into tools/*_lib.py).

    & $py tests/outside_chain.py                 # build %TEMP%/dtk_outside_chain/<dataset>/data, run tool 3bis
    $env:DTK_TEST_ROOT = "$env:TEMP\dtk_outside_chain"; $env:DTK_REFACTOR_ROOT = "$env:TEMP\dtk_outside_saved"
    & $py tests/refactor_check.py save before    # then, after the move: run this script again + compare before

The inputs (recordings, hypnograms, event companions, config_param) come from the chain folders the test suite
builds (%TEMP%/dtk_tests/<dataset>/data): run the suite once first. Every run starts from a fresh copy, so the
outputs never mix two runs.
"""
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from nbdriver import NotebookSession  # noqa: E402

CHAIN_ROOT = Path(os.path.realpath(os.path.join(tempfile.gettempdir(), 'dtk_tests')))
ROOT = Path(os.path.realpath(os.path.join(tempfile.gettempdir(), 'dtk_outside_chain')))
T3BIS = 'tools/3bis_autoscoring_gssc_voila.ipynb'
INPUT_SUFFIXES = ('.edf', '.txt', '.csv', '.xml', '.tsv')


def build(dataset):
    """Fresh copy of the chain dataset's inputs in ROOT/<dataset>/data; returns that folder."""
    src = CHAIN_ROOT / dataset / 'data'
    data = ROOT / dataset / 'data'
    if data.exists():
        shutil.rmtree(data)
    data.mkdir(parents=True)
    for f in src.iterdir():
        if f.is_file() and f.suffix.lower() in INPUT_SUFFIXES:
            shutil.copy2(f, data / f.name)
    shutil.copytree(src / 'config_param', data / 'config_param')
    return data


def run_tool3bis(data):
    """Score every participant with the defaults (one electrode, C3, EOG L-R when declared)."""
    with NotebookSession(T3BIS) as nb:
        nb.pick('fc_data', data)
        nb.click('btn_scan')
        nb.click('btn_run')
        return nb.widget_text()


if __name__ == '__main__':
    for dataset in ('synthetic', 'real'):
        if not (CHAIN_ROOT / dataset / 'data').is_dir():
            print(f'[{dataset}] no chain folder: run the test suite first')
            continue
        t0 = time.time()
        text = run_tool3bis(build(dataset))
        bad = 'Traceback (most recent call last)' in text
        print(f'[{dataset}] tool 3bis done in {time.time() - t0:.0f} s' + ('  ⚠ TRACEBACK in its log' if bad else ''))
        if bad:
            print(text[-3000:])
