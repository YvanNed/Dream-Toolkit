r"""Shared fixtures. Run from the repo root (the three levels, see TESTS.md):
    & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" -m pytest -m quick         # ~2 min
    & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" -m pytest -m "not full"    # standard, ~9 min (~5 with the chain reused)
    & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" -m pytest                  # full, ~24 min
Only the tests of some tools: -m "tool9 and not full" (markers tool1 ... tool9, see pytest.ini).
--fresh-chain forces the chain to run from scratch instead of reusing its last run (chaincache.py).

The chain runs once per dataset and per session (fixture `chain`), reusing the last run for the tools whose
code did not change; every test reads its outputs. Run folders are kept under %TEMP%/dtk_tests/<dataset>/.
The 'real' dataset (marker `full`) is skipped on a machine without tools/test_data.
"""
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import chaincache  # noqa: E402
import minidb      # noqa: E402
import pipeline    # noqa: E402
import snapshot    # noqa: E402

GOLDEN = HERE / 'golden'
DATASETS = ['synthetic', 'real']
# A short fixed root: the toolkit writes deep output trees and Windows caps paths at 260 characters.
RUN_ROOT = Path(os.path.realpath(os.environ.get('DTK_TEST_ROOT',
                                                os.path.join(tempfile.gettempdir(), 'dtk_tests'))))


def pytest_addoption(parser):
    parser.addoption('--fresh-chain', action='store_true', default=False,
                     help='run the test chain from scratch instead of reusing its last run')


def run_dataset_chain(dataset, root, log=print):
    """Build `dataset` in `root` from scratch, run the whole chain, snapshot the outputs (make_golden.py)."""
    if root.exists():
        shutil.rmtree(root)
    t0 = time.time()
    data = minidb.build(root, dataset)
    logs = pipeline.run_chain(data, minidb.manual_review_participant(dataset),
                              log=lambda m: log(f'[{dataset} {time.time() - t0:5.0f}s] {m}'))
    snap = snapshot.take(data, root / 'snapshot')
    return {'dataset': dataset, 'root': root, 'data': data, 'snapshot': snap, 'logs': logs}


@pytest.fixture(scope='session', params=[pytest.param('synthetic'),
                                         pytest.param('real', marks=pytest.mark.full)])
def chain(request):
    dataset = request.param
    if not minidb.available(dataset):
        pytest.skip(f'{dataset} test data not available on this machine')
    root = RUN_ROOT / dataset
    manual = minidb.manual_review_participant(dataset)
    data, logs = chaincache.run_cached_chain(
        dataset, root, build=lambda r: minidb.build(r, dataset),
        steps_for=lambda d: pipeline.chain_steps(d, manual),
        fresh=request.config.getoption('--fresh-chain'), log=lambda m: print(m, flush=True))
    snap = snapshot.take(data, root / 'snapshot')
    return {'dataset': dataset, 'root': root, 'data': data, 'snapshot': snap, 'logs': logs}


@pytest.fixture(scope='session')
def golden(chain):
    folder = GOLDEN / chain['dataset']
    if not folder.is_dir():
        pytest.skip(f'no golden snapshot for {chain["dataset"]}: run tests/make_golden.py on trusted code')
    return folder
