r"""Shared fixtures. Run from the repo root:
    & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" -m pytest
    ... -m pytest -k synthetic        # the synthetic dataset only
    ... -m pytest -k real             # the real local dataset only

The full chain runs once per dataset and per session (fixture `chain`); every test reads its outputs.
Run folders are kept under %TEMP%/dtk_tests/<dataset>/ for inspection after a failure. The 'real'
dataset is skipped on a machine without tools/test_data.
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

import minidb      # noqa: E402
import pipeline    # noqa: E402
import snapshot    # noqa: E402

GOLDEN = HERE / 'golden'
DATASETS = ['synthetic', 'real']
# A short fixed root: the toolkit writes deep output trees and Windows caps paths at 260 characters.
RUN_ROOT = Path(os.path.realpath(os.environ.get('DTK_TEST_ROOT',
                                                os.path.join(tempfile.gettempdir(), 'dtk_tests'))))


def run_dataset_chain(dataset, root, log=print):
    """Build `dataset` in `root`, run 5 -> 5bis -> 6 -> 7 -> 8bis -> 8 -> 9, snapshot the outputs."""
    if root.exists():
        shutil.rmtree(root)
    t0 = time.time()
    data = minidb.build(root, dataset)
    logs = pipeline.run_chain(data, minidb.manual_review_participant(dataset),
                              log=lambda m: log(f'[{dataset} {time.time() - t0:5.0f}s] {m}'))
    snap = snapshot.take(data, root / 'snapshot')
    return {'dataset': dataset, 'root': root, 'data': data, 'snapshot': snap, 'logs': logs}


@pytest.fixture(scope='session', params=DATASETS)
def chain(request):
    dataset = request.param
    if not minidb.available(dataset):
        pytest.skip(f'{dataset} test data not available on this machine')
    return run_dataset_chain(dataset, RUN_ROOT / dataset)


@pytest.fixture(scope='session')
def golden(chain):
    folder = GOLDEN / chain['dataset']
    if not folder.is_dir():
        pytest.skip(f'no golden snapshot for {chain["dataset"]}: run tests/make_golden.py on trusted code')
    return folder
