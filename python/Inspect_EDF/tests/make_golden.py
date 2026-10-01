r"""(Re)build a golden snapshot: run the full chain on a test dataset with the CURRENT code and store
the outputs in tests/golden/<dataset>/.

Run it only on code you trust (before a change), from the repo root:
    & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" tests/make_golden.py --dataset synthetic
    & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" tests/make_golden.py --dataset real

tests/golden/synthetic/ is versioned; tests/golden/real/ is git-ignored (real participant data).
The run folder (%TEMP%/dtk_golden/<dataset>) is kept so the outputs can be inspected.
"""
import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import conftest    # noqa: E402
import minidb      # noqa: E402
import snapshot    # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', choices=conftest.DATASETS, required=True)
    args = ap.parse_args()
    if not minidb.available(args.dataset):
        sys.exit(f'{args.dataset} test data not available on this machine')
    root = Path(os.path.realpath(os.path.join(tempfile.gettempdir(), 'dtk_golden', args.dataset)))
    t0 = time.time()
    run = conftest.run_dataset_chain(args.dataset, root, log=lambda m: print(m, flush=True))
    (root / 'logs').mkdir(exist_ok=True)
    for name, text in run['logs'].items():
        (root / 'logs' / f'tool{name}.txt').write_text(text, encoding='utf-8')
    dest = snapshot.take(run['data'], HERE / 'golden' / args.dataset)
    n = sum(1 for f in dest.rglob('*') if f.is_file())
    print(f'golden snapshot: {n} files in {dest} ({time.time() - t0:.0f}s, run folder {root})')


if __name__ == '__main__':
    main()
