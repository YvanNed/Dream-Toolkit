r"""Before / after check for a refactor that must not change any output (tool 0 redesign: code moved into
shared libraries, see tools/plan_tool0_redesign.md).

The golden snapshot (snapshot.py) compares the machine-readable tables and summarises each .fif. Moving code
must leave EVERY output unchanged, so this check compares the whole output tree of the test chain, file by
file, the human outputs included:

    .tsv / .csv          every value (strings exactly, numbers exactly, else within the golden tolerance)
    .json                every value (same rule)
    -epo.fif             signal array, metadata table, channel names / types, bads, sfreq, tmin, meas_date
    .npz                 every array
    .xlsx                every sheet, as tables
    .html / .txt / .md   text, once the clock times written in it are masked (see MASKS)
    anything else        bytes (.png ...)

The recordings themselves (.edf ...) are inputs that no tool modifies: they are neither copied nor compared.

Usage, from the repo root:
    $py = "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe"
    & $py tests/refactor_check.py save before        # keep a copy of the current chain outputs
    ... move the code, run the tests: the chain reruns from the first changed tool ...
    & $py tests/refactor_check.py compare before     # every difference; exit code 1 when any

The chain outputs are read in %TEMP%/dtk_tests/<dataset>/data (the folder conftest.py uses), the copies are
kept in %TEMP%/dtk_refactor/<label>/<dataset>/. A tool the chain did not rerun has the very same files, so
only the rerun tools are really compared: rerun the ones the change touches.
"""
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from snapshot import VOLATILE  # noqa: E402  (columns / keys that legitimately differ between two runs)

DATASETS = ['synthetic', 'real']
RUN_ROOT = Path(os.path.realpath(os.environ.get('DTK_TEST_ROOT',
                                                os.path.join(tempfile.gettempdir(), 'dtk_tests'))))
SAVE_ROOT = Path(os.path.realpath(os.environ.get('DTK_REFACTOR_ROOT',
                                                 os.path.join(tempfile.gettempdir(), 'dtk_refactor'))))

# Recordings: inputs, never written by a tool.
RECORDING_SUFFIXES = {'.edf', '.bdf', '.cdt', '.dpo', '.dap', '.rs3', '.ceo'}
# Same tolerance as the golden comparison (snapshot.compare_tables).
RTOL, ATOL = 1e-6, 1e-9
# Clock times written into the human outputs (report dates, run durations). Each mask is a deliberate,
# listed exception: anything else that differs between two runs is reported.
MASKS = [
    (re.compile(r'\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}(:\d{2}(\.\d+)?)?'), '<DATETIME>'),
]


# ---------------------------------------------------------------------------------------------------------
# Copy
# ---------------------------------------------------------------------------------------------------------
def output_files(data):
    """Every file of the chain folder except the recordings, as relative paths."""
    data = Path(data)
    return sorted(f.relative_to(data) for f in data.rglob('*')
                  if f.is_file() and f.suffix.lower() not in RECORDING_SUFFIXES)


def save(label):
    for dataset in DATASETS:
        data = RUN_ROOT / dataset / 'data'
        if not data.is_dir():
            print(f'[{dataset}] no chain outputs in {data}: skipped')
            continue
        dest = SAVE_ROOT / label / dataset
        if dest.exists():
            shutil.rmtree(dest)
        files = output_files(data)
        for rel in files:
            (dest / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(data / rel, dest / rel)
        print(f'[{dataset}] {len(files)} files saved to {dest}')


# ---------------------------------------------------------------------------------------------------------
# Compare: each function returns (status, details) with status 'same', 'close' (numbers within the
# tolerance only) or 'different'.
# ---------------------------------------------------------------------------------------------------------
def _compare_arrays(old, new, where):
    old, new = np.asarray(old), np.asarray(new)
    if old.shape != new.shape:
        return 'different', [f'{where}: shape {old.shape} -> {new.shape}']
    if old.dtype.kind in 'fc' or new.dtype.kind in 'fc':
        if np.array_equal(old, new, equal_nan=True):
            return 'same', []
        close = np.isclose(old, new, rtol=RTOL, atol=ATOL, equal_nan=True)
        diff = np.abs(old.astype(float) - new.astype(float))
        detail = f'{where}: {int((~np.isclose(old, new, rtol=0, atol=0, equal_nan=True)).sum())} value(s) ' \
                 f'differ, max |diff| = {np.nanmax(diff):.3g}'
        return ('close' if close.all() else 'different'), [detail]
    if np.array_equal(old, new):
        return 'same', []
    return 'different', [f'{where}: {int((old != new).sum())} value(s) differ']


def _compare_frames(old, new, where):
    if list(old.columns) != list(new.columns):
        return 'different', [f'{where}: columns {list(old.columns)} -> {list(new.columns)}']
    if len(old) != len(new):
        return 'different', [f'{where}: row count {len(old)} -> {len(new)}']
    status, details = 'same', []
    for c in old.columns:
        if c in VOLATILE:
            continue
        a, b = old[c], new[c]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            s, d = _compare_arrays(a.to_numpy(dtype=float), b.to_numpy(dtype=float), f'{where} column {c!r}')
        else:
            av = a.fillna('<NA>').astype(str).to_numpy()
            bv = b.fillna('<NA>').astype(str).to_numpy()
            bad = av != bv
            if not bad.any():
                continue
            i = int(np.flatnonzero(bad)[0])
            s, d = 'different', [f'{where} column {c!r}: {int(bad.sum())} value(s) differ, '
                                 f'first at row {i}: {av[i]!r} -> {bv[i]!r}']
        if s == 'different':
            status = 'different'
        elif s == 'close' and status == 'same':
            status = 'close'
        details += d
    return status, details


def _compare_json(old, new, where=''):
    if isinstance(old, dict) and isinstance(new, dict):
        status, details = 'same', []
        if set(old) - VOLATILE != set(new) - VOLATILE:
            return 'different', [f'{where}: keys {sorted(set(old) - VOLATILE)} -> {sorted(set(new) - VOLATILE)}']
        for k in old:
            if k in VOLATILE:
                continue
            s, d = _compare_json(old[k], new[k], f'{where}/{k}')
            status = 'different' if 'different' in (status, s) else ('close' if 'close' in (status, s) else 'same')
            details += d
        return status, details
    if isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            return 'different', [f'{where}: length {len(old)} -> {len(new)}']
        status, details = 'same', []
        for i, (a, b) in enumerate(zip(old, new)):
            s, d = _compare_json(a, b, f'{where}[{i}]')
            status = 'different' if 'different' in (status, s) else ('close' if 'close' in (status, s) else 'same')
            details += d
        return status, details
    if isinstance(old, float) and isinstance(new, (int, float)) and not isinstance(new, bool):
        if old == new or (np.isnan(old) and np.isnan(new)):
            return 'same', []
        if np.isclose(old, new, rtol=RTOL, atol=ATOL):
            return 'close', [f'{where}: {old!r} -> {new!r}']
    if old == new:
        return 'same', []
    return 'different', [f'{where}: {old!r} -> {new!r}']


def _read_table(path):
    sep = '\t' if path.suffix.lower() == '.tsv' else ','
    return pd.read_csv(path, sep=sep, keep_default_na=True)


def _compare_epochs(old_path, new_path):
    import mne
    old = mne.read_epochs(str(old_path), preload=True, verbose=False)
    new = mne.read_epochs(str(new_path), preload=True, verbose=False)
    checks = [('ch_names', list(old.ch_names), list(new.ch_names)),
              ('channel types', old.get_channel_types(), new.get_channel_types()),
              ('bads', list(old.info['bads']), list(new.info['bads'])),
              ('sfreq', float(old.info['sfreq']), float(new.info['sfreq'])),
              ('tmin', float(old.tmin), float(new.tmin)),
              ('meas_date', str(old.info['meas_date']), str(new.info['meas_date']))]
    for name, a, b in checks:
        if a != b:
            return 'different', [f'{name}: {a!r} -> {b!r}']
    status, details = _compare_arrays(old.get_data(), new.get_data(), 'signal')
    if (old.metadata is None) != (new.metadata is None):
        return 'different', details + ['metadata present in one file only']
    if old.metadata is not None:
        s, d = _compare_frames(old.metadata.reset_index(drop=True), new.metadata.reset_index(drop=True),
                               'metadata')
        status = 'different' if 'different' in (status, s) else ('close' if 'close' in (status, s) else 'same')
        details += d
    return status, details


def _compare_npz(old_path, new_path):
    with np.load(old_path, allow_pickle=False) as old, np.load(new_path, allow_pickle=False) as new:
        if sorted(old.files) != sorted(new.files):
            return 'different', [f'arrays {sorted(old.files)} -> {sorted(new.files)}']
        status, details = 'same', []
        for k in old.files:
            s, d = _compare_arrays(old[k], new[k], k)
            status = 'different' if 'different' in (status, s) else ('close' if 'close' in (status, s) else 'same')
            details += d
        return status, details


def _compare_xlsx(old_path, new_path):
    old = pd.read_excel(old_path, sheet_name=None)
    new = pd.read_excel(new_path, sheet_name=None)
    if list(old) != list(new):
        return 'different', [f'sheets {list(old)} -> {list(new)}']
    status, details = 'same', []
    for sheet in old:
        s, d = _compare_frames(old[sheet], new[sheet], f'sheet {sheet!r}')
        status = 'different' if 'different' in (status, s) else ('close' if 'close' in (status, s) else 'same')
        details += d
    return status, details


def _masked(text):
    for pattern, repl in MASKS:
        text = pattern.sub(repl, text)
    return text


def _compare_text(old_path, new_path):
    a = _masked(old_path.read_text(encoding='utf-8', errors='replace'))
    b = _masked(new_path.read_text(encoding='utf-8', errors='replace'))
    if a == b:
        return 'same', []
    la, lb = a.splitlines(), b.splitlines()
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            return 'different', [f'line {i + 1}: {x[:150]!r} -> {y[:150]!r}']
    return 'different', [f'line count {len(la)} -> {len(lb)}']


def compare_file(old_path, new_path):
    name = new_path.name.lower()
    suffix = new_path.suffix.lower()
    if old_path.read_bytes() == new_path.read_bytes():
        return 'same', []
    if name.endswith('-epo.fif'):
        return _compare_epochs(old_path, new_path)
    if suffix in ('.tsv', '.csv'):
        return _compare_frames(_read_table(old_path), _read_table(new_path), 'table')
    if suffix == '.json':
        return _compare_json(json.loads(old_path.read_text(encoding='utf-8')),
                             json.loads(new_path.read_text(encoding='utf-8')))
    if suffix == '.npz':
        return _compare_npz(old_path, new_path)
    if suffix == '.xlsx':
        return _compare_xlsx(old_path, new_path)
    if suffix in ('.html', '.txt', '.md'):
        return _compare_text(old_path, new_path)
    return 'different', ['bytes differ']


def compare(label):
    """Print the report; return the number of files that differ (or are missing / added)."""
    n_bad = 0
    for dataset in DATASETS:
        saved = SAVE_ROOT / label / dataset
        data = RUN_ROOT / dataset / 'data'
        if not saved.is_dir():
            print(f'[{dataset}] nothing saved under {saved}: skipped')
            continue
        if not data.is_dir():
            print(f'[{dataset}] no chain outputs in {data}: skipped')
            continue
        old_files, new_files = set(output_files(saved)), set(output_files(data))
        counts = {'same': 0, 'close': 0, 'different': 0}
        lines = []
        for rel in sorted(old_files - new_files):
            lines.append(f'  MISSING   {rel.as_posix()}')
        for rel in sorted(new_files - old_files):
            lines.append(f'  ADDED     {rel.as_posix()}')
        for rel in sorted(old_files & new_files):
            try:
                status, details = compare_file(saved / rel, data / rel)
            except Exception as e:
                status, details = 'different', [f'could not compare: {e!r}']
            counts[status] += 1
            if status != 'same':
                tag = 'DIFFERENT' if status == 'different' else 'CLOSE    '
                lines.append(f'  {tag} {rel.as_posix()}')
                lines += [f'      {d}' for d in details[:5]]
                if len(details) > 5:
                    lines.append(f'      ... {len(details) - 5} more')
        n_bad_ds = counts['different'] + len(old_files ^ new_files)
        n_bad += n_bad_ds
        print(f'[{dataset}] {len(old_files & new_files)} files compared: {counts["same"]} identical, '
              f'{counts["close"]} within the tolerance (rtol {RTOL}, atol {ATOL}), {counts["different"]} '
              f'different, {len(old_files - new_files)} missing, {len(new_files - old_files)} added')
        for line in lines:
            print(line)
    print('RESULT:', 'no difference' if n_bad == 0 else f'{n_bad} file(s) differ, missing or added')
    return n_bad


if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] not in ('save', 'compare'):
        print(__doc__)
        sys.exit(2)
    if sys.argv[1] == 'save':
        save(sys.argv[2])
    else:
        sys.exit(1 if compare(sys.argv[2]) else 0)
