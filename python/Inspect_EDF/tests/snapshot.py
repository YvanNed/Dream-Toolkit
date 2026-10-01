"""Snapshot the outputs of a pipeline run, and compare a snapshot to the golden one.

A snapshot keeps the machine-readable outputs only (TSV + JSON, gzipped) and summarises every MNE
`.fif` as a small JSON (epoch count, channels, bads, sfreq) + its metadata table. HTML / xlsx / PNG are
human outputs and are not compared. Absolute paths of the run folder are replaced by `<DATA>` so two
runs in different folders compare equal.

The golden snapshot (tests/golden/) is taken once on the code BEFORE a change: comparing a later run
to it proves the change left the existing numbers alone.
"""
import gzip
import io
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

# Columns / JSON keys that legitimately differ between two identical runs.
VOLATILE = {'validated_at', 'decided_at', 'reviewed_at', 'created_at', 'processed_at'}
INPUT_CONFIG = {'remap_reref_persubject.json', 'event_remap.json'}


def _is_output(rel):
    first = rel.parts[0]
    if first == 'derivatives' or first.startswith('reports_'):
        return True
    return first == 'config_param' and rel.name not in INPUT_CONFIG


def _replace_paths(text, data):
    # Both the given form and the resolved one: on Windows a temp folder may be passed in its short
    # 8.3 form (YVAN~1.NED) while the tools write the long one (yvan.nedelec).
    variants = set()
    for d in {Path(data), Path(os.path.realpath(str(data)))}:
        variants |= {str(d), d.as_posix(), str(d).replace('\\', '\\\\')}
    for variant in sorted(variants, key=len, reverse=True):
        text = text.replace(variant, '<DATA>')
    return text


def take(data, dest):
    """Write the snapshot of the outputs found under `data` into `dest` (replaced)."""
    import mne
    data, dest = Path(data), Path(dest)
    if dest.exists():
        # Delete the files, then the folders on a best-effort basis: OneDrive briefly locks a
        # folder it is syncing, and an empty leftover folder is harmless.
        for f in dest.rglob('*'):
            if f.is_file():
                f.unlink()
        shutil.rmtree(dest, ignore_errors=True)
    for f in sorted(data.rglob('*')):
        if not f.is_file():
            continue
        rel = f.relative_to(data)
        if not _is_output(rel):
            continue
        out = dest / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        if f.suffix in ('.tsv', '.json'):
            text = _replace_paths(f.read_text(encoding='utf-8'), data)
            with gzip.open(str(out) + '.gz', 'wt', encoding='utf-8', newline='') as fh:
                fh.write(text)
        elif f.name.endswith('-epo.fif'):
            ep = mne.read_epochs(str(f), preload=False, verbose=False)
            summary = {'n_epochs': len(ep), 'ch_names': list(ep.ch_names), 'bads': list(ep.info['bads']),
                       'sfreq': float(ep.info['sfreq']), 'tmin': float(ep.tmin), 'tmax': float(ep.tmax)}
            (Path(str(out) + '.summary.json')).write_text(json.dumps(summary, indent=1), encoding='utf-8')
            if ep.metadata is not None:
                with gzip.open(str(out) + '.metadata.tsv.gz', 'wt', encoding='utf-8', newline='') as fh:
                    ep.metadata.to_csv(fh, sep='\t', index=False)
    return dest


def _read_table(path):
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        return pd.read_csv(io.StringIO(fh.read()), sep='\t')


def _read_json(path):
    if str(path).endswith('.gz'):
        with gzip.open(path, 'rt', encoding='utf-8') as fh:
            return json.load(fh)
    return json.loads(Path(path).read_text(encoding='utf-8'))


def compare_tables(new, old, rtol=1e-6, atol=1e-9, ignore=()):
    """Differences between two DataFrames, as readable strings. Columns ADDED in `new` are allowed
    (additive changes); a column removed, a row-count change or a changed value is reported."""
    diffs = []
    missing = [c for c in old.columns if c not in new.columns]
    if missing:
        diffs.append(f'columns removed: {missing}')
    if len(new) != len(old):
        return diffs + [f'row count {len(old)} -> {len(new)}']
    for c in old.columns:
        if c in missing or c in VOLATILE or c in ignore:
            continue
        a, b = old[c], new[c]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            av, bv = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
            bad = ~np.isclose(av, bv, rtol=rtol, atol=atol, equal_nan=True)
        else:
            bad = (a.fillna('<NA>').astype(str).to_numpy() != b.fillna('<NA>').astype(str).to_numpy())
        if bad.any():
            i = int(np.flatnonzero(bad)[0])
            diffs.append(f'column {c!r}: {int(bad.sum())} value(s) differ, first at row {i}: '
                         f'{a.iloc[i]!r} -> {b.iloc[i]!r}')
    return diffs


def _compare_json(new, old, where=''):
    diffs = []
    if isinstance(old, dict) and isinstance(new, dict):
        for k in old:
            if k in VOLATILE:
                continue
            if k not in new:
                diffs.append(f'{where}/{k}: key removed')
            else:
                diffs += _compare_json(new[k], old[k], f'{where}/{k}')
    elif isinstance(old, float) and isinstance(new, (int, float)):
        if not np.isclose(old, new, rtol=1e-6, atol=1e-9, equal_nan=True):
            diffs.append(f'{where}: {old!r} -> {new!r}')
    elif old != new:
        diffs.append(f'{where}: {old!r} -> {new!r}')
    return diffs


def compare(new_dir, golden_dir, prefixes=None, expect_missing=()):
    """Compare snapshot `new_dir` to `golden_dir`. Returns {relative file: [differences]} (empty dict =
    identical). `prefixes` restricts to golden files whose relative path starts with one of them;
    `expect_missing` lists relative paths (or prefixes) a change is ALLOWED to remove."""
    new_dir, golden_dir = Path(new_dir), Path(golden_dir)
    report = {}
    for g in sorted(golden_dir.rglob('*')):
        if not g.is_file():
            continue
        rel = g.relative_to(golden_dir).as_posix()
        if prefixes and not any(rel.startswith(p) for p in prefixes):
            continue
        n = new_dir / rel
        if not n.exists():
            if not any(rel.startswith(p) for p in expect_missing):
                report[rel] = ['file missing in the new run']
            continue
        if rel.endswith('.tsv.gz'):
            d = compare_tables(_read_table(n), _read_table(g))
        else:
            d = _compare_json(_read_json(n), _read_json(g))
        if d:
            report[rel] = d
    return report


def format_report(report):
    return '\n'.join(f'{f}:\n  ' + '\n  '.join(d) for f, d in report.items())
