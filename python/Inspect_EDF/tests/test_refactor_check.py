"""The before / after checker (refactor_check.py) must see every kind of change, and nothing else."""
import json

import numpy as np
import pandas as pd
import pytest

import refactor_check as RC

pytestmark = pytest.mark.quick


def _write_tree(root, table, meta, html, png=b'\x89PNG-same'):
    data = root / 'data'
    (data / 'reports_x').mkdir(parents=True, exist_ok=True)
    table.to_csv(data / 'reports_x' / 'p1_metrics.tsv', sep='\t', index=False)
    (data / 'reports_x' / 'p1_params.json').write_text(json.dumps(meta), encoding='utf-8')
    (data / 'reports_x' / 'p1_report.html').write_text(html, encoding='utf-8')
    (data / 'reports_x' / 'p1_fig.png').write_bytes(png)
    (data / 'p1.edf').write_bytes(b'recording')          # an input: never copied nor compared
    np.savez(data / 'reports_x' / 'p1_psd.npz', psd=np.array([1.0, 2.0, 3.0]))
    return data


@pytest.fixture
def roots(tmp_path, monkeypatch):
    monkeypatch.setattr(RC, 'RUN_ROOT', tmp_path / 'run')
    monkeypatch.setattr(RC, 'SAVE_ROOT', tmp_path / 'saved')
    monkeypatch.setattr(RC, 'DATASETS', ['synthetic'])
    return tmp_path / 'run' / 'synthetic'


TABLE = pd.DataFrame({'channel': ['C3', 'O1'], 'std_uV': [12.5, 30.25], 'processed_at': ['t1', 't1']})
META = {'threshold': 0.15, 'methods_run': ['amplitude', 'flat']}
HTML = '<p>Report written 2026-10-09 14:03:12</p><p>12 epochs</p>'


def test_identical_rerun_gives_no_difference(roots):
    _write_tree(roots, TABLE, META, HTML)
    RC.save('before')
    # Same content, other clock time in the report and in a volatile column: still identical.
    _write_tree(roots, TABLE.assign(processed_at='t2'), META, HTML.replace('14:03:12', '15:44:01'))
    assert RC.compare('before') == 0


def test_saved_copy_skips_the_recordings(roots):
    _write_tree(roots, TABLE, META, HTML)
    RC.save('before')
    assert not (RC.SAVE_ROOT / 'before' / 'synthetic' / 'p1.edf').exists()


@pytest.mark.parametrize('change', ['value', 'tiny', 'json', 'html', 'png', 'npz', 'missing', 'added'])
def test_every_change_is_seen(roots, change, capsys):
    data = _write_tree(roots, TABLE, META, HTML)
    RC.save('before')
    if change == 'value':
        _write_tree(roots, TABLE.assign(std_uV=[12.5, 31.0]), META, HTML)
    elif change == 'tiny':                 # within the golden tolerance: reported as CLOSE, not a failure
        _write_tree(roots, TABLE.assign(std_uV=[12.5, 30.25 * (1 + 1e-9)]), META, HTML)
    elif change == 'json':
        _write_tree(roots, TABLE, dict(META, threshold=0.2), HTML)
    elif change == 'html':
        _write_tree(roots, TABLE, META, HTML.replace('12 epochs', '13 epochs'))
    elif change == 'png':
        _write_tree(roots, TABLE, META, HTML, png=b'\x89PNG-other')
    elif change == 'npz':
        np.savez(data / 'reports_x' / 'p1_psd.npz', psd=np.array([1.0, 2.0, 3.5]))
    elif change == 'missing':
        (data / 'reports_x' / 'p1_fig.png').unlink()
    elif change == 'added':
        (data / 'reports_x' / 'p1_new.tsv').write_text('a\n1\n', encoding='utf-8')
    n_bad = RC.compare('before')
    out = capsys.readouterr().out
    if change == 'tiny':
        assert n_bad == 0 and 'CLOSE' in out
    else:
        assert n_bad == 1, out
