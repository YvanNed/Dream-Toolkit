"""Reuse the last run of the test chain: only the tools whose code changed (and the ones after them) run again.

The chain 5 -> 5bis -> 6 -> 7 -> 8bis -> 8 -> 9 is long (minutes). When a change touches tool 9 only, re-running
5 to 8 proves nothing new. The last run is kept in its folder with a fingerprint of each tool (its notebook +
the shared modules it uses) and of the 'base' (test driver, test-data generator). Next time:

    base changed, no state, unreadable state, --fresh-chain   -> full run from scratch
    otherwise                                                 -> rerun from the FIRST tool whose fingerprint
                                                                 changed: its outputs and those of every later
                                                                 tool are deleted first (+ their registry rows)

The state file is updated after each tool succeeds, so an interrupted run resumes where it stopped. The logs of
the tools that are not rerun are kept, so 'no traceback in any tool' still covers them.
"""
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
TOOLS = REPO / 'tools'
sys.path.insert(0, str(TOOLS))
import participant_selection_lib as PSL     # noqa: E402

STATE_NAME = 'chain_state.json'
STATE_VERSION = 1

# What a change to which file means for the chain
BASE_DEPS = [HERE / 'nbdriver.py', HERE / 'pipeline.py', HERE / 'minidb.py', HERE / 'synthetic.py',
             HERE / 'chaincache.py', TOOLS / 'generate_test_data.py']
LIB_PSL, LIB_QC = TOOLS / 'participant_selection_lib.py', TOOLS / 'qc_rejected_epochs_lib.py'
TOOL_DEPS = {
    '5': [TOOLS / '5_compute_sleep_macrostructure_voila.ipynb', LIB_PSL],
    '5bis': [TOOLS / '5bis_exclude_sleep_macrostructure_voila.ipynb', LIB_PSL],
    '6': [TOOLS / '6_quality_overview_voila.ipynb', LIB_PSL],
    '7': [TOOLS / '7_preprocessing_voila.ipynb', LIB_PSL],
    '8bis': [TOOLS / '8bis_reject_automatically_voila.ipynb', LIB_PSL, LIB_QC],
    '8': [TOOLS / '8_reject_manually_voila.ipynb', LIB_PSL, LIB_QC],
    '9': [TOOLS / '9_spectral_features_voila.ipynb', LIB_PSL],
}
# What each tool writes in the data folder (deleted before it reruns)
TOOL_OUTPUTS = {
    '5': ['derivatives/features_macrostructure', 'reports_features_macrostructure/global_sleep_metrics.tsv',
          'reports_features_macrostructure/global_sleep_metrics_checks.tsv',
          'reports_features_macrostructure/sleep_metrics_database.xlsx',
          'reports_features_macrostructure/sleep_metrics_database_report.html',
          'reports_features_macrostructure/sleep_metrics_failed.tsv'],
    '5bis': ['reports_features_macrostructure/participant_selection.tsv',
             'reports_features_macrostructure/participant_selection_criteria.tsv',
             'reports_features_macrostructure/participant_selection_report.html'],
    '6': ['reports_quality_overview'],
    '7': ['derivatives/raw_epo', 'reports_preprocessing'],
    '8bis': ['derivatives/rejection_auto', 'reports_rejection_auto'],
    '8': ['derivatives/rejection_manual', 'reports_rejection_manual'],
    '9': ['derivatives/features_spectral', 'reports_features_spectral'],
}
# The registry rows each tool writes (removed before it reruns)
TOOL_REGISTRY_SOURCES = {'5bis': ['5bis_criteria', '5bis_manual'], '7': ['7_manual'],
                         '8bis': ['8bis_auto'], '8': ['8_manual']}


def _hash(paths):
    h = hashlib.sha256()
    for p in paths:
        h.update(p.name.encode())
        h.update(p.read_bytes() if p.is_file() else b'<missing>')
    return h.hexdigest()


def fingerprints(dataset):
    base = _hash(BASE_DEPS) + dataset
    return base, {tool: _hash(deps) for tool, deps in TOOL_DEPS.items()}


def _load_state(root):
    try:
        state = json.loads((Path(root) / STATE_NAME).read_text(encoding='utf-8'))
        return state if state.get('version') == STATE_VERSION else None
    except Exception:
        return None


def _save_state(root, state):
    tmp = Path(root) / (STATE_NAME + '.tmp')
    tmp.write_text(json.dumps(state, indent=1), encoding='utf-8')
    tmp.replace(Path(root) / STATE_NAME)


def first_tool_to_rerun(order, state, base, tool_fp):
    """Index in `order` of the first tool that must run again (len(order) = nothing to rerun), or None when
    the whole chain must restart from scratch (no usable state, base changed)."""
    if state is None or state.get('base') != base:
        return None
    done = state.get('tools', {})
    for i, tool in enumerate(order):
        if done.get(tool) != tool_fp[tool]:
            return i
    return len(order)


def rewind(data, tools):
    """Delete the outputs and the registry rows of `tools` (the ones about to rerun)."""
    data = Path(data)
    for tool in tools:
        for rel in TOOL_OUTPUTS[tool]:
            p = data / rel
            if p.is_dir():
                shutil.rmtree(p)
            elif p.exists():
                p.unlink()
        for source in TOOL_REGISTRY_SOURCES.get(tool, []):
            if PSL.registry_path(data).is_file():
                PSL.write_registry_rows(data, source, [])


def run_cached_chain(dataset, root, build, steps_for, fresh=False, log=print):
    """Run the chain of `dataset` in `root`, reusing the last run when possible.

    build(root) -> data folder (fresh test database); steps_for(data) -> ordered [(tool, callable)]."""
    root = Path(root)
    base, tool_fp = fingerprints(dataset)
    state = None if fresh else _load_state(root)
    probe_steps = [t for t, _ in steps_for(root / 'data')]
    start = first_tool_to_rerun(probe_steps, state, base, tool_fp)
    t0 = time.time()
    if start is None:
        if root.exists():
            shutil.rmtree(root)
        data = build(root)
        state = {'version': STATE_VERSION, 'base': base, 'tools': {}, 'logs': {}}
        start = 0
        log(f'[{dataset}] chain: full run')
    else:
        data = root / 'data'
        rerun = probe_steps[start:]
        log(f'[{dataset}] chain: reusing {", ".join(probe_steps[:start]) or "nothing"}; '
            f'rerunning {", ".join(rerun) or "nothing"}')
        rewind(data, rerun)
        for tool in rerun:
            state['tools'].pop(tool, None)
    for tool, run in steps_for(data)[start:]:
        log(f'[{dataset} {time.time() - t0:5.0f}s] --- tool {tool}')
        state['logs'][tool] = run()
        state['tools'][tool] = tool_fp[tool]
        _save_state(root, state)
    return data, dict(state['logs'])
