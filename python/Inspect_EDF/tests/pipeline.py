"""Run the toolkit chain 5 -> 5bis -> 6 -> 7 -> 8bis -> 8 -> 9 on the test database, each tool in
its own kernel, driven through its widgets like a user would (see nbdriver.py).

Every `run_toolX(data, ...)` returns the text the tool printed, so a test can assert on warnings.
Paths follow the toolkit's default layout: outputs inside the data folder.
"""
from pathlib import Path

from nbdriver import NotebookSession

T5 = 'tools/5_compute_sleep_macrostructure_voila.ipynb'
T5BIS = 'tools/5bis_exclude_sleep_macrostructure_voila.ipynb'
T6 = 'tools/6_quality_overview_voila.ipynb'
T7 = 'tools/7_preprocessing_voila.ipynb'
T8 = 'tools/8_reject_manually_voila.ipynb'
T8BIS = 'tools/8bis_reject_automatically_voila.ipynb'
T9 = 'tools/9_spectral_features_voila.ipynb'


def run_tool5(data):
    with NotebookSession(T5) as nb:
        nb.pick('fc_data', data)
        nb.pick('fc_subj', Path(data) / 'participants.tsv')
        nb.set('dd_join_col', 'participant_id')
        nb.click('btn_scan')
        nb.click('btn_run')
        return nb.widget_text()


def run_tool5bis(data, criteria=None, comment='test selection'):
    """criteria = {metric: (min or None, max or None)}, e.g. {'age': (None, 60)}."""
    with NotebookSession(T5BIS) as nb:
        nb.pick('fc_data', data)
        nb.click('btn_load')
        for metric, (lo, hi) in (criteria or {}).items():
            if lo is not None:
                nb.run(f'ROWS[{metric!r}]["tmin"].value = {str(lo)!r}')
            if hi is not None:
                nb.run(f'ROWS[{metric!r}]["tmax"].value = {str(hi)!r}')
            nb.run(f'ROWS[{metric!r}]["cb"].value = True')
        nb.set('txt_comment', comment)
        nb.click('btn_validate')
        return nb.widget_text() + nb.value('lbl_validate.value')


def run_tool6(data, subset=None, skip=True):
    with NotebookSession(T6) as nb:
        nb.pick('fc_folder', data)
        nb.pick('fc_config', Path(data) / 'config_param' / 'remap_reref_persubject.json')
        if subset is not None:
            select_subset(nb, subset, skip)
        nb.click('btn_run')
        return nb.widget_text()


def select_subset(nb, ids, skip):
    """Put a tool's shared participant selector in 'Subset only' mode on `ids`."""
    nb.run("participant_selector.dd_mode.value = 'Subset only'")
    nb.run(f'participant_selector.tags.value = {list(ids)!r}')
    nb.set('participant_selector.cb_skip', skip)


def run_tool7(data, events=True, subset=None, skip=True, channels=None, before_run=None):
    """subset: run only these participants; channels: {file_id: {channel: keep}} edits applied in the
    per-participant editor; before_run(nb): any extra driving before Run."""
    data = Path(data)
    with NotebookSession(T7) as nb:
        nb.pick('fc_edf', data)
        nb.pick('fc_quality', data / 'reports_quality_overview' / 'quality_summary.tsv')
        nb.pick('fc_config', data / 'config_param' / 'remap_reref_persubject.json')
        nb.pick('fc_events', data / 'config_param' / 'event_remap.json')
        nb.pick('fc_out', data)
        if events:
            nb.set('cb_event_reject', True)
        nb.click('btn_load')
        if subset is not None:
            select_subset(nb, subset, skip)
        for fid, chans in (channels or {}).items():
            for ch, keep in chans.items():
                nb.run(f'participant_channels[{fid!r}][{ch!r}] = {keep!r}')
        if before_run is not None:
            before_run(nb)
        nb.click('btn_run')
        return nb.widget_text()


def run_tool8bis(data):
    data = Path(data)
    with NotebookSession(T8BIS) as nb:
        nb.pick('fc_raw', data / 'derivatives' / 'raw_epo')
        nb.pick('fc_reports', data / 'reports_preprocessing')
        nb.click('btn_scan')
        nb.click('btn_run')
        return nb.widget_text()


def run_tool8(data, file_id):
    """Load one participant and save its decision with no manual override (the default decision)."""
    data = Path(data)
    with NotebookSession(T8) as nb:
        nb.pick('fc_data', data)
        nb.pick('fc_raw', data / 'derivatives' / 'raw_epo')
        nb.pick('fc_reports', data / 'reports_preprocessing')
        nb.set('dd_part', file_id)
        nb.click('btn_load')
        nb.click('btn_save')
        return nb.widget_text()


def run_tool9(data, raw_dir, decision_dir=None, subset=None, skip=True, reaggregate=False):
    """decision_dir: a tool 8 / 8bis rejection_* folder (None = no decision). reaggregate: click
    'Re-aggregate only' instead of Run."""
    data = Path(data)
    with NotebookSession(T9) as nb:
        nb.pick('fc_data', data)
        nb.pick('fc_raw', raw_dir)
        if decision_dir is not None:
            nb.pick('fc_decision', decision_dir)
        nb.pick('fc_subj', data / 'participants.tsv')
        nb.set('dd_join_col', 'participant_id')
        nb.click('btn_scan')
        nb.set('cb_avg_lin', True)          # export both averaging spaces (more columns covered)
        nb.set('cb_thirds', True)           # and the night-thirds companion table
        if subset is not None:
            select_subset(nb, subset, skip)
        nb.click('btn_reagg' if reaggregate else 'btn_run')
        return nb.widget_text()


def chain_steps(data, manual_participant):
    """The chain, in order, as [(tool, callable)], each tool with its default parameters. Tool 8 (interactive,
    one participant at a time) is driven on `manual_participant`."""
    data = Path(data)
    return [('5', lambda: run_tool5(data)),
            ('5bis', lambda: run_tool5bis(data, {'age': (None, 60)})),
            ('6', lambda: run_tool6(data)),
            ('7', lambda: run_tool7(data)),
            ('8bis', lambda: run_tool8bis(data)),
            ('8', lambda: run_tool8(data, manual_participant)),
            ('9', lambda: run_tool9(data, data / 'derivatives' / 'raw_epo',
                                    data / 'derivatives' / 'rejection_auto'))]


def run_chain(data, manual_participant, log=print):
    """The full chain from scratch (used by make_golden.py; the test fixture reuses the last run instead,
    see chaincache.py)."""
    logs = {}
    for name, fn in chain_steps(data, manual_participant):
        log(f'--- tool {name}')
        logs[name] = fn()
    return logs
