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


def run_tool6(data):
    with NotebookSession(T6) as nb:
        nb.pick('fc_folder', data)
        nb.pick('fc_config', Path(data) / 'config_param' / 'remap_reref_persubject.json')
        nb.click('btn_run')
        return nb.widget_text()


def run_tool7(data, events=True):
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


def run_tool9(data, clean_dir):
    data = Path(data)
    with NotebookSession(T9) as nb:
        nb.pick('fc_data', data)
        nb.pick('fc_clean', clean_dir)
        nb.pick('fc_subj', data / 'participants.tsv')
        nb.set('dd_join_col', 'participant_id')
        nb.click('btn_scan')
        nb.set('cb_avg_lin', True)          # export both averaging spaces (more columns covered)
        nb.set('cb_thirds', True)           # and the night-thirds companion table
        nb.click('btn_run')
        return nb.widget_text()


def run_chain(data, manual_participant, log=print):
    """The full chain, in order, with the default parameters of each tool. Tool 8 (interactive, one
    participant at a time) is driven on `manual_participant`."""
    data = Path(data)
    logs = {}
    for name, fn in [('5', lambda: run_tool5(data)),
                     ('5bis', lambda: run_tool5bis(data, {'age': (None, 60)})),
                     ('6', lambda: run_tool6(data)),
                     ('7', lambda: run_tool7(data)),
                     ('8bis', lambda: run_tool8bis(data)),
                     ('8', lambda: run_tool8(data, manual_participant)),
                     ('9', lambda: run_tool9(data, data / 'derivatives' / 'clean_epo_auto'))]:
        log(f'--- tool {name}')
        logs[name] = fn()
    return logs
