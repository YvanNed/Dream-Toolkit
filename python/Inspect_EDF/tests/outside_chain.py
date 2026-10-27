r"""Tools that are not part of the test chain (1, 3, 3bis), or a tool with options the chain leaves at their
defaults (6variant = tool 6 with resampling + high-pass; 7variant = tool 7 with resampling, 10 s epochs, no notch,
LOWESS-only smoothing), run on a copy of the INPUTS of a chain dataset, for the before / after check of the
tool-0 redesign (code moved into tools/*_lib.py).

    & $py tests/outside_chain.py [1] [3] [3bis] [6variant] [7variant]   # build %TEMP%/dtk_outside_chain/<dataset>/data
    $env:DTK_TEST_ROOT = "$env:TEMP\dtk_outside_chain"; $env:DTK_REFACTOR_ROOT = "$env:TEMP\dtk_outside_saved"
    & $py tests/refactor_check.py save before    # then, after the move: run this script again + compare before

Without an argument, tool 3bis only.

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
T1 = 'tools/1_inspect_edf_voila.ipynb'
T3 = 'tools/3_remap_hypno_voila.ipynb'
T3BIS = 'tools/3bis_autoscoring_voila.ipynb'
T6 = 'tools/6_quality_overview_voila.ipynb'
T7 = 'tools/7_preprocessing_voila.ipynb'
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


def run_tool1(data):
    """Tool 1 on every recording, with its defaults (EEG + EOG inspected)."""
    with NotebookSession(T1) as nb:
        nb.pick('chooser', data)
        nb.click('run_button')
        return nb.widget_text()


def run_tool3(data):
    """Tool 3 on the raw hypnogram exports when there are some (else the remapped hypnograms), with its default
    rules, written under its own suffix (the chain's hypnograms stay untouched)."""
    suffix = '_Hypnogram_Export.txt' if any(Path(data).glob('*_Hypnogram_Export.txt')) else '_Hypnogram_remapped.txt'
    with NotebookSession(T3) as nb:
        nb.pick('chooser', data)
        nb.set('suffix_box', suffix)
        nb.set('output_suffix_box', '_Hypnogram_tool3.txt')
        nb.click('btn_scan')
        nb.click('btn_run_remap')
        nb.click_labelled('Save remapping')
        nb.click('btn_run_save')
        nb.click_labelled('Save files')
        nb.click('btn_run_verify')
        return nb.widget_text()


def run_tool3bis(data):
    """Score every participant with the defaults (GSSC alone, consensus over C3/C4/F3/F4 present: C3 on the test
    datasets, EOG L-R when declared, contralateral mastoid reference)."""
    with NotebookSession(T3BIS) as nb:
        nb.pick('fc_data', data)
        nb.click('btn_scan')
        nb.click('btn_run')
        return nb.widget_text()


def run_tool6_variant(data):
    """Tool 6 with its two optional signal steps on (resampling to 128 Hz, high-pass 0.5 Hz)."""
    with NotebookSession(T6) as nb:
        nb.pick('fc_folder', data)
        nb.pick('fc_config', Path(data) / 'config_param' / 'remap_reref_persubject.json')
        nb.set('cb_resample', True)
        nb.set('txt_target_freq', 128)
        nb.set('hp_check', True)
        nb.set('hp_freq', 0.5)
        nb.click('btn_run')
        return nb.widget_text()


def run_tool7_variant(data):
    """Tool 7 with the options the chain leaves at their defaults: resampling to 128 Hz, 10 s epochs, no notch,
    PSD smoothing by LOWESS only, events flagged. Needs no quality table (channels from the configuration)."""
    data = Path(data)
    with NotebookSession(T7) as nb:
        nb.pick('fc_edf', data)
        nb.pick('fc_config', data / 'config_param' / 'remap_reref_persubject.json')
        nb.pick('fc_events', data / 'config_param' / 'event_remap.json')
        nb.pick('fc_out', data)
        nb.set('cb_event_reject', True)
        nb.set('cb_resample', True)
        nb.set('txt_target_freq', 128)
        nb.set('dd_epoch_len', 10)
        nb.set('cb_notch', False)
        nb.set('cb_ss_median', False)
        nb.click('btn_load')
        nb.click('btn_run')
        return nb.widget_text()


RUNNERS = {'1': run_tool1, '3': run_tool3, '3bis': run_tool3bis, '6variant': run_tool6_variant,
           '7variant': run_tool7_variant}


if __name__ == '__main__':
    tools = sys.argv[1:] or ['3bis']
    for dataset in ('synthetic', 'real'):
        if not (CHAIN_ROOT / dataset / 'data').is_dir():
            print(f'[{dataset}] no chain folder: run the test suite first')
            continue
        data = build(dataset)
        for tool in tools:
            t0 = time.time()
            text = RUNNERS[tool](data)
            bad = 'Traceback (most recent call last)' in text
            print(f'[{dataset}] tool {tool} done in {time.time() - t0:.0f} s' + ('  ⚠ TRACEBACK in its log' if bad else ''))
            if bad:
                print(text[-3000:])
