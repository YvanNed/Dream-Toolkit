"""The rules of tools/edf_io_lib.py, lifted from the batch tools (tool-0 redesign): the hypnogram length
rule and the companion-suffix choices of tool 7, the context-channel guess of tool 2, tool 1's header table."""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import edf_io_lib as EIO     # noqa: E402

pytestmark = [pytest.mark.quick, pytest.mark.tool0]


def test_hypnogram_length_rule_of_tool7():
    assert EIO.hypno_length_status(480, 480) == 'ok'
    assert EIO.hypno_length_status(481, 480) == 'trim'        # recording ends inside its last epoch
    assert EIO.hypno_length_status(482, 480) == 'mismatch'    # tool 7 skips the participant
    assert EIO.hypno_length_status(479, 480) == 'mismatch'


def test_mid_recording_question_marks_of_tool3():
    hyp = ['?'] * 3 + ['W'] * 20 + ['?'] + ['N2'] * 20 + ['?'] * 2
    assert EIO.mid_recording_unscored(hyp) == [23]         # the ones near both ends are lights off / on


def test_hypnogram_suffix_prefers_the_most_specific_and_never_an_event_export():
    edfs = [Path(f'p{i}.edf') for i in range(4)]
    txt = ([Path(f'p{i}_Hypnogram_Export.txt') for i in range(4)]
           + [Path(f'p{i}_Hypnogram_remapped.txt') for i in range(3)]
           + [Path(f'p{i}_ScoredEvents_Export.txt') for i in range(4)])
    counts = EIO.companion_suffix_counts(edfs, txt)
    assert EIO.select_hypno_suffix(counts) == ('_Hypnogram_remapped.txt', 3)   # >= 50 % of 4, longest
    assert EIO.select_hypno_suffix({'_ScoredEvents_Export.txt': 4}) == ('_ScoredEvents_Export.txt', 4)
    assert EIO.select_hypno_suffix({}) == (None, 0)
    assert EIO.select_event_suffix({'_event_xml.csv': 4, '_ev.csv': 4, '_x.csv': 1}) == ('_ev.csv', 4)


def test_context_channels_guess_of_tool2():
    df = pd.DataFrame({'channel': ['C3', 'E2', 'E1', 'Menton', 'ECG'],
                       'transducer_type': ['', '', '', '', '']})
    assert EIO.detect_context_channels(df) == {'eog_left': 'E1', 'eog_right': 'E2', 'emg': 'Menton', 'ecg': 'ECG'}
    # no transducer column (Curry): names only
    assert EIO.detect_context_channels(df[['channel']])['emg'] == 'Menton'


def test_header_table_of_tool1():
    hdr = {'channel': ['C3', 'EOG'], 'transducer_type': ['', ''], 'dimension': ['uV', 'mV'],
           'physical_min': ['-250', '-1'], 'physical_max': ['250', '1'], 'digital_min': ['-32768', '-32768'],
           'digital_max': ['32767', '32767'], 'prefiltering': ['HP:0.3Hz LP:35Hz', ''],
           'samples_per_record': ['256', '256'], 'sampling_frequency': ['256', '256']}
    df = EIO.header_channel_table(hdr)
    assert df['res_theoretical'].round(6).tolist() == [round(500 / 65536, 6), round(2 / 65536 * 1000, 6)]
    assert df.loc[0, 'highpass'] == 0.3 and df.loc[0, 'lowpass'] == 35.0
