"""The signal-preparation label of tools/quality_lib.py (tool 6, tool 0's step 5): what it says, and when a skipped
participant's stored label is reported as different from the run's settings."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'tools'))
import quality_lib as QL     # noqa: E402

pytestmark = [pytest.mark.quick, pytest.mark.tool6]


def test_signal_preparation_label():
    assert QL.signal_preparation_label() == 'none'
    assert QL.signal_preparation_label(resampled_to=128.0) == 'resampled to 128 Hz'
    assert QL.signal_preparation_label(highpassed_at=0.5) == 'high-pass 0.5 Hz'
    assert QL.signal_preparation_label(128.0, 0.5) == 'resampled to 128 Hz, high-pass 0.5 Hz'


def test_signal_preparation_mismatch():
    # same settings, or a table written before the column existed: nothing to report
    assert QL.signal_preparation_mismatch('none', False, None, False, None) == ''
    assert QL.signal_preparation_mismatch('high-pass 0.5 Hz', False, None, True, 0.5) == ''
    assert QL.signal_preparation_mismatch(float('nan'), False, None, True, 0.5) == ''
    assert QL.signal_preparation_mismatch('', False, None, True, 0.5) == ''
    # a different high-pass is always reported
    assert (QL.signal_preparation_mismatch('high-pass 0.5 Hz', False, None, False, None)
            == 'made with high-pass 0.5 Hz, this run: none')
    assert QL.signal_preparation_mismatch('none', False, None, True, 0.3) == 'made with none, this run: high-pass 0.3 Hz'
    # a stored resampling is compared with the run's
    assert QL.signal_preparation_mismatch('resampled to 128 Hz', True, 128, False, None) == ''
    assert QL.signal_preparation_mismatch('resampled to 128 Hz', False, None, False, None) != ''
    assert QL.signal_preparation_mismatch('resampled to 128 Hz', True, 256, False, None) != ''
    # not resampled while this run resamples: the file may already have been at the target rate (not told)
    assert QL.signal_preparation_mismatch('none', True, 128, False, None) == ''
