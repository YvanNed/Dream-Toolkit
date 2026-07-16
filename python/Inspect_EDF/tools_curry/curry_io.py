"""
curry_io.py — signal, hypnogram and events loading for Curry 9 .cdt files.

Three public functions:
    read_curry_signal(cdt_path, include=None, preload=False)
        → mne.io.RawCurry   (wraps mne.io.read_raw_curry)

    load_hypnogram_curry(txt_path)
        → list[str]         (one label per epoch, same format as EDF tools)

    load_events_curry(txt_path, rec_start_dt)
        → pd.DataFrame with columns [Name, Start, Duration]
          where Start and Duration are in seconds from recording start.
          Returns an empty DataFrame on failure.
"""

import re
import datetime
import pandas as pd
import mne


# ---------------------------------------------------------------------------
# Signal loading
# ---------------------------------------------------------------------------

def read_curry_signal(cdt_path: str,
                      include: list | None = None,
                      preload: bool = False) -> mne.io.BaseRaw:
    """
    Load a Curry .cdt file via MNE 1.12+ (requires curryreader).

    Unlike the EDF loader, all Curry channels share a single sampling rate,
    so channel selection is a simple raw.pick() — no include-at-read-time
    workaround needed.

    Parameters
    ----------
    cdt_path : str
        Path to the .cdt file.
    include : list[str] | None
        Channel names to keep (remapped or original). If None, all channels
        are kept. Applied after loading.
    preload : bool
        Whether to load signal data into memory immediately.

    Returns
    -------
    mne.io.BaseRaw
        Raw object with the requested channels. Signal is in MNE's native
        unit (V); multiply by 1e6 to get µV as usual.
    """
    raw = mne.io.read_raw_curry(str(cdt_path), preload=preload, verbose="ERROR")
    if include is not None:
        present = [ch for ch in include if ch in raw.ch_names]
        raw.pick(present)
    return raw


# ---------------------------------------------------------------------------
# Hypnogram loading
# ---------------------------------------------------------------------------

def load_hypnogram_curry(txt_path: str) -> list:
    """
    Load a Curry hypnogram export (one label per line, plain text).

    Identical format to the EDF-side hypnogram tools. Labels are returned
    as-is (e.g. "W", "1", "2", "3", "R") — remapping to AASM is done by
    tool 3 as usual.

    Parameters
    ----------
    txt_path : str
        Path to the *_Hypnogram_Export.txt file.

    Returns
    -------
    list[str]
    """
    with open(txt_path, "r", encoding="utf-8", errors="replace") as fh:
        labels = [line.strip() for line in fh if line.strip()]
    return labels


# ---------------------------------------------------------------------------
# Events loading — French text export
# ---------------------------------------------------------------------------

# Duration format: "M:SS" or "M:SS.s" — left = minutes, right = seconds
_DUR_RE = re.compile(r"^(\d+):(\d+(?:\.\d+)?)$")


def load_events_curry(txt_path: str,
                      rec_start_dt: datetime.datetime) -> pd.DataFrame:
    """
    Parse a Curry scored-events text export (*_ScoredEvents_Export.txt).

    File format (comma-separated, no header):
        HH:MM:SS , epoch# , stage_FR , event_label_FR , M:SS[.s] , - , - , position

    Encoding varies between exports: UTF-16 (with BOM) on some, plain UTF-8/ANSI on others.
    The BOM is sniffed below -- do not assume UTF-16.

    The clock-time column is converted to seconds-from-recording-start with
    midnight rollover (events after midnight are on the next calendar day).

    Parameters
    ----------
    txt_path : str
        Path to the *_ScoredEvents_Export.txt file.
    rec_start_dt : datetime.datetime
        Recording start as a timezone-naive or timezone-aware datetime.
        Obtainable from curry_header: datetime.datetime.strptime(
            hdr["start_datetime"], "%Y-%m-%d %H:%M:%S.%f")

    Returns
    -------
    pd.DataFrame
        Columns: Name (str), Start (float, seconds), Duration (float, seconds).
        Empty DataFrame on failure or if the file is absent.
    """
    try:
        return _parse_events_txt(txt_path, rec_start_dt)
    except Exception as exc:
        print(f"⚠ load_events_curry({txt_path}): {exc}")
        return pd.DataFrame(columns=["Name", "Start", "Duration"])


def _parse_events_txt(txt_path: str,
                      rec_start_dt: datetime.datetime) -> pd.DataFrame:
    # Decode UTF-16 (with or without BOM)
    raw_bytes = open(txt_path, "rb").read()
    if raw_bytes[:2] in (b"\xff\xfe", b"\xfe\xff"):
        text = raw_bytes.decode("utf-16")
    else:
        text = raw_bytes.decode("utf-8", errors="replace")

    # Strip timezone from rec_start_dt for naive arithmetic
    if hasattr(rec_start_dt, "tzinfo") and rec_start_dt.tzinfo is not None:
        rec_start_naive = rec_start_dt.replace(tzinfo=None)
    else:
        rec_start_naive = rec_start_dt

    rec_date = rec_start_naive.date()

    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 5:
            continue

        clock_str = parts[0]   # "HH:MM:SS"
        label_fr  = parts[3]   # French event label
        dur_str   = parts[4]   # "M:SS" or "M:SS.s"

        # --- Parse clock time ---
        # fromisoformat requires zero-padded HH:MM:SS; the export sometimes
        # writes single-digit hours (e.g. "0:05:11", "9:30:00"), so we split manually.
        try:
            hh, mm, ss = clock_str.split(":")
            clock_t = datetime.time(int(hh), int(mm), int(ss))
        except (ValueError, TypeError):
            continue
        event_dt = datetime.datetime.combine(rec_date, clock_t)

        # Midnight rollover: if the event time is before recording start
        # (by more than 1 hour, to tolerate sub-minute rounding), it's next day.
        if (rec_start_naive - event_dt).total_seconds() > 3600:
            event_dt += datetime.timedelta(days=1)

        start_sec = (event_dt - rec_start_naive).total_seconds()

        # --- Parse duration "M:SS[.s]" ---
        m = _DUR_RE.match(dur_str)
        if m:
            dur_sec = int(m.group(1)) * 60 + float(m.group(2))
        else:
            dur_sec = 0.0

        rows.append({"Name": label_fr, "Start": start_sec, "Duration": dur_sec})

    return pd.DataFrame(rows, columns=["Name", "Start", "Duration"])


# ---------------------------------------------------------------------------
# Utility: build recording-start datetime from curry_header dict
# ---------------------------------------------------------------------------

def rec_start_from_header(hdr: dict) -> datetime.datetime | None:
    """
    Build a timezone-naive datetime from a curry_header dict.
    Returns None if the header has no start_datetime.
    """
    dt_str = hdr.get("start_datetime", "")
    if not dt_str:
        return None
    try:
        return datetime.datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S.%f")
    except ValueError:
        try:
            return datetime.datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")
        except ValueError:
            return None
