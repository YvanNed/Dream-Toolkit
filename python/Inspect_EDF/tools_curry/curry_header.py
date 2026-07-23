"""
curry_header.py — lightweight ASCII parser for Curry 9 .cdt.dpo parameter files.

Reads channel labels, sampling frequency, start date/time, data unit, impedances,
and sensor positions WITHOUT loading the binary .cdt signal data.

The .cdt.dpo is a plain-text file (~7 kB) structured as named blocks:
    BLOCK_NAME START
       Key = Value
    BLOCK_NAME END

    BLOCK_NAME START_LIST   # Do not edit!
       <rows of data>
    BLOCK_NAME END_LIST

Usage
-----
    from curry_header import read_curry_header
    hdr = read_curry_header("path/to/subject.cdt")
    # hdr['sfreq'], hdr['n_channels'], hdr['ch_labels'], hdr['impedances'], …
"""

import os
import re
import chardet


def read_curry_header(cdt_path: str) -> dict:
    """
    Parse a Curry .cdt.dpo parameter file and return a header dict.

    Parameters
    ----------
    cdt_path : str
        Path to the .cdt data file (the .cdt.dpo companion is found automatically).

    Returns
    -------
    dict with keys:
        file_id        : str   — stem of the .cdt file
        dpo_path       : str   — resolved .cdt.dpo path
        sfreq          : float — global sampling frequency (Hz)
        n_samples      : int   — total number of samples
        n_channels     : int   — total channel count (EEG + others)
        duration_sec   : float — n_samples / sfreq
        n_epochs_30s   : int   — floor(duration_sec / 30)
        start_datetime : str   — "YYYY-MM-DD HH:MM:SS.mmm" (empty string if absent)
        data_unit      : str   — e.g. "uV" (from first DEVICE_PARAMETERS block)
        file_version   : str   — FileVersion field
        byte_order     : str   — DataByteOrder field
        data_format    : str   — DataFormat field
        ch_labels      : list[str]  — EEG channel labels (from LABELS block)
        ch_labels_other: list[str]  — other-group channel labels (LABELS_OTHERS)
        all_ch_labels  : list[str]  — ch_labels + ch_labels_other (full order)
        eeg_group_size : int   — number of channels in the EEG group
        other_group_size: int  — number of channels in the other group
        impedances     : list[list[float]]  — shape (n_channels, n_sessions);
                         -1 = not measured, -2 = disabled/not connected
        sensor_xyz     : list[tuple]  — (x, y, z) mm for each EEG channel (may be [])
        raw_params     : dict  — all key=value pairs from DATA_PARAMETERS block
    """
    cdt_path = str(cdt_path)
    dpo_path = cdt_path + ".dpo"
    if not os.path.isfile(dpo_path):
        raise FileNotFoundError(f"Curry parameter file not found: {dpo_path}")

    file_id = os.path.splitext(os.path.basename(cdt_path))[0]
    # strip second extension (.cdt) if file_id still ends with it
    if file_id.endswith(".cdt"):
        file_id = file_id[:-4]

    text = _read_text(dpo_path)
    blocks = _parse_blocks(text)

    # --- DATA_PARAMETERS (global metadata) ---
    params = blocks.get("DATA_PARAMETERS", {})
    sfreq        = float(params.get("SampleFreqHz",  0) or 0)
    n_samples    = int(params.get("NumSamples",      0) or 0)
    n_channels   = int(params.get("NumChannels",     0) or 0)
    file_version = str(params.get("FileVersion", "")).strip()
    byte_order   = str(params.get("DataByteOrder", "")).strip()
    data_format  = str(params.get("DataFormat", "")).strip()

    start_datetime = _parse_start_datetime(params)
    duration_sec   = (n_samples / sfreq) if sfreq else 0.0
    n_epochs_30s   = int(duration_sec // 30)

    # --- DEVICE_PARAMETERS (units for EEG group) ---
    dev_params = blocks.get("DEVICE_PARAMETERS", {})
    data_unit  = str(dev_params.get("DataUnit", "")).strip()
    eeg_group_size = int(dev_params.get("NumChanThisGroup", 0) or 0)

    dev_other  = blocks.get("DEVICE_PARAMETERS_OTHERS", {})
    other_group_size = int(dev_other.get("NumChanThisGroup", 0) or 0)

    # --- Channel labels ---
    ch_labels       = blocks.get("LABELS_list",        [])
    ch_labels_other = blocks.get("LABELS_OTHERS_list", [])
    all_ch_labels   = ch_labels + ch_labels_other

    # --- Impedances ---
    impedances = blocks.get("IMPEDANCE_VALUES_list", [])

    # --- Sensor positions (EEG group) ---
    sensor_xyz = blocks.get("SENSORS_list", [])

    return {
        "file_id":         file_id,
        "dpo_path":        dpo_path,
        "sfreq":           sfreq,
        "n_samples":       n_samples,
        "n_channels":      n_channels,
        "duration_sec":    duration_sec,
        "n_epochs_30s":    n_epochs_30s,
        "start_datetime":  start_datetime,
        "data_unit":       data_unit,
        "file_version":    file_version,
        "byte_order":      byte_order,
        "data_format":     data_format,
        "ch_labels":       ch_labels,
        "ch_labels_other": ch_labels_other,
        "all_ch_labels":   all_ch_labels,
        "eeg_group_size":  eeg_group_size,
        "other_group_size":other_group_size,
        "impedances":      impedances,
        "sensor_xyz":      sensor_xyz,
        "raw_params":      params,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_text(path: str) -> str:
    """Read a file, auto-detecting encoding (UTF-8-BOM, UTF-16, latin-1 fallback)."""
    raw_bytes = open(path, "rb").read()
    # BOM detection first (faster and more reliable than chardet for these files)
    if raw_bytes[:3] == b"\xef\xbb\xbf":
        return raw_bytes[3:].decode("utf-8")
    if raw_bytes[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw_bytes.decode("utf-16")
    # chardet for anything else
    detected = chardet.detect(raw_bytes)
    enc = detected.get("encoding") or "latin-1"
    try:
        return raw_bytes.decode(enc)
    except (UnicodeDecodeError, LookupError):
        return raw_bytes.decode("latin-1", errors="replace")


def _parse_blocks(text: str) -> dict:
    """
    Parse the .dpo block structure into a dict.

    Key=Value blocks → dict of {block_name: {key: value}}.
    List blocks      → dict of {block_name + "_list": parsed_rows}.

    Handles:
      LABELS / LABELS_OTHERS         → list of strings (one label per row)
      SENSORS                        → list of (x, y, z) float tuples
      IMPEDANCE_VALUES               → list of lists of floats (rows = sessions)
    """
    result = {}

    # --- Key=Value blocks (BLOCK START … BLOCK END) ---
    kv_pattern = re.compile(
        r"^(\w[\w\s]*?)\s+START\s*$(.+?)^\1\s+END\s*$",
        re.MULTILINE | re.DOTALL,
    )
    for m in kv_pattern.finditer(text):
        block_name = m.group(1).strip()
        body       = m.group(2)
        kv = {}
        for line in body.splitlines():
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, _, v = line.partition("=")
                kv[k.strip()] = v.strip()
        # Last writer wins for repeated block names (e.g. DEVICE_PARAMETERS appears twice)
        # We store the first (EEG group) and the *_OTHERS variant separately via their name
        if block_name not in result:
            result[block_name] = kv
        # keep DEVICE_PARAMETERS_OTHERS distinct even if the block regex fires again
        # (it has a different name so this branch is only for exact duplicates)

    # --- List blocks (BLOCK START_LIST … BLOCK END_LIST) ---
    list_pattern = re.compile(
        r"^(\w[\w\s]*?)\s+START_LIST.*?$(.+?)^\1\s+END_LIST\s*$",
        re.MULTILINE | re.DOTALL,
    )
    for m in list_pattern.finditer(text):
        block_name = m.group(1).strip()
        body       = m.group(2)
        rows = [ln.strip() for ln in body.splitlines()
                if ln.strip() and not ln.strip().startswith("#")]

        if block_name in ("LABELS", "LABELS_OTHERS"):
            result[block_name + "_list"] = rows

        elif block_name == "SENSORS":
            coords = []
            for row in rows:
                parts = row.split()
                if len(parts) >= 3:
                    try:
                        coords.append(tuple(float(p) for p in parts[:3]))
                    except ValueError:
                        pass
            result["SENSORS_list"] = coords

        elif block_name == "IMPEDANCE_VALUES":
            sessions = []
            for row in rows:
                try:
                    sessions.append([float(v) for v in row.split()])
                except ValueError:
                    pass
            result["IMPEDANCE_VALUES_list"] = sessions

    return result


def _parse_start_datetime(params: dict) -> str:
    """Build an ISO-ish datetime string from DATA_PARAMETERS StartYear/Month/…"""
    try:
        year  = int(params.get("StartYear",  0) or 0)
        month = int(params.get("StartMonth", 0) or 0)
        day   = int(params.get("StartDay",   0) or 0)
        hour  = int(params.get("StartHour",  0) or 0)
        minute= int(params.get("StartMin",   0) or 0)
        sec   = int(params.get("StartSec",   0) or 0)
        ms    = int(params.get("StartMillisec", 0) or 0)
        if year and month and day:
            return f"{year:04d}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{sec:02d}.{ms:03d}"
    except (ValueError, TypeError):
        pass
    return ""


# ---------------------------------------------------------------------------
# Convenience: impedance summary per channel
# ---------------------------------------------------------------------------

def get_impedance_summary(hdr: dict) -> list[dict]:
    """
    Return per-channel impedance summary from a parsed header dict.

    Each entry: {'label': str, 'group': 'EEG'|'other',
                 'last_kOhm': float|None,   # last measured value (-1/-2 = sentinel)
                 'all_kOhm': list[float]}    # all sessions
    Sentinel values: -1 = not measured at that session, -2 = disabled/not connected.
    """
    labels = hdr["all_ch_labels"]
    n_eeg  = hdr["eeg_group_size"]
    sessions = hdr["impedances"]  # list of rows; each row = one session, n_channels values

    # Transpose: list-of-sessions → list-of-channels
    if sessions:
        per_channel = list(zip(*sessions))  # tuple per channel across sessions
    else:
        per_channel = [() for _ in labels]

    summary = []
    for i, label in enumerate(labels):
        group = "EEG" if i < n_eeg else "other"
        all_vals = list(per_channel[i]) if i < len(per_channel) else []
        # last non-sentinel value; raw file unit is Ohms — convert to kOhm for usability
        # (sentinel: -1 = not measured at that session, -2 = disabled/not connected)
        measured = [v for v in all_vals if v not in (-1.0, -2.0)]
        last_ohm = measured[-1] if measured else None
        last_kohm = (last_ohm / 1000.0) if last_ohm is not None else None
        all_kohm = [(v / 1000.0 if v not in (-1.0, -2.0) else v) for v in all_vals]
        summary.append({
            "label":      label,
            "group":      group,
            "last_kOhm":  last_kohm,
            "all_kOhm":   all_kohm,
        })
    return summary
