"""
Generate 2_select&remap_channels_curry_voila.ipynb from the EDF original.
Works at the parsed-JSON level (join source lines, replace, re-split).

Run from Inspect_EDF root:
    & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" tools_curry/_make_tool2_curry.py
"""
import json, os, sys

SRC = "tools/2_select&remap_channels_edf_voila.ipynb"
DST = "tools_curry/2_select&remap_channels_curry_voila.ipynb"

with open(SRC, encoding="utf-8") as f:
    nb = json.load(f)

errors = []

def get_src(cell):
    """Return cell source as a single joined string."""
    s = cell["source"]
    return "".join(s) if isinstance(s, list) else s

def set_src(cell, text):
    """Store source back, preserving original type (list vs string)."""
    orig = cell["source"]
    if isinstance(orig, list):
        cell["source"] = text.splitlines(keepends=True)
    else:
        cell["source"] = text

def replace_in_cell(cell, old, new, label):
    src = get_src(cell)
    if old not in src:
        errors.append(label)
        print(f"  ⚠ NOT FOUND: {label}")
        return False
    count = src.count(old)
    if count > 1:
        print(f"  ⚠ MULTIPLE ({count}) matches for: {label} — replacing all")
    set_src(cell, src.replace(old, new))
    print(f"  ✓ {label}")
    return True

def replace_all_cells(old, new, label):
    found = False
    for cell in nb["cells"]:
        src = get_src(cell)
        if old in src:
            set_src(cell, src.replace(old, new))
            found = True
    if not found:
        errors.append(label)
        print(f"  ⚠ NOT FOUND in any cell: {label}")
    else:
        print(f"  ✓ {label}")


def replace_between_anchors(cell, start_needle, end_needle, new_text, label):
    """Replace from the START of the line containing start_needle up to and including
    end_needle. Robust to internal whitespace/unicode (the EDF scan loop has \\u202f
    narrow no-break spaces in comments that break exact-match)."""
    src = get_src(cell)
    a = src.find(start_needle)
    if a == -1:
        errors.append(label); print(f"  ⚠ start anchor NOT FOUND: {label}"); return False
    line_start = src.rfind("\n", 0, a) + 1
    b = src.find(end_needle, a)
    if b == -1:
        errors.append(label); print(f"  ⚠ end anchor NOT FOUND: {label}"); return False
    set_src(cell, src[:line_start] + new_text + src[b + len(end_needle):])
    print(f"  ✓ {label}")
    return True

# Work on the code cell (there's typically 1 large code cell — find it)
code_cells = [c for c in nb["cells"] if c["cell_type"] == "code"]
main_cell = max(code_cells, key=lambda c: len(get_src(c)))

# ---------------------------------------------------------------------------
# 1. Add curry_header import (after 'import chardet')
# ---------------------------------------------------------------------------
replace_in_cell(main_cell,
    "    import chardet\n",
    "    import chardet\n"
    "    import sys as _sys\n"
    "    # curry shared modules — found whether Voila is launched from the repo root or tools_curry/\n"
    "    _here = os.getcwd()\n"
    "    for _cand in (_here, os.path.join(_here, 'tools_curry'),\n"
    "                  os.path.join(os.path.dirname(_here), 'tools_curry')):\n"
    "        if os.path.isfile(os.path.join(_cand, 'curry_header.py')):\n"
    "            _cand = os.path.abspath(_cand)\n"
    "            if _cand not in _sys.path:\n"
    "                _sys.path.insert(0, _cand)\n"
    "            break\n"
    "    from curry_header import read_curry_header\n",
    "import curry_header"
)

# ---------------------------------------------------------------------------
# 2. File discovery: .edf → .cdt (3 locations)
# ---------------------------------------------------------------------------
replace_in_cell(main_cell,
    "stems = [f.stem for f in folder.rglob('*') if f.suffix.lower() == '.edf' and not f.name.startswith('._')]",
    "stems = [f.stem for f in folder.rglob('*') if f.suffix == '.cdt' and not f.name.startswith('._')]",
    "stems discovery"
)
replace_in_cell(main_cell,
    "        chooser.edf_files = [\n"
    "            f for f in Path(chooser.folder_path).rglob('*')\n"
    "            if f.suffix.lower() == '.edf' and not f.name.startswith('._')\n"
    "            ]\n"
    "        if not chooser.edf_files:\n"
    "            print(f\"⚠️ There is no .edf file in your folder\")\n"
    "        else:\n"
    "            print(f\"\\nThere is {len(chooser.edf_files)} .edf files in your folder!\")",
    "        chooser.edf_files = [\n"
    "            f for f in Path(chooser.folder_path).rglob('*')\n"
    "            if f.suffix == '.cdt' and not f.name.startswith('._')\n"
    "            ]\n"
    "        if not chooser.edf_files:\n"
    "            print(f\"⚠️ There is no .cdt file in your folder\")\n"
    "        else:\n"
    "            print(f\"\\nThere is {len(chooser.edf_files)} .cdt files in your folder!\")",
    "chooser.edf_files discovery"
)
# Section 6 test discovery
replace_in_cell(main_cell,
    "                    if f.suffix.lower() == '.edf' and not f.name.startswith('._')\n",
    "                    if f.suffix == '.cdt' and not f.name.startswith('._')\n",
    "Section 6 discovery"
)

# ---------------------------------------------------------------------------
# 3. Scan loop body: replace read_edf_header_custom block
# ---------------------------------------------------------------------------
OLD_SCAN = (
    "                    edf_header = read_edf_header_custom(edf_path) \n"
    "                    \n"
    "                    # get subject name (corresponding to file_name)\n"
    "                    sub_name = edf_path.stem\n"
    "                    \n"
    "                    # get subject group (from the parent folder because in the ICEBERG database subfolders were created per patient group)\n"
    "                    sub_folder = edf_path.parent.name # get the parent folder of the subject file (path)\n"
    "                    \n"
    "                    # create df from signal info\n"
    "                    df = pd.DataFrame(edf_header)\n"
    "                        \n"
    "                    # theoretical resolution (edf are 16bit files so the eeg signal can take 2^16 values within the dynamic range)\n"
    "                    df['res_theoretical'] = (abs(pd.to_numeric(df['physical_min']))+abs(pd.to_numeric(df['physical_max'])))/pow(2,16)\n"
    "                    # turn theoretical resolution to uV if dimension is mV (if no dimension, it is a mess)\n"
    "                    df.loc[df['dimension'].str.contains('mv', case=False, na=False), 'res_theoretical'] *= 1000\n"
    "                    \n"
    "                    # get filtering info in different columns\n"
    "                    df['lowpass']   = df['prefiltering'].apply(lambda x: extract_filter_value(x, 'LP'))\n"
    "                    df['highpass']  = df['prefiltering'].apply(lambda x: extract_filter_value(x, 'HP'))\n"
    "                    df['notch']  = df['prefiltering'].apply(lambda x: extract_filter_value(x, 'NOTCH'))\n"
    "                    \n"
    "                    # add subject info in the dataframe\n"
    "                    df['subject'] = sub_name\n"
    "                    df['sub_folder'] = sub_folder\n"
    "                    df['group'] = np.nan # initialyze column 'group' with NaN\n"
    "                    # get group from participants table if any (else group will be inferred from subfolder or filename extension later)\n"
    "                    if found_group:\n"
    "                        df['group'] = subj_table.loc[subj_table['participant_id'] == sub_name, 'group'].iloc[0]\n"
    "        \n"
    "                    # extract filename component before and after subject number (so we assume subject name contains at least incrementing numbers that are at the beginning of the file name)  \n"
    "                    #   ^       → start of string  \n"
    "                    # (.*?)     → group 1: as few chars as possible, up to the first digit  \n"
    "                    # (\\d+)     → group 2: the number itself  \n"
    "                    # (.*)      → group 3: the rest of the string  \n"
    "                    # $         → end of string\n"
    "                    pre_comp = sub_num = post_comp = np.nan\n"
    "                    pattern = re.compile(r'^(.*?)(\\d+)(.*)$')\n"
    "                    m = pattern.match(sub_name)\n"
    "                    if m:\n"
    "                        pre_comp = m.group(1) or np.nan\n"
    "                        sub_num = m.group(2) or np.nan\n"
    "                        post_comp = m.group(3) or np.nan\n"
    "                    df['pre_fn_comp'] = pre_comp\n"
    "                    df['post_fn_comp'] = post_comp\n"
    "                    df['sub_num'] = sub_num\n"
    "                    \n"
    "                    df['path'] = str(edf_path)\n"
    "                    df['session'] = np.nan # session will be inferred later from file name component\n"
    "                    \n"
    "                    # select only the columns of interest\n"
    "                    df = df[['subject', 'group', 'session', 'path', 'sub_folder', 'sub_num', 'pre_fn_comp', 'post_fn_comp', 'channel', 'transducer_type', 'dimension', 'sampling_frequency', \n"
    "                         'highpass', 'lowpass', 'notch', 'physical_min', 'physical_max', 'res_theoretical']]\n"
    "                    \n"
    "                    # store subject data\n"
    "                    df_list.append(df)"
)
NEW_SCAN = (
    "                    hdr = read_curry_header(str(edf_path))\n"
    "                    sub_name = hdr['file_id']\n"
    "                    sub_folder = edf_path.parent.name\n"
    "                    # Build one row per channel using group info from the .dpo header\n"
    "                    rows = []\n"
    "                    for _i, _label in enumerate(hdr['all_ch_labels']):\n"
    "                        _grp = 'EEG' if _i < hdr['eeg_group_size'] else 'other'\n"
    "                        rows.append({'channel': _label, 'channel_group': _grp,\n"
    "                                     'sfreq': hdr['sfreq'], 'data_unit': hdr['data_unit']})\n"
    "                    df = pd.DataFrame(rows)\n"
    "                    df['subject'] = sub_name\n"
    "                    df['sub_folder'] = sub_folder\n"
    "                    df['group'] = np.nan\n"
    "                    if found_group:\n"
    "                        df['group'] = subj_table.loc[subj_table['participant_id'] == sub_name, 'group'].iloc[0]\n"
    "                    pre_comp = sub_num = post_comp = np.nan\n"
    "                    _m = re.compile(r'^(.*?)(\\d+)(.*)$').match(sub_name)\n"
    "                    if _m:\n"
    "                        pre_comp = _m.group(1) or np.nan\n"
    "                        sub_num  = _m.group(2) or np.nan\n"
    "                        post_comp = _m.group(3) or np.nan\n"
    "                    df['pre_fn_comp'] = pre_comp\n"
    "                    df['post_fn_comp'] = post_comp\n"
    "                    df['sub_num'] = sub_num\n"
    "                    df['path'] = str(edf_path)\n"
    "                    df['session'] = np.nan\n"
    "                    df_list.append(df)"
)
# Anchor-based (robust to the   narrow no-break spaces in the EDF comments,
# which broke the exact-match OLD_SCAN and silently reverted this on regeneration).
replace_between_anchors(main_cell,
    "edf_header = read_edf_header_custom(edf_path)",
    "df_list.append(df)",
    NEW_SCAN, "scan loop body")

# ---------------------------------------------------------------------------
# 4. Channel mask → Curry group-based
# ---------------------------------------------------------------------------
OLD_MASK = (
    "        # select only EEG and EOGs channels and return a warning if the number of participant is smaller/higher\n"
    "        mask = (\n"
    "            df_full['transducer_type'].str.contains(\n"
    "                r'\\bEEG\\b|\\bAGAGCL ELECTRODE\\b|\\bEOG\\b', case=False, na=False\n"
    "            )\n"
    "            | df_full['channel'].str.contains(r'EOG', case=False, na=False)\n"
    "            | df_full['channel'].str.contains(KNOWN_EEG_CHANNEL_RE, na=False)\n"
    "        )\n"
    "        df_ch = df_full[mask]\n"
    "        # remove the emg/ecg channels that were captured with the AGAGCL ELECTRODE transducer type \n"
    "        df_ch = df_ch[~df_ch['channel'].str.contains(r'emg|ecg', case=False, na=False)] # the ~ allows to not select the selection (like ! in matlab)"
)
NEW_MASK = (
    "        # For Curry, the EEG group is authoritative from the .dpo LABELS block\n"
    "        df_ch = df_full[df_full['channel_group'] == 'EEG'].copy()"
)
replace_in_cell(main_cell, OLD_MASK, NEW_MASK, "channel mask")

# ---------------------------------------------------------------------------
# 5. Section 6 signal load
# ---------------------------------------------------------------------------
OLD_LOAD = (
    "                                raw = mne.io.read_raw_edf(edf_path, preload=True, include=selected_channels) # we need to preload to use  re-ref; old param from Thomas encoding=\"latin-1\", \n"
    "                            except Exception as e:\n"
    "                                err = f\"\\t❌ Unexpected problem in loading edf file with mne: {e}\\n\""
)
NEW_LOAD = (
    "                                raw = mne.io.read_raw_curry(str(edf_path), preload=False, verbose='ERROR')\n"
    "                                _present = [ch for ch in selected_channels if ch in raw.ch_names]\n"
    "                                raw.pick(_present)\n"
    "                                raw.load_data()\n"
    "                                raw.rename_channels({k: v for k, v in sub_config['remap'].items() if k in raw.ch_names})\n"
    "                            except Exception as e:\n"
    "                                err = f\"\\t❌ Unexpected problem in loading Curry file with mne: {e}\\n\""
)
replace_in_cell(main_cell, OLD_LOAD, NEW_LOAD, "Section 6 signal load")

# ---------------------------------------------------------------------------
# 6. failed_edf_read.tsv → failed_cdt_read.tsv (all cells)
# ---------------------------------------------------------------------------
replace_all_cells("failed_edf_read.tsv", "failed_cdt_read.tsv", "failed_cdt_read.tsv")

# ---------------------------------------------------------------------------
# 7. Section 6 info text about loading
# ---------------------------------------------------------------------------
replace_in_cell(main_cell,
    "display(HTML(f\"<p>- <code>raw = mne.io.read_raw_edf(edf_path, preload=True, include=list(sub_config['remap'].keys()))</code></p>\"))",
    "display(HTML(f\"<p>- <code>raw = mne.io.read_raw_curry(cdt_path, preload=False)</code> then <code>raw.pick(selected_channels)</code></p>\"))",
    "Section 6 info text"
)

# ---------------------------------------------------------------------------
# 8. chooser title
# ---------------------------------------------------------------------------
replace_in_cell(main_cell,
    "Choose your study folder",
    "Choose your study folder (containing .cdt files)",
    "chooser title"
)

# ---------------------------------------------------------------------------
# Validate and write
# ---------------------------------------------------------------------------
try:
    # Re-serialize and re-parse to confirm valid JSON
    out_str = json.dumps(nb, ensure_ascii=False, indent=1)
    json.loads(out_str)
    print(f"\nJSON valid. n_cells={len(nb['cells'])}")
except Exception as e:
    print(f"\n✗ JSON invalid after edits: {e}")
    sys.exit(1)

if errors:
    print(f"\n⚠ {len(errors)} pattern(s) not found: {errors}")

with open(DST, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write("\n")
print(f"Written: {DST}")
