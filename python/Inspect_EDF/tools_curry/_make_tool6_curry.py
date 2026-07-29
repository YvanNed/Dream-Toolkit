"""
_make_tool6_curry.py — Generate 6_preprocessing_curry_voila.ipynb from the EDF original.

Adaptations from the EDF version:
- File discovery: .edf -> .cdt (bare suffix)
- Signal loading: mne.io.read_raw_edf(include=...) -> mne.io.read_raw_curry + pick + load_data
- Remove MNE suffix-duplicate helpers (drop_suffix_duplicates, adapt_remap_dict_to_suffixes)
  — Curry has a single global sampling rate and no -0/-1 duplicate channels.
- Events: Compumedics CSV/XML companions -> Curry French text export (*_ScoredEvents_Export.txt),
  parsed via curry_io.load_events_curry using the .cdt.dpo recording-start datetime.
- Curry shared-module imports (curry_header, curry_io).

Everything else (rejection methods, heatmap, per-stage summaries, skip/merge, custom stages,
sidecar JSON read by tool 7) is format-agnostic and kept byte-for-byte.

Run from Inspect_EDF root:
    & "$env:LOCALAPPDATA\\miniforge3\\envs\\inspect_edf\\python.exe" tools_curry/_make_tool6_curry.py
"""

import json, os, sys

SRC = "tools/6_preprocessing_voila.ipynb"
DST = "tools_curry/6_preprocessing_curry_voila.ipynb"

with open(SRC, encoding="utf-8") as f:
    nb = json.load(f)

errors = []


def get_src(cell):
    s = cell["source"]
    return "".join(s) if isinstance(s, list) else s


def set_src(cell, text):
    orig = cell["source"]
    if isinstance(orig, list):
        cell["source"] = text.splitlines(keepends=True)
    else:
        cell["source"] = text


def replace_in_cell(cell, old, new, label):
    src = get_src(cell)
    if old not in src:
        errors.append(label)
        print(f"  X NOT FOUND: {label}")
        return False
    count = src.count(old)
    if count > 1:
        print(f"  ! MULTIPLE ({count}) matches: {label} — replacing all")
    set_src(cell, src.replace(old, new))
    print(f"  ok {label}")
    return True


def replace_all_cells(old, new, label, required=True):
    found = False
    for cell in nb["cells"]:
        src = get_src(cell)
        if old in src:
            set_src(cell, src.replace(old, new))
            found = True
    if not found and required:
        errors.append(label)
        print(f"  X NOT FOUND in any cell: {label}")
    else:
        print(f"  ok {label}")


cells = nb["cells"]
# Index map (from inspection): 0=md, 1=imports, 2=shared fns, 3=md, 4=paths UI,
# 5=md, 6=params UI, 7=md, 8=run
md0     = cells[0]
imports = cells[1]
shared  = cells[2]
ui_paths = cells[4]
ui_param = cells[6]
run_cell = cells[8]

# ===========================================================================
print("=== Cell 0: markdown intro ===")
replace_in_cell(md0, "# Preprocessing — Phase 2",
                "# Preprocessing — Phase 2 — Curry 9 (.cdt)", "md title")
replace_in_cell(md0, "Select the **paths** (EDF folder,",
                "Select the **paths** (.cdt folder,", "md paths line")
replace_in_cell(md0, "[group sub-folders mirroring the EDF folder structure, if any]",
                "[group sub-folders mirroring the .cdt folder structure, if any]", "md output structure")

# ===========================================================================
print("\n=== Cell 1: imports ===")
# Keep the xml.etree.ElementTree import: _events_df_from_xml now passes through unchanged, so the
# Curry twin can also read a *.cdt.XML companion if one exists (mirrors the full TXT/CSV/XML chain).
# NB: the EDF tool-6 import block is wrapped in a try/except ImportError, so these lines are
# indented 4 spaces; the injected curry imports are kept inside the try so they are guarded too.
replace_in_cell(imports,
    "    from specparam import SpectralModel\nexcept ImportError as e:",
    "    from specparam import SpectralModel\n"
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
    "    from curry_header import read_curry_header\n"
    "    from curry_io import rec_start_from_header\n"
    "except ImportError as e:",
    "curry imports")

# ===========================================================================
print("\n=== Cell 2: shared functions ===")

# 2a. Remove MNE suffix-duplicate helpers (both functions + trailing blanks)
replace_in_cell(shared,
    "def drop_suffix_duplicates(raw):\n"
    "    \"\"\"Keep only the -0 variant when MNE creates -0/-1 duplicates for repeated channel names.\"\"\"\n"
    "    groups = {}\n"
    "    for ch in raw.ch_names:\n"
    "        if ch.endswith('-0') or ch.endswith('-1'):\n"
    "            base = ch.rsplit('-', 1)[0]\n"
    "            groups.setdefault(base, []).append(ch)\n"
    "    to_drop = []\n"
    "    for base, ch_list in groups.items():\n"
    "        if any(c.endswith('-0') for c in ch_list):\n"
    "            to_drop.extend(c for c in ch_list if not c.endswith('-0'))\n"
    "    if to_drop:\n"
    "        raw.drop_channels(to_drop)\n"
    "    return raw, to_drop\n"
    "\n\n"
    "def adapt_remap_dict_to_suffixes(raw, remap_dict):\n"
    "    \"\"\"Handle MNE's -0 suffix when the remap config uses the base channel name.\"\"\"\n"
    "    ch_set = set(raw.ch_names)\n"
    "    new_remap = {}\n"
    "    for base_label, target in remap_dict.items():\n"
    "        if base_label in ch_set:\n"
    "            new_remap[base_label] = target\n"
    "        elif f'{base_label}-0' in ch_set:\n"
    "            new_remap[f'{base_label}-0'] = target\n"
    "    return new_remap\n"
    "\n\n",
    "",
    "remove suffix helpers")

# 2b. Swap ONLY the recording-start reader: EDF fixed-header offsets -> Curry .cdt header.
# The rest of the event section (event_companion_paths, _events_df_from_txt/csv/xml, and the
# load_events TXT->CSV->XML dispatcher) is format-agnostic and passes through unchanged, so the
# Curry twin mirrors the same TXT-first / CSV / XML fallback chain. The function name is kept so
# the passed-through load_events dispatcher (which calls read_edf_start_datetime) still resolves.
OLD_STARTDT = (
    "def read_edf_start_datetime(edf_path):\n"
    "    \"\"\"Read the EDF recording-start datetime from the fixed header (offset 168 = date\n"
    "    'dd.mm.yy', 176 = time 'hh.mm.ss'), applying the EDF 2-digit-year clipping\n"
    "    (00-84 -> 20xx, 85-99 -> 19xx). Header-only read; returns a datetime or None on failure.\"\"\"\n"
    "    try:\n"
    "        with open(edf_path, 'rb') as f:\n"
    "            f.seek(168)\n"
    "            date_str = f.read(8).decode('ascii', 'replace').strip()   # dd.mm.yy\n"
    "            time_str = f.read(8).decode('ascii', 'replace').strip()   # hh.mm.ss\n"
    "        dd, mm, yy = (int(x) for x in date_str.split('.'))\n"
    "        hh, mi, ss = (int(x) for x in time_str.split('.'))\n"
    "        year = 2000 + yy if yy <= 84 else 1900 + yy\n"
    "        return datetime.datetime(year, mm, dd, hh, mi, ss)\n"
    "    except Exception:\n"
    "        return None"
)
NEW_STARTDT = (
    "def read_edf_start_datetime(cdt_path):\n"
    "    \"\"\"Curry recording-start datetime, read from the .cdt header (curry_header) rather than the\n"
    "    EDF fixed header. Returns a datetime or None on failure. (Name kept for the shared\n"
    "    load_events dispatcher, which converts the .txt clock times to seconds using it.)\"\"\"\n"
    "    try:\n"
    "        hdr = read_curry_header(str(cdt_path))\n"
    "        return rec_start_from_header(hdr)\n"
    "    except Exception:\n"
    "        return None"
)
replace_in_cell(shared, OLD_STARTDT, NEW_STARTDT, "start-datetime reader")

# ===========================================================================
print("\n=== Cell 4: paths UI ===")

# 4a. First discovery line (in _detect_hypno_suffixes; note the double space after edf_files)
replace_in_cell(ui_paths,
    "edf_files  = [f for f in sorted(edf_folder.rglob('*')) if f.suffix.lower() == '.edf' and not f.name.startswith('._')]",
    "curry_files  = [f for f in sorted(curry_folder.rglob('*')) if f.suffix == '.cdt' and not f.name.startswith('._')]",
    "cell4 discovery #1")

# 4b. Second discovery line (in _update_existing_reports_info; single space)
replace_in_cell(ui_paths,
    "edf_files = [f for f in sorted(edf_folder.rglob('*')) if f.suffix.lower() == '.edf' and not f.name.startswith('._')]",
    "curry_files = [f for f in sorted(curry_folder.rglob('*')) if f.suffix == '.cdt' and not f.name.startswith('._')]",
    "cell4 discovery #2")

# 4c. Event detection: the EDF source now has TWO passthrough blocks (TXT export + CSV) plus two
# suffix widgets, both suitable for Curry as-is (TXT default _ScoredEvents_Export.txt, CSV default
# _event_xml.csv — some Curry datasets also ship a CSV). Only the outer loop variable is migrated
# (edf -> cdt) so each block scans the .cdt stems; the file scans (edf_folder.rglob) and file
# lists (edf_files) are handled by the global renames below.
replace_in_cell(ui_paths,
    "        for edf in edf_files:\n"
    "            for t in all_txt_evt:\n"
    "                if os.path.normcase(t.name).startswith(os.path.normcase(edf.stem)):\n"
    "                    suf = t.name[len(edf.stem):]\n",
    "        for cdt in curry_files:\n"
    "            for t in all_txt_evt:\n"
    "                if os.path.normcase(t.name).startswith(os.path.normcase(cdt.stem)):\n"
    "                    suf = t.name[len(cdt.stem):]\n",
    "cell4 TXT detection loop var")
replace_in_cell(ui_paths,
    "        for edf in edf_files:\n"
    "            for c in all_csv:\n"
    "                if os.path.normcase(c.name).startswith(os.path.normcase(edf.stem)):\n"
    "                    suf = c.name[len(edf.stem):]\n",
    "        for cdt in curry_files:\n"
    "            for c in all_csv:\n"
    "                if os.path.normcase(c.name).startswith(os.path.normcase(cdt.stem)):\n"
    "                    suf = c.name[len(cdt.stem):]\n",
    "cell4 CSV detection loop var")
# Both event-detection else branches say "next to the EDF files" -> ".cdt files".
replace_all_cells("next to the EDF files", "next to the .cdt files",
                  "event detection EDF-files text", required=False)

# 4d. Hypno suffix detection: migrate loop var names only.
# The event-export exclusion is NOT done here any more: the EDF source now excludes 'event'
# suffixes from the AUTO-SELECTION itself (and keeps every .txt in the displayed list), so it
# passes through to the twin. Filtering all_txt here would additionally hide the event export
# from the displayed candidate list, which we want visible. See SPEC, Hypnogram-suffix
# auto-detection.
replace_in_cell(ui_paths,
    "        all_txt = [f for f in edf_folder.rglob('*') if f.suffix.lower() == '.txt']\n"
    "        suffix_counts = {}\n"
    "        for edf in edf_files:\n"
    "            for txt in all_txt:\n"
    "                if os.path.normcase(txt.name).startswith(os.path.normcase(edf.stem)):\n"
    "                    suffix = txt.name[len(edf.stem):]\n"
    "                    suffix_counts[suffix] = suffix_counts.get(suffix, 0) + 1\n",
    "        all_txt = [f for f in curry_folder.rglob('*') if f.suffix.lower() == '.txt']\n"
    "        suffix_counts = {}\n"
    "        for cdt in curry_files:\n"
    "            for txt in all_txt:\n"
    "                if os.path.normcase(txt.name).startswith(os.path.normcase(cdt.stem)):\n"
    "                    suffix = txt.name[len(cdt.stem):]\n"
    "                    suffix_counts[suffix] = suffix_counts.get(suffix, 0) + 1\n",
    "cell4 hypno detection")

# 4d2. Docstring mentioning "EDF folder"
replace_in_cell(ui_paths,
    '"""Auto-detect hypnogram .txt suffixes in the EDF folder and populate the widget."""',
    '"""Auto-detect hypnogram .txt suffixes in the .cdt folder and populate the widget."""',
    "cell4 docstring")

# 4e. fc_config title (from select&remap_channels_edf -> curry)
replace_in_cell(ui_paths,
    "fc_config.title = '<b>remap_reref_persubject.json</b> (from select&amp;remap_channels_edf) :'",
    "fc_config.title = '<b>remap_reref_persubject.json</b> (from select&amp;remap_channels_curry) :'",
    "cell4 fc_config title")

# 4f. Suffix widgets pass through unchanged: the EDF 'Event TXT suffix' (_ScoredEvents_Export.txt)
# and 'Event CSV suffix' (_event_xml.csv) defaults both suit Curry, so the twin mirrors both fields.

# ===========================================================================
print("\n=== Cell 6: params UI (event counter) ===")

# 6a. Discovery line in _count_affected_epochs
replace_in_cell(ui_param,
    "edf_cand = [f for f in edf_folder.rglob('*')\n"
    "                        if os.path.normcase(f.stem) == os.path.normcase(fid) and f.suffix.lower() == '.edf']",
    "cdt_cand = [f for f in curry_folder.rglob('*')\n"
    "                        if os.path.normcase(f.stem) == os.path.normcase(fid) and f.suffix == '.cdt']",
    "cell6 discovery")

# 6b. "Select the EDF data folder" message
replace_in_cell(ui_param,
    "print('Select the EDF data folder in Section 1 first.')",
    "print('Select the .cdt data folder in Section 1 first.')",
    "cell6 EDF msg")

# ===========================================================================
print("\n=== Cell 8: run ===")

# 8a-pre. Ground-truth recording discovery in on_load_participants (.edf -> .cdt)
replace_in_cell(run_cell,
    "    rec_paths = [f for f in sorted(edf_folder.rglob('*'))\n"
    "                 if f.suffix.lower() == '.edf' and not f.name.startswith('._')]",
    "    rec_paths = [f for f in sorted(curry_folder.rglob('*'))\n"
    "                 if f.suffix == '.cdt' and not f.name.startswith('._')]",
    "cell8 on_load recording discovery")

# 8a. Discovery line
replace_in_cell(run_cell,
    "edf_candidates = [f for f in edf_folder.rglob('*') if os.path.normcase(f.stem) == os.path.normcase(file_id) and f.suffix.lower() == '.edf']",
    "cdt_candidates = [f for f in curry_folder.rglob('*') if os.path.normcase(f.stem) == os.path.normcase(file_id) and f.suffix == '.cdt']",
    "cell8 discovery")

# 8b. Load + rename block
OLD_LOAD = (
    "            # [load] preload=False : on ne lit que l'en-tête pour l'instant.\n"
    "            # include= utilise les noms d'ORIGINE du JSON (remap.keys()). Passer include= à la\n"
    "            # LECTURE (et non un pick paresseux après coup) est important : (a) ça exclut\n"
    "            # d'emblée les canaux hors-montage à fréquence plus élevée (ex. ECG 512 Hz), ce qui\n"
    "            # évite un AssertionError du lecteur EDF de MNE sur la lecture partielle ET préserve\n"
    "            # la fréquence native de l'EEG ; (b) ça évite tout dictionnaire de remap inversé.\n"
    "            try:\n"
    "                raw = mne.io.read_raw_edf(\n"
    "                    str(edf_path), preload=False, encoding='latin-1',\n"
    "                    include=list(sub_config.get('remap', {}).keys()), verbose=False\n"
    "                )\n"
    "            except Exception as e:\n"
    "                with out_run:\n"
    "                    print(f'[{file_id}] Error loading EDF: {e}')\n"
    "                failed.append({'file_id': file_id, 'reason': f'EDF loading: {e}'})\n"
    "                progress.value = idx + 1\n"
    "                continue\n"
    "\n"
    "            raw, _ = drop_suffix_duplicates(raw)\n"
    "\n"
    "            # [B] Renommage canaux origine -> remappé — non-fatal en soi, mais la sélection\n"
    "            # qui suit en dépend (cf. garde-fou 'present' juste après).\n"
    "            try:\n"
    "                remap_adapted = adapt_remap_dict_to_suffixes(raw, sub_config['remap'])\n"
    "                raw.rename_channels(remap_adapted)\n"
    "            except Exception as e:\n"
    "                with out_run:\n"
    "                    print(f'[{file_id}] ⚠ Error renaming channels: {e} — channels not renamed.')\n"
    "\n"
    "            # selected_channels porte les noms *remappés* (UI / quality_summary), même\n"
    "            # namespace que raw après le renommage. On retire les canaux désélectionnés AVANT\n"
    "            # load_data() pour ne lire sur disque que les canaux finalement gardés.\n"
    "            present = [ch for ch in selected_channels if ch in raw.ch_names]\n"
    "            if not present:\n"
    "                with out_run:\n"
    "                    print(f'[{file_id}] No selected channels found after renaming — skipped.')\n"
    "                failed.append({'file_id': file_id, 'reason': 'no channels after rename'})\n"
    "                progress.value = idx + 1\n"
    "                continue\n"
    "            to_drop = [ch for ch in raw.ch_names if ch not in present]\n"
    "            if to_drop:\n"
    "                raw.drop_channels(to_drop)\n"
    "            try:\n"
    "                raw.load_data()   # ne lit sur disque que les canaux gardés\n"
    "            except Exception as e:\n"
    "                with out_run:\n"
    "                    print(f'[{file_id}] Error reading EDF data: {e}')\n"
    "                failed.append({'file_id': file_id, 'reason': f'EDF data loading: {e}'})\n"
    "                progress.value = idx + 1\n"
    "                continue\n"
)
NEW_LOAD = (
    "            # [load] Read the Curry header lazily (preload=False). All Curry channels share a\n"
    "            # single sampling rate, so the EDF include=-at-read trick is unnecessary: we pick the\n"
    "            # montage channels (original names from remap.keys()) then load only those from disk.\n"
    "            try:\n"
    "                raw = mne.io.read_raw_curry(str(cdt_path), preload=False, verbose='ERROR')\n"
    "                _montage = [ch for ch in sub_config.get('remap', {}).keys() if ch in raw.ch_names]\n"
    "                raw.pick(_montage)\n"
    "            except Exception as e:\n"
    "                with out_run:\n"
    "                    print(f'[{file_id}] Error loading Curry file: {e}')\n"
    "                failed.append({'file_id': file_id, 'reason': f'Curry loading: {e}'})\n"
    "                progress.value = idx + 1\n"
    "                continue\n"
    "\n"
    "            # [B] Rename original channel names -> remapped names — non-fatal in itself, but the\n"
    "            # selection just below depends on it (cf. the 'present' guard).\n"
    "            try:\n"
    "                raw.rename_channels({k: v for k, v in sub_config['remap'].items() if k in raw.ch_names})\n"
    "            except Exception as e:\n"
    "                with out_run:\n"
    "                    print(f'[{file_id}] ⚠ Error renaming channels: {e} — channels not renamed.')\n"
    "\n"
    "            # selected_channels carries the *remapped* names (UI / quality_summary), same\n"
    "            # namespace as raw after renaming. We drop the deselected channels BEFORE\n"
    "            # load_data() so only the kept channels are read from disk.\n"
    "            present = [ch for ch in selected_channels if ch in raw.ch_names]\n"
    "            if not present:\n"
    "                with out_run:\n"
    "                    print(f'[{file_id}] No selected channels found after renaming — skipped.')\n"
    "                failed.append({'file_id': file_id, 'reason': 'no channels after rename'})\n"
    "                progress.value = idx + 1\n"
    "                continue\n"
    "            to_drop = [ch for ch in raw.ch_names if ch not in present]\n"
    "            if to_drop:\n"
    "                raw.drop_channels(to_drop)\n"
    "            try:\n"
    "                raw.load_data()   # reads only the kept channels from disk\n"
    "            except Exception as e:\n"
    "                with out_run:\n"
    "                    print(f'[{file_id}] Error reading Curry data: {e}')\n"
    "                failed.append({'file_id': file_id, 'reason': f'Curry data loading: {e}'})\n"
    "                progress.value = idx + 1\n"
    "                continue\n"
)
replace_in_cell(run_cell, OLD_LOAD, NEW_LOAD, "cell8 load block")

# 8b2. Validation ERROR message listing required paths (user-facing)
replace_in_cell(run_cell,
    "print('ERROR: Please select the required paths (EDF folder, JSON config, output). '",
    "print('ERROR: Please select the required paths (.cdt folder, JSON config, output). '",
    "cell8 validation msg")

# 8c. "EDF not found" -> ".cdt not found" (print + reason)
replace_in_cell(run_cell,
    "print(f'[{file_id}] EDF not found in {edf_folder}')",
    "print(f'[{file_id}] .cdt not found in {curry_folder}')",
    "cell8 not-found print")
replace_in_cell(run_cell,
    "failed.append({'file_id': file_id, 'reason': 'EDF not found'})",
    "failed.append({'file_id': file_id, 'reason': '.cdt not found'})",
    "cell8 not-found reason")

# ===========================================================================
# Per-channel rejection (memory) + progress feedback: now implemented in the EDF source
# (cell 2 per-channel compute_rejection_masks; cell 8 del raw / drop epochs_data_uV /
# reject call + _rej_progress callback). It is format-agnostic and passes through into the
# Curry twin unchanged — no transformation needed here. See tools/6_preprocessing_voila.ipynb.

# ===========================================================================
print("\n=== Cell 8: [H] context-channels companion (Curry reader) ===")
# The [H] block reads the declared EOG/EMG/ECG context channels to persist a *_context-epo.fif.
# For Curry, swap the EDF reader for read_raw_curry and drop the (removed) suffix-dedup helper.
# Run BEFORE the global edf_path -> cdt_path rename so these OLD strings still match.
replace_in_cell(run_cell,
    "                    # Read only the declared context channels (include=), exactly like the EEG read:\n"
    "                    # otherwise MNE upsamples every channel in the file to its max rate (e.g. a fast\n"
    "                    # ECG), which was raising 'bad allocation' on mixed-rate montages.\n"
    "                    ctx_probe = mne.io.read_raw_edf(str(edf_path), preload=False, encoding='latin-1',\n"
    "                                                    include=list(orig_to_role.keys()), verbose=False)\n"
    "                    ctx_probe, _ = drop_suffix_duplicates(ctx_probe)\n",
    "                    ctx_probe = mne.io.read_raw_curry(str(edf_path), preload=False, verbose='ERROR')\n",
    "cell8 [H] context reader")

# ===========================================================================
print("\n=== Global identifier renames (ordered, longest-first) ===")
# Order matters: replace the longer identifiers before their prefixes.
for old, new in [
    ("edf_candidates", "cdt_candidates"),   # before edf_cand
    ("edf_rel_str",    "cdt_rel_str"),      # before edf_rel
    ("edf_folder",     "curry_folder"),
    ("edf_files",      "curry_files"),
    ("edf_cand",       "cdt_cand"),
    ("edf_path",       "cdt_path"),
    ("edf_rel",        "cdt_rel"),
    ("fc_edf",         "fc_curry"),
]:
    replace_all_cells(old, new, f"{old} -> {new}", required=False)

# Remaining user-facing "No EDF files found" text (2 spots in cell 4)
replace_all_cells("No EDF files found in selected folder",
                  "No .cdt files found in selected folder",
                  "No EDF files text", required=False)

# ===========================================================================
print("\n=== Validation ===")
try:
    out_str = json.dumps(nb, ensure_ascii=False, indent=1)
    json.loads(out_str)
    print(f"JSON valid. n_cells={len(nb['cells'])}")
except Exception as e:
    print(f"X JSON invalid after edits: {e}")
    sys.exit(1)

# Compile every code cell as a syntax gate
import ast
for i, cell in enumerate(nb["cells"]):
    if cell["cell_type"] == "code":
        src = get_src(cell)
        try:
            ast.parse(src)
        except SyntaxError as e:
            print(f"X SyntaxError in code cell {i}: {e}")
            errors.append(f"syntax cell {i}")

if errors:
    print(f"\n! {len(errors)} issue(s): {errors}")

with open(DST, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write("\n")
print(f"Written: {DST}")
