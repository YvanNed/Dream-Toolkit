"""
_make_tool5_curry.py — Generate 5_quality_overview_curry_voila.ipynb from the EDF original.

Changes from EDF version:
- File discovery: .edf → .cdt (bare suffix)
- Signal loading: mne.io.read_raw_edf → mne.io.read_raw_curry with pick + load_data
- No EDF physical bounds: remove bounds_pct metric, get_phys_bounds_uV, thresh_bounds widget
- No MNE suffix-duplicate helpers: remove drop_suffix_duplicates, adapt_remap_dict_to_suffixes
- Curry shared-module imports added (curry_header, curry_io)
- Time-series / histogram amplitude axes capped at a wide physiological ceiling
  (DISPLAY_YLIM_UV = 500 uV): DC-coupled data has no export clipping, so drift can
  blow up the p99.9 autoscale and crush the real EEG

Run from Inspect_EDF root:
    & "$env:LOCALAPPDATA\\miniforge3\\envs\\inspect_edf\\python.exe" tools_curry/_make_tool5_curry.py
"""

import json, os, sys

SRC = "tools/5_quality_overview_voila.ipynb"
DST = "tools_curry/5_quality_overview_curry_voila.ipynb"

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


def replace_in_cell(cell, old, new, label, allow_missing=False):
    src = get_src(cell)
    if old not in src:
        if not allow_missing:
            errors.append(label)
            print(f"  ⚠ NOT FOUND: {label}")
        return False
    count = src.count(old)
    if count > 1:
        print(f"  ⚠ MULTIPLE ({count}) matches: {label} — replacing all")
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


code_cells = [c for c in nb["cells"] if c["cell_type"] == "code"]
cell0 = code_cells[0]  # imports
cell1 = code_cells[1]  # shared functions
cell3 = code_cells[2]  # UI
cell4 = code_cells[3]  # main processing

print("=== Cell 0: imports ===")
replace_in_cell(cell0,
    "import yasa",
    "import yasa\nimport sys as _sys\n"
    "# curry shared modules — located next to this notebook\n"
    "_here = os.path.dirname(os.path.abspath('__file__'))\n"
    "if _here not in _sys.path:\n"
    "    _sys.path.insert(0, _here)\n"
    "from curry_header import read_curry_header\n"
    "from curry_io import read_curry_signal, load_hypnogram_curry",
    "cell0 curry imports",
)

print("\n=== Cell 1: shared functions ===")

# 1a. Remove drop_suffix_duplicates
replace_in_cell(cell1,
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
    "\n\n",
    "",
    "remove drop_suffix_duplicates",
)

# 1b. Remove adapt_remap_dict_to_suffixes
replace_in_cell(cell1,
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
    "remove adapt_remap_dict_to_suffixes",
)

# 1c. Remove get_phys_bounds_uV
replace_in_cell(cell1,
    "def get_phys_bounds_uV(extras, ch_idx):\n"
    "    \"\"\"\n"
    "    Reconstruct physical_min/max (in uV) from MNE _raw_extras.\n"
    "    MNE 1.9 stores physical_max but not physical_min explicitly, and it keeps\n"
    "    physical_max / offset in the channel's *native* EDF physical unit (uV, mV,\n"
    "    V...), not in uV. extras['units'] is the multiplier from that native unit to\n"
    "    volts (1e-6 for uV, 1e-3 for mV, 1.0 for V), so multiplying by units * 1e6\n"
    "    brings the bounds into uV, matching raw.get_data() * 1e6. Without this, a\n"
    "    channel declared in mV (e.g. Compumedics EOG/EMG/ECG, physical_max=1.0 mV)\n"
    "    is compared against a 1.0 uV bound -- 1000x too small -- flagging ~100% of\n"
    "    samples as 'at EDF bounds'.\n"
    "    For symmetric 16-bit EDF: physical_min = -physical_max + 2*offset.\n"
    "    Returns (phys_min_uV, phys_max_uV).\n"
    "    \"\"\"\n"
    "    to_uV = float(extras['units'][ch_idx]) * 1e6\n"
    "    phys_max = float(extras['physical_max'][ch_idx]) * to_uV\n"
    "    offset = float(extras['offsets'][ch_idx]) * to_uV\n"
    "    return -phys_max + 2 * offset, phys_max\n"
    "\n\n",
    "",
    "remove get_phys_bounds_uV",
)

# 1d. compute_signal_metrics: remove phys bounds parameters
replace_in_cell(cell1,
    "def compute_signal_metrics(sig_uV, phys_min_uV, phys_max_uV):\n"
    "    \"\"\"Compute all quality metrics for one EEG channel (unfiltered signal, in uV).\"\"\"",
    "def compute_signal_metrics(sig_uV):\n"
    "    \"\"\"Compute all quality metrics for one EEG channel (unfiltered signal, in uV).\"\"\"",
    "compute_signal_metrics signature",
)

# 1e. Remove bounds_pct computation block
replace_in_cell(cell1,
    "    # EDF physical bounds: fraction of samples at/near the header declared limits.\n"
    "    # Detects saturation at the EDF dynamic range boundary (hard clipping).\n"
    "    if not (np.isnan(phys_min_uV) or np.isnan(phys_max_uV)):\n"
    "        bounds_pct = float(\n"
    "            ((sig_uV <= phys_min_uV + 0.5) | (sig_uV >= phys_max_uV - 0.5)).mean()\n"
    "        ) * 100\n"
    "    else:\n"
    "        bounds_pct = 0.0\n"
    "\n"
    "    # Histogram + Savitzky-Golay",
    "    # Histogram + Savitzky-Golay",
    "remove bounds_pct computation",
)

# 1f. Remove bounds_pct from return dict
replace_in_cell(cell1,
    "        'bounds_pct': bounds_pct,\n",
    "",
    "remove bounds_pct from return dict",
)

# 1g. Remove bounds_pct check from flag_channel
replace_in_cell(cell1,
    "    if metrics['bounds_pct'] > thresholds['bounds_pct']:\n"
    "        reasons.append(f\"bounds_pct={metrics['bounds_pct']:.2f}% > {thresholds['bounds_pct']}%\")\n",
    "",
    "remove bounds_pct flag_channel check",
)

# 1h. Remove bounds_pct from OVERVIEW_KEY_METRICS (shared by generate_dataset_overview and the
# per-participant "All electrodes" figures; 4-space indent — the literal KEY_METRICS/LABELS blocks
# were lifted to these module-level constants).
replace_in_cell(cell1,
    "    \"std_uV\", \"flat_pct\", \"bounds_pct\", \"hist_extreme_pct\",",
    "    \"std_uV\", \"flat_pct\", \"hist_extreme_pct\",",
    "OVERVIEW_KEY_METRICS remove bounds_pct",
)

# 1i. Remove bounds_pct from ALL_NUMERIC
replace_in_cell(cell1,
    "        \"flat_pct\", \"bounds_pct\", \"hist_extreme_pct\",",
    "        \"flat_pct\", \"hist_extreme_pct\",",
    "ALL_NUMERIC remove bounds_pct",
)

# 1j. Remove bounds_pct from OVERVIEW_METRIC_LABELS (4-space indent — lifted to a module constant)
replace_in_cell(cell1,
    "    \"bounds_pct\": \"At EDF bounds (%)\",\n",
    "",
    "OVERVIEW_METRIC_LABELS remove bounds_pct",
)

# 1k. Remove bounds_pct from STAGE_METRICS (inside generate_dataset_overview)
replace_in_cell(cell1,
    "STAGE_METRICS = [\"mean_uV\", \"std_uV\", \"flat_pct\", \"bounds_pct\", \"hist_extreme_pct\", \"p99_abs_uV\", \"p999_abs_uV\"]",
    "STAGE_METRICS = [\"mean_uV\", \"std_uV\", \"flat_pct\", \"hist_extreme_pct\", \"p99_abs_uV\", \"p999_abs_uV\"]",
    "STAGE_METRICS remove bounds_pct",
)

# 1l. Update comment referencing bounds_pct in flag_channel
replace_in_cell(cell1,
    "    # the declared EDF physical range, undetectable via bounds_pct alone).",
    "    # the declared EDF physical range.",
    "comment bounds_pct alone",
)

print("\n=== Cell 3: UI ===")

# 3a. File discovery: .edf → .cdt
replace_in_cell(cell3,
    "edf_files = [f for f in sorted(data_folder.rglob('*')) if f.suffix.lower() == '.edf' and not f.name.startswith('._')]\n"
    "    n_total = len(edf_files)\n"
    "    if n_total == 0:\n"
    "        existing_reports_info.value = '<small style=\"color:#888;\">No EDF files found in selected folder (recursive scan).</small>'\n",
    "cdt_files = [f for f in sorted(data_folder.rglob('*')) if f.suffix == '.cdt' and not f.name.startswith('._')]\n"
    "    n_total = len(cdt_files)\n"
    "    if n_total == 0:\n"
    "        existing_reports_info.value = '<small style=\"color:#888;\">No .cdt files found in selected folder (recursive scan).</small>'\n",
    "cell3 cdt discovery",
)

# 3b. (removed) The n_existing report count is now a plain `for f in edf_files` loop that
# also checks the per-file _quality_metrics.tsv; the global edf_files → cdt_files rename below
# translates it, so no targeted replacement is needed here.

# 3c. Hypno suffix detection loop
replace_in_cell(cell3,
    "    for edf in edf_files:\n"
    "        for txt in all_txt:\n"
    "            if os.path.normcase(txt.name).startswith(os.path.normcase(edf.stem)):\n"
    "                suffix = txt.name[len(edf.stem):]\n"
    "                suffix_counts[suffix] = suffix_counts.get(suffix, 0) + 1",
    "    for cdt in cdt_files:\n"
    "        for txt in all_txt:\n"
    "            if os.path.normcase(txt.name).startswith(os.path.normcase(cdt.stem)):\n"
    "                suffix = txt.name[len(cdt.stem):]\n"
    "                suffix_counts[suffix] = suffix_counts.get(suffix, 0) + 1",
    "cell3 hypno suffix loop",
)

# 3d. EDF files matched text
replace_in_cell(cell3,
    "f'&nbsp;— {best_count}/{n_total} EDF files matched</small>'",
    "f'&nbsp;— {best_count}/{n_total} .cdt files matched</small>'",
    "cell3 EDF files matched text",
)

# 3e. fc_config title
replace_in_cell(cell3,
    "fc_config.title = '<b>Select remap/reref config JSON</b> (from select&amp;remap_channels_edf):'",
    "fc_config.title = '<b>Select remap/reref config JSON</b> (from select&amp;remap_channels_curry):'",
    "cell3 fc_config title",
)

# 3f. Remove thresh_bounds widget (between thresh_flat and thresh_peaks)
replace_in_cell(cell3,
    ")\nthresh_bounds = widgets.BoundedFloatText(\n"
    "    value=1.0, min=0.0, max=100.0, step=0.1,\n"
    "    description='bounds_pct (%) >',\n"
    "    style={'description_width': '160px'},\n"
    "    layout=widgets.Layout(width='300px')\n"
    ")\nthresh_peaks = widgets.BoundedIntText(",
    ")\nthresh_peaks = widgets.BoundedIntText(",
    "cell3 remove thresh_bounds widget",
)

# 3g. Remove thresh_row(thresh_bounds,...) from VBox
replace_in_cell(cell3,
    ",\n        thresh_row(thresh_bounds,\n"
    "                   'Fraction at EDF physical-range limits — detects hard saturation (declared range)'),\n"
    "        thresh_row(thresh_peaks,",
    ",\n        thresh_row(thresh_peaks,",
    "cell3 remove thresh_bounds VBox row",
)

print("\n=== Cell 4: main processing ===")

# 4a. edf_files discovery and progress setup
replace_in_cell(cell4,
    "edf_files = [f for f in sorted(data_folder.rglob('*')) if f.suffix.lower() == '.edf' and not f.name.startswith('._')]\n"
    "        if not edf_files:\n"
    "            with out:\n"
    "                print(f'No EDF files found in {data_folder}')\n"
    "            return\n"
    "\n"
    "        progress.max = len(edf_files)\n"
    "        progress.value = 0\n",
    "cdt_files = [f for f in sorted(data_folder.rglob('*')) if f.suffix == '.cdt' and not f.name.startswith('._')]\n"
    "        if not cdt_files:\n"
    "            with out:\n"
    "                print(f'No .cdt files found in {data_folder}')\n"
    "            return\n"
    "\n"
    "        progress.max = len(cdt_files)\n"
    "        progress.value = 0\n",
    "cell4 cdt discovery",
)

# 4b. Remove bounds_pct from thresholds dict
replace_in_cell(cell4,
    "thresholds = {\n"
    "            'flat_pct': thresh_flat.value,\n"
    "            'bounds_pct': thresh_bounds.value,\n"
    "            'n_peaks': thresh_peaks.value,\n",
    "thresholds = {\n"
    "            'flat_pct': thresh_flat.value,\n"
    "            'n_peaks': thresh_peaks.value,\n",
    "cell4 remove bounds_pct from thresholds",
)

# 4c. Remove 'At EDF bounds (%)' from interpretations.update
replace_in_cell(cell4,
    "            'At EDF bounds (%)': f'flag if&nbsp;&gt;&nbsp;{thresholds[\"bounds_pct\"]:g}%&nbsp;&mdash; threshold defined manually in the notebook ; typical threshold is ??',\n",
    "",
    "cell4 remove At EDF bounds interpretations",
)

# 4d. Replace EDF loading block + phys bounds + drop suffix + rename
replace_in_cell(cell4,
    "# --- Load EDF ---\n"
    "            try:\n"
    "                raw = mne.io.read_raw_edf(\n"
    "                    str(edf_path), preload=True, encoding='latin-1',\n"
    "                    include=selected_channels, verbose=False\n"
    "                )\n"
    "            except Exception as e:\n"
    "                with out:\n"
    "                    print(f'ERROR loading {file_id}: {e}')\n"
    "                failed.append({'file_id': file_id, 'reason': f'EDF loading: {e}'})\n"
    "                continue\n"
    "\n"
    "            # --- Save physical bounds before channel manipulation ---\n"
    "            extras = raw._raw_extras[0]\n"
    "            phys_bounds_by_name = {}\n"
    "            for idx, ch in enumerate(raw.ch_names):\n"
    "                base = ch.rsplit('-', 1)[0] if (ch.endswith('-0') or ch.endswith('-1')) else ch\n"
    "                phys_bounds_by_name[base] = get_phys_bounds_uV(extras, idx)\n"
    "\n"
    "            # --- Drop suffix duplicates, build remap, rename ---\n"
    "            raw, _ = drop_suffix_duplicates(raw)\n"
    "            remap_adapted = adapt_remap_dict_to_suffixes(raw, sub_config['remap'])\n"
    "            final_to_orig = {}\n"
    "            for orig, new in remap_adapted.items():\n"
    "                base = orig.rsplit('-', 1)[0] if (orig.endswith('-0') or orig.endswith('-1')) else orig\n"
    "                final_to_orig[new] = base\n"
    "            raw.rename_channels(remap_adapted)\n"
    "            sf = raw.info['sfreq']\n",
    "# --- Load Curry signal (pick selected channels before loading to limit memory use) ---\n"
    "            try:\n"
    "                raw = mne.io.read_raw_curry(str(edf_path), preload=False, verbose='ERROR')\n"
    "                _present = [ch for ch in selected_channels if ch in raw.ch_names]\n"
    "                raw.pick(_present)\n"
    "                raw.load_data()\n"
    "            except Exception as e:\n"
    "                with out:\n"
    "                    print(f'ERROR loading {file_id}: {e}')\n"
    "                failed.append({'file_id': file_id, 'reason': f'Curry loading: {e}'})\n"
    "                continue\n"
    "\n"
    "            raw.rename_channels({k: v for k, v in sub_config['remap'].items() if k in raw.ch_names})\n"
    "            sf = raw.info['sfreq']\n",
    "cell4 replace load block",
)

# --- Curry-only: the high-pass control/apply now lives natively in the EDF notebook (inherited
# here). DC-coupled Curry data benefits from it, so flip its default ON and adjust the label.
# The default high-pass value is also lowered (0.1 Hz) since Curry has no hardware high-pass. ---
replace_in_cell(cell3,
    "    value=False, description='Activate high-pass',\n",
    "    value=True, description='Activate high-pass (DC-coupled Curry data)',\n",
    "cell3 high-pass default (DC-coupled ON)")

replace_in_cell(cell3,
    "hp_freq = widgets.BoundedFloatText(\n"
    "    value=0.5, min=0.01, max=5.0, step=0.05,\n",
    "hp_freq = widgets.BoundedFloatText(\n"
    "    value=0.1, min=0.01, max=5.0, step=0.05,\n",
    "cell3 high-pass default freq (Curry)")

# 3h. Acquisition scan: the EDF byte parser (read_edf_sf_highpass) does not apply to .cdt, so read
# the Curry header via MNE instead. Curry is DC-coupled → no acquisition high-pass (none/DC).
# (read_edf_sf_highpass stays defined but unused in the Curry twin — harmless.)
replace_in_cell(cell3,
    "                _sf, _hp = read_edf_sf_highpass(_edf, _keep or None)\n",
    "                _r = mne.io.read_raw_curry(str(_edf), preload=False, verbose='ERROR')\n"
    "                _sf, _hp = f\"{_r.info['sfreq']:.0f} Hz\", 'none/DC'\n",
    "cell3 acq scan reader (Curry)")

replace_in_cell(cell3,
    "            '<small style=\"color:#555;\">Selected channels (EDF header): '\n",
    "            '<small style=\"color:#555;\">Selected channels (Curry header): '\n",
    "cell3 acq scan label (Curry)")

# 4e. Remove EDF-specific comment before selected_channels
replace_in_cell(cell4,
    "sub_config = config_dict[file_id]\n"
    "            # include= utilise les noms d'ORIGINE de l'EDF (les clés du remap), évalué à la\n"
    "            # LECTURE. On ne garde ainsi que les canaux du montage : la fréquence native de\n"
    "            # l'EEG est préservée (pas de suréchantillonnage vers un canal hors-montage plus\n"
    "            # rapide comme un ECG 512 Hz), et on évite un AssertionError de lecture partielle\n"
    "            # du lecteur EDF de MNE. Même motif que 6_preprocessing_voila.\n"
    "            selected_channels",
    "sub_config = config_dict[file_id]\n"
    "            selected_channels",
    "cell4 remove EDF comment",
)

# 4f. Remove phys_bounds lookup + update compute_signal_metrics call (main channel loop)
replace_in_cell(cell4,
    "                orig_base = final_to_orig.get(ch, ch)\n"
    "                phys_min, phys_max = phys_bounds_by_name.get(orig_base, (np.nan, np.nan))\n"
    "                m = compute_signal_metrics(sig_uV, phys_min, phys_max)\n",
    "                m = compute_signal_metrics(sig_uV)\n",
    "cell4 remove phys_bounds call",
)

# 4f2. Update compute_signal_metrics call in per-stage loop
replace_in_cell(cell4,
    "                        sm = compute_signal_metrics(stage_sig, phys_min, phys_max)\n",
    "                        sm = compute_signal_metrics(stage_sig)\n",
    "cell4 remove phys_bounds from stage call",
)

# 4g. Remove bounds_pct from rows_summary
replace_in_cell(cell4,
    "                    'p999_abs_uV': m['p999_abs_uV'],\n"
    "                    'flat_pct': m['flat_pct'],\n"
    "                    'bounds_pct': m['bounds_pct'],\n"
    "                    'hist_extreme_pct': m['hist_extreme_pct'],\n",
    "                    'p999_abs_uV': m['p999_abs_uV'],\n"
    "                    'flat_pct': m['flat_pct'],\n"
    "                    'hist_extreme_pct': m['hist_extreme_pct'],\n",
    "cell4 remove bounds_pct from rows_summary",
)

# 4h. Remove bounds_pct from rows_stage_summary
replace_in_cell(cell4,
    "                            'flat_pct': sm['flat_pct'],\n"
    "                            'bounds_pct': sm['bounds_pct'],\n"
    "                            'hist_extreme_pct': sm['hist_extreme_pct'],\n",
    "                            'flat_pct': sm['flat_pct'],\n"
    "                            'hist_extreme_pct': sm['hist_extreme_pct'],\n",
    "cell4 remove bounds_pct from rows_stage_summary",
)

# 4i. Remove 'At EDF bounds (%)' from metrics table in report
replace_in_cell(cell4,
    "                        ('At EDF bounds (%)', f\"{m['bounds_pct']:.3f}%\"),\n"
    "                    ",
    "                    ",
    "cell4 remove At EDF bounds from metrics table",
)

# 4j. progress.value at end of loop
replace_in_cell(cell4,
    "        progress.value = len(edf_files)\n",
    "        progress.value = len(cdt_files)\n",
    "cell4 progress.value end",
)

# 4k. Curry-only: cap the shared time-series / histogram amplitude limit at a wide
# physiological ceiling. DC-coupled Curry data (no export clipping) can carry large
# drift/artifacts that blow up the p99.9 autoscale and crush the real EEG; capping keeps
# ~all physiological signal visible while clean low-amplitude channels still zoom in below.
replace_in_cell(cell4,
    "            _p999 = [m['p999_abs_uV'] for m in ch_metrics.values()]\n"
    "            y_lim_ts = float(max(_p999)) if _p999 else 1.0\n",
    "            # DC-coupled data carries no export clipping, so drift/artifacts can push the\n"
    "            # p99.9 autoscale far past physiological range and crush the real EEG. Cap the\n"
    "            # shared time-series / histogram amplitude limit at a wide physiological ceiling\n"
    "            # (uV); clean low-amplitude channels still auto-zoom below the cap.\n"
    "            DISPLAY_YLIM_UV = 500.0\n"
    "            _p999 = [m['p999_abs_uV'] for m in ch_metrics.values()]\n"
    "            y_lim_ts = min(float(max(_p999)), DISPLAY_YLIM_UV) if _p999 else 1.0\n",
    "cell4 cap time-series y-limit at physiological ceiling",
)

print("\n=== Global replacements ===")

# Global: edf_path → cdt_path (only in code cells, cell 4 effectively)
replace_all_cells("edf_path", "cdt_path", "edf_path → cdt_path")

# Global: edf_files → cdt_files (cells 3 and 4 residual after targeted replacements)
replace_all_cells("edf_files", "cdt_files", "edf_files → cdt_files")

# Global: edf_relative → cdt_relative (variable naming consistency)
replace_all_cells("edf_relative", "cdt_relative", "edf_relative → cdt_relative")

# Global: edf_reports_dir → cdt_reports_dir
replace_all_cells("edf_reports_dir", "cdt_reports_dir", "edf_reports_dir → cdt_reports_dir")

# Update markdown title cell
for cell in nb["cells"]:
    if cell["cell_type"] == "markdown":
        src = get_src(cell)
        if "quality" in src.lower() or "overview" in src.lower():
            new_src = src.replace(
                "EDF Quality Overview",
                "Quality Overview — Curry 9 (.cdt)"
            ).replace(
                "quality_overview_voila",
                "5_quality_overview_curry_voila"
            )
            if new_src != src:
                set_src(cell, new_src)
                print("  ✓ markdown title updated")

# ---------------------------------------------------------------------------
# Validate and write
# ---------------------------------------------------------------------------
print("\n=== Validation ===")
try:
    out_str = json.dumps(nb, ensure_ascii=False, indent=1)
    json.loads(out_str)
    print(f"JSON valid. n_cells={len(nb['cells'])}")
except Exception as e:
    print(f"✗ JSON invalid after edits: {e}")
    sys.exit(1)

if errors:
    print(f"\n⚠ {len(errors)} pattern(s) not found: {errors}")

with open(DST, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write("\n")
print(f"Written: {DST}")
