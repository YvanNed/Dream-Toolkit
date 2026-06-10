# Dream-Toolkit — Inspect_EDF — Project Specification

## Project overview

Research toolkit for sleep scientists to **inspect, validate, and preprocess EEG/PSG databases in EDF format** without requiring programming expertise. The target users are sleep researchers who work with polysomnography (PSG) data and need to perform quality control and preprocessing across multi-subject, multi-session datasets.

The toolkit provides two delivery modes for each tool:
- **Voila notebooks** (`.ipynb` launched with `voila`): hide all code and expose a clean interactive GUI for non-programmers
- **Standard Jupyter notebooks** (`.ipynb`): show all code for debugging and customization
- **Python scripts** (`.py`): batch processing on full datasets

## Repository structure

```
Inspect_EDF/
├── environment.yml              # Conda environment definition (inspect_edf)
├── CLAUDE.md                    # Development rules and design decisions
├── SPEC.md                      # This file — project specification
├── tools/
│   ├── inspect_edf.ipynb                      # EDF parameter inspector (Jupyter)
│   ├── inspect_edf_voila.ipynb                # EDF parameter inspector (Voila GUI)
│   ├── inspect_edf_perdataset.py              # Batch EDF inspection per dataset
│   ├── inspect_edf_perparticipant.py          # Batch EDF inspection per participant
│   ├── select&remap_channels_edf.ipynb        # Channel selection & harmonization (Jupyter)
│   ├── select&remap_channels_edf_voila.ipynb  # Channel selection & harmonization (Voila)
│   ├── check_hypno_config.py                  # Hypnogram validation (legacy script)
│   ├── remap_hypno.ipynb                      # Hypnogram label remapping (Jupyter)
│   ├── remap_hypno_voila.ipynb                # Hypnogram label remapping (Voila GUI)
│   ├── quality_overview_voila.ipynb           # Phase 1 quality overview (Voila GUI)
│   ├── preprocessing_voila.ipynb             # Phase 2 preprocessing + epoch rejection (Voila GUI)
│   ├── SpectralPower_&_AperiodicFit_PSG.py    # Spectral analysis pipeline
│   ├── generate_test_data.py                  # Inject controlled defects into a clean EDF (test fixtures)
│   ├── test_data/                             # Real EDF fixtures + generated defective files + manifest
│   ├── images/                                # Reference images for quality checks
│   ├── preprocessing_phase1_example_scripts/  # Draft/example scripts used during Phase 1 development
│   └── old/                                   # Versioned development notebooks (archive)
```

**Sibling directory** `../Check_EDF/` contains exploratory notebooks used during development (not production tools).

## Conda environment

Defined in `environment.yml`. Key packages:
- **Python 3.12.10**
- **MNE 1.9** — EDF reading, epoching, signal processing
- **YASA 0.6** — sleep staging, hypnogram handling, spectral helpers
- **pandas 2.2, numpy 2.2** — data manipulation
- **voila 0.5, ipywidgets 8.1, ipyfilechooser 0.6** — interactive GUI layer
- **chardet 5.2** — encoding detection for EDF headers
- **edfio** — EDF read/write (used directly for export with per-channel physical range control)
- **specparam** — aperiodic/periodic spectral decomposition (1/f fitting)

## How to run the tools

### Interactive Voila apps (no-code mode)
```bash
conda activate inspect_edf
voila tools/inspect_edf_voila.ipynb
voila "tools/select&remap_channels_edf_voila.ipynb"
voila tools/remap_hypno_voila.ipynb
voila tools/quality_overview_voila.ipynb
voila tools/preprocessing_voila.ipynb
```

### Standard Jupyter notebooks
```bash
conda activate inspect_edf
jupyter notebook tools/inspect_edf.ipynb
```

### Batch Python scripts
```bash
conda activate inspect_edf
python tools/inspect_edf_perdataset.py
python tools/inspect_edf_perparticipant.py
```

## Tool descriptions

### 1. EDF Inspector (`inspect_edf_voila.ipynb`, `inspect_edf_perdataset.py`, `inspect_edf_perparticipant.py`)

Inspects EDF file parameters across an entire dataset **without loading signal data**. Reads EDF headers using a custom binary parser to handle encoding edge cases robustly (see design decisions in CLAUDE.md).

**Checks performed for EEG, EOG, and ECG channels:**
- Channel configuration and montage consistency across participants
- Sampling frequency consistency
- Filter settings consistency
- Signal units
- Inverted polarity (physical_min > physical_max)
- Signal clipping (dynamic range ≤ 500 µV)
- Poor resolution (dynamic range ≥ 0.1 µV per digital unit)

**Outputs** (all written to `<study_folder>/summary_inspection/`, created on first run with a `README.md`):
- `FULL_summary_table_edf.tsv` — all parameters for all channels/files
- `EEG_summary_table.tsv`, `EOG_summary_table.tsv`, `ECG_summary_table.tsv`
- `EEG_inverted_polarity_edf.tsv`, `EEG_bad_dynamic_range_edf.tsv`, `EEG_bad_resolution_edf.tsv`
- `EDF_inspection_report.html`, `EDF_perParticipant_report.html`
- `failed_edf_read.tsv` — files that could not be read
- `README.md` — describes each output file and which tool generates it

### 2. Channel Selection & Remapping (`select&remap_channels_edf_voila.ipynb`)

Interactive tool to select channels of interest and harmonize their labels across a heterogeneous dataset. Produces a JSON remapping configuration consumed by downstream analysis tools.

**Channel detection robustness**: The EDF scan section uses a 3-condition mask to identify EEG/EOG channels:
1. `transducer_type` contains `EEG`, `AGAGCL ELECTRODE`, or `EOG` (standard acquisition systems)
2. Channel name contains `EOG` (fallback for EOG channels)
3. Channel name matches `KNOWN_EEG_CHANNEL_RE` — an anchored regex covering the full 10-10 system, mastoids (M1/M2, A1/A2), and common EOG labels (LOC, ROC, E1/E2)

Condition 3 is required for EDFs exported by `mne.export.export_raw()`, which writes an empty `transducer_type` field.

**Section 5 — JSON save**: section 5 exposes a single "Preview & Save" button. Clicking it generates the per-participant JSON dict, writes it immediately to `<data_folder>/config_param/remap_reref_persubject.json`, and displays the preview and usage instructions. There is no separate save step. `mne_reref_plan.json` is no longer generated (it was redundant with `remap_reref_persubject.json` which already carries `ref_channels` per participant).

### 3. Hypnogram Label Remapping (`remap_hypno_voila.ipynb`, `remap_hypno.ipynb`)

Interactive tool to harmonize sleep stage labels across a heterogeneous database, converting different scoring conventions (e.g. `0,1,2,3,4` or `W,S1,S2,S3,S4`) to the standard AASM format (`W`, `N1`, `N2`, `N3`, `R`).

**Workflow (5 sections):**
1. **Scan** — Select data folder, hypnogram suffix, and output suffix → auto-detects files recursively; reports unique label configurations, flags boundary `?` epochs, and highlights suspicious labels (any label not in `DEFAULT_MAPPING`, including short alphanumeric labels like `U` or `M`). Also reports: how many remapped files already exist with the output suffix; a list of `.txt` files matching neither suffix (informational, only when both suffixes are defined); automatically exports `mid_uncertain_epochs_to_verify.tsv` if mid-recording `?` epochs or unexpected labels are found. If the TSV already exists (e.g. edited manually between sessions), it is loaded automatically and corrections are pre-applied in memory — the file is never overwritten at scan time. Optional checkbox to exclude participants whose remapped file already exists from all downstream processing.
2. **Uncertain and unexpected epoch review** — Two independent sub-widgets, each activated when the relevant issue is found:
   - **2a — Mid-recording `?` correction**: one-by-one navigation through all mid-recording `?` epochs; displays ±N context epochs (adjustable slider, default 5); per-epoch assignment with auto-advance; warns user to cross-check in scoring software. Clicking "Corrections done" writes the unified `mid_uncertain_epochs_to_verify.tsv` with all current corrections.
   - **2b — Unexpected label review**: navigation through all epochs carrying a label not in `DEFAULT_MAPPING` (e.g. `U`, `M`); "Jump to label" dropdown to filter by label type; context table identical to 2a; **Apply to this epoch** for per-epoch assignment (auto-advances); **Apply to all in this file** to batch-replace all occurrences of the current label within the current participant. Corrections update `hypno_data` in memory; fully corrected labels disappear from Section 3 configs on next run. Corrections are persisted in the same `mid_uncertain_epochs_to_verify.tsv` when "Corrections done" is clicked.
3. **Remap labels** — Per-configuration accordion widget with combobox suggestions pre-filled from `DEFAULT_MAPPING`; suspicious labels are highlighted in red with inline epoch context; warns if the final mapping leaves non-AASM labels; confirmation required before proceeding
4. **Save** — Writes remapped hypnograms next to originals using the output suffix defined in Section 1; end message confirms completion and recalls the suffix used
5. **Verify** — Before/after configuration summary; verdict fails only if non-AASM labels remain (multiple configurations with valid AASM labels are acceptable — e.g. insomnia patients legitimately missing N3)

**Hypnogram suffix auto-detection**: when the data folder is selected, the tool scans `.txt` files next to each EDF, counts candidate suffixes (all files per EDF are matched — no early break), and auto-fills the `Hypnogram suffix:` widget with the most common suffix. In case of equal counts, the **shortest** suffix is preferred — the goal is to select the raw (unremapped) hypnogram, not an already-processed version. This is the reverse of the strategy used in `quality_overview_voila` and `preprocessing_voila`, which prefer the longest suffix among candidates appearing for ≥50% of the maximum count. All detected suffixes and their counts are shown in a colour-coded info label below the widget (green = all EDFs matched, orange = partial match or no files found).

**Key constants:**
- `DEFAULT_MAPPING`: `0→W`, `1→N1`, `2→N2`, `3→N3`, `4→N3`, `5→R`, `?→W`, `S1→N1`…
- `STANDARD_LABELS`: `{W, N1, N2, N3, R}`
- `ACCEPTABLE_LABELS`: `STANDARD_LABELS | {MT}` — used for AASM compliance warnings (MT = movement time is tolerated)

**Outputs**:
- One `.txt` file per participant with the output suffix (e.g. `_Hypnogram_remapped.txt`), one label per line
- `mid_uncertain_epochs_to_verify.tsv` — written to `<data_folder>/` at scan time when mid-recording `?` epochs or unexpected labels are found; columns: `participant_id`, `epoch_index`, `epoch_time_sec`, `total_epochs`, `original_label` (`?` for mid-recording unscored epochs, or the raw unexpected label e.g. `U`, `M`), `context` (±5 epochs), `corrected_label`; updated with all corrections (both types) after "Corrections done" is clicked. If the file already exists at scan time it is loaded and applied in memory instead of being overwritten — corrections whose `original_label` no longer matches the current hypnogram value (re-scored since the TSV was written) are flagged as conflicts and shown in a warning panel in Section 2.

`check_hypno_config.py` is the legacy script that preceded this notebook; kept for reference.

### 4. Spectral Analysis (`SpectralPower_&_AperiodicFit_PSG.py`)

Full PSG spectral pipeline: epoch rejection → PSD (Welch, 4 s windows) → aperiodic fit (SpecParam) → frequency band power extraction (Delta, Theta, Alpha, Sigma, Beta) → group-level statistics. Reads the channel remapping JSON produced by tool #2.

## Planned modules (in development)

### Preprocessing pipeline

A preprocessing module is being developed to follow the inspection and channel-harmonization steps. It is structured in three phases.

**Phase 1 — Quality overview (per participant, per channel)**
Implemented as `tools/quality_overview_voila.ipynb`. Produces one `mne.Report` HTML per participant. For each EEG channel: signal amplitude histogram with Savitzky-Golay smooth + peak detection, time series, metrics table, and a YASA hypnospectrogram (0.1–40 Hz bandpass applied per-channel just before plotting). Flags suspect channels for priority inspection. At the end of each run, generates `dataset_overview.html` — a single-page dataset-level summary with statistics and distribution plots per electrode, consumed by Phase 2 to identify channels to exclude.

**EDF scan**: recursive (`rglob('*.edf')`), so datasets organized in subfolders (e.g. `group1/`, `group2/`) are fully covered without needing to run the tool per subfolder.

**Hypnogram suffix auto-detection**: when the data folder is selected, the tool scans `.txt` files next to each EDF, counts candidate suffixes (all files per EDF are matched — no early break), and auto-fills the `Hypno suffix:` widget. Selection rule: among suffixes whose count is ≥ 50% of the maximum count, the **longest** is preferred (more specific = remapped/processed version); count is the tiebreaker when lengths are equal. The 50% threshold prevents rare or accidental long suffixes from winning when remapping is far from complete. All detected suffixes and their counts are shown in a colour-coded info label below the widget (green = all EDFs matched, orange = partial match or no files found). The same rule applies in `preprocessing_voila`.

**Live output warnings**: if a hypnogram is not found, fails to load, has a length mismatch with the EDF, or contains unrecognised stage labels, a plain-text `⚠` warning is printed in the notebook output area immediately after the per-participant result line (in addition to the yellow banner already shown in the HTML report's Overview section).

**Hypnogram label validation**: after mapping labels to YASA integers (`W→0`, `N1→1`, `N2→2`, `N3→3`, `R→4`), the tool checks for unrecognised labels. `MT` (movement time) is silently tolerated — YASA treats it as an artifact epoch (NaN). Any other unrecognised label triggers a warning. Two severity levels: if > 10% of epochs are unrecognised, `hypno_vec` is set to `None` (spectrogram skipped, error-level warning pointing to `remap_hypno_voila`); if ≤ 10%, `hypno_vec` is kept but a warning lists the unrecognised labels and their count. This catches hypnograms that were not yet remapped to AASM convention (e.g. `S1/S2/S3/S4` or raw numeric labels).

**Spectrogram fault tolerance**: the `yasa.hypno_upsample_to_data()` + `yasa.plot_spectrogram()` calls are wrapped in a `try/except`. If YASA raises an exception, the channel's spectrogram section in the HTML report shows the error message and processing continues to the next channel and participant without interruption.

**Report structure**: the MNE Report has one section per channel (e.g. "C3", "Fp1") plus an "Overview" section with the flag summary. Each channel section contains: histogram → time series → metrics table → spectrogram. The selector widget shows how many participants already have an existing report before the "Skip participants with an existing report" checkbox.

**Metrics table — interpretation column**: for `Flat signal (%)`, `At EDF bounds (%)`, and `Extreme histogram (%)`, the interpretation cell shows the actual threshold value used for that run (read from the widget at run time, e.g. "flag if > 3.5%"), making each report self-documenting. The three threshold-based entries are generated dynamically inside `run_analysis` via a local `interpretations` dict that overrides the static `METRIC_INTERPRETATIONS` dict.

**Histogram axis scaling**: the shared Y-axis for histograms is computed from non-suspect channels only (those not flagged for `flat_pct` or `std_uV`). A flat or dead channel concentrates all samples in 1–2 bins and would otherwise crush healthy distributions. Fallback to global max if all channels are suspect. The shared X-axis is `±x_lim_hist` where `x_lim_hist = p99.9 of |amplitude| across all channels` (no cap) — this allows cross-channel comparison: a clipped channel appears compressed, a dead channel appears as a narrow central spike.

**Flagging criteria and default thresholds:**

| Metric | Threshold | Detects |
|--------|-----------|---------|
| `flat_pct` | > 3.5% | Flat segments, dead channels |
| `bounds_pct` | > 1.0% | Saturation at EDF physical-range limits |
| `n_peaks` | ≥ 2 | Bimodal (DC drift) or multimodal (quantization) distribution |
| `std_uV` | < 5.0 µV | Dead / near-dead channels |
| `hist_extreme_pct` | > 1.0% | In-range clipping (saturation within declared EDF range) |

`flat_pct` = fraction of consecutive sample pairs with \|diff\| < `max(2×ADC_step, 0.06 µV)`. `bounds_pct` = fraction of samples within 0.5 µV of the EDF physical_min/max header limits. `hist_extreme_pct` = fraction of samples in the outermost histogram bins. **Kurtosis is intentionally NOT used as a flagging criterion** — normal PSG EEG has physiologically high kurtosis (spindles, K-complexes produce values of 100–500), making it unreliable without per-subject normalization. `p99_abs_uV` and `p999_abs_uV` (99th and 99.9th percentile of |amplitude|) are recorded as informational metrics in `quality_summary.tsv` but are not currently used for flagging; they are useful for cross-channel and cross-dataset amplitude comparison.

**`dataset_overview.html` — dataset-level summary**: generated at the end of every run from the full cumulative `quality_summary.tsv` (reflects all participants processed to date, not just the current run). Contains two levels:
- **Global section (all electrodes pooled)**: stats table (mean / median / p5 / p25 / p75 / p95 per metric), n_peaks frequency table by electrode (flags DC drift and quantization cases), grouped boxplots (one subplot per key metric, x-axis = electrode — compares electrodes side by side).
- **Per-electrode sections**: stats table for that electrode only, boxplots of each metric's distribution across participants with individual data points overlaid (reveals outlier participants for that channel).

Key metrics shown in plots: `std_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `p99_abs_uV`, `p999_abs_uV`. All numeric metrics (`mean_uV`, `kurtosis`, `skewness` included) appear in the stats tables. HTML is standalone (figures embedded as base64 PNG, no external dependencies).

**Outputs per run**:
- `<data_folder>/reports_quality_overview/<relative_subfolder>/<file_id>_quality_overview.html` — HTML reports mirror the EDF subfolder structure under `reports_quality_overview/`
- `<data_folder>/reports_quality_overview/quality_summary.tsv` — all numeric metrics for all channels, cumulative across runs (always at root); columns: `file_id`, `channel`, `mean_uV`, `std_uV`, `kurtosis`, `skewness`, `p99_abs_uV`, `p999_abs_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `n_peaks`, `suspect_reason`, `exclude` (last two columns)
- `<data_folder>/reports_quality_overview/dataset_overview.html` — dataset-level statistics and plots (always at root, regenerated each run)
- `<data_folder>/reports_quality_overview/failed_files.tsv` — files that could not be read (at root)

**End-of-run summary** (printed in the notebook output): participants processed, participants with ≥1 flagged channel, total flagged channels, files failed to load, path to `dataset_overview.html`.

**Phase 2 — Preprocessing + epoch rejection** (`tools/preprocessing_voila.ipynb`)

Implemented as a Voila notebook with four sections: (1) path configuration, (2) preprocessing and rejection parameters, (3) participant selection, (4) processing loop.

**Inputs**:
- `quality_summary.tsv` from Phase 1 — `exclude` column identifies channels to drop before preprocessing
- `remap_reref_persubject.json` from `select&remap_channels_edf` — drives channel remapping and re-referencing per participant
- Raw EDF files and remapped hypnograms (default suffix `_Hypnogram_remapped.txt`)

**Channel-name handling when loading EDF data from notebook outputs** (critical):

The EDF files on disk carry their **original** acquisition channel names (e.g. `Fp1, C3, O1, A2`), usually alongside non-montage channels (ECG, respiration, SpO2, position) that may be sampled at a **higher rate** than the EEG (e.g. ECG at 512 Hz while EEG is at 256 Hz). The config/report files refer to channels by their **harmonized (remapped)** names:
- `remap_reref_persubject.json` stores the mapping original → remapped, e.g. `{"Fp1":"F","C3":"C","O1":"O","A2":"M"}`. Its **keys are the original EDF names**, its values the harmonized names.
- `quality_summary.tsv` lists each channel under its **remapped** name, because `quality_overview_voila` renames the channels (`raw.rename_channels`) *before* computing per-channel metrics. So the `exclude` flags and the channel checkboxes derived from it are expressed in remapped names.

Both tools load the EDF with the same robust pattern — `include=` evaluated **at read time** using the **original** names taken straight from the remap keys, with `preload=False` so no signal is read yet:
```python
raw = mne.io.read_raw_edf(str(edf_path), preload=False, encoding='latin-1',
                          include=list(sub_config['remap'].keys()), verbose=False)
raw, _ = drop_suffix_duplicates(raw)
raw.rename_channels(adapt_remap_dict_to_suffixes(raw, sub_config['remap']))
```
Why `include=` at read time (and **not** a lazy `raw.pick(...)` after loading):
- **Avoids an MNE EDF-reader bug**: reading a channel *subset* lazily (`read_raw_edf(preload=False)` with no `include`, then `pick`, then `load_data`) raises a bare `AssertionError` (`max(n_smp_read) == smp_exp` in `mne/io/edf/edf.py`) whenever the subset **excludes the file's highest-sampling-rate channel** — exactly the case for a PSG where the ECG (512 Hz) is dropped and only EEG (256 Hz) is kept. Passing `include=` at read time rebuilds the record structure for the included set, so the assertion never fires. (This is why a bare, message-less error surfaced during development.)
- **Preserves the native sampling rate**: with only EEG channels included, MNE's common sampling rate stays at the EEG rate (256 Hz). Loading *all* channels first (`preload=True` on the full montage) would upsample the EEG to the file max (512 Hz) — a silent change of semantics.
- **No inverse-remap dictionary**: the include list is `list(remap.keys())` directly, identical to `quality_overview_voila` — the two tools stay consistent.

The two tools then differ only in what follows:
- `quality_overview_voila` analyses **every** remap channel, so it keeps them all (uses `preload=True`, no further selection).
- `preprocessing_voila` lets the user deselect channels in the UI (selection in **remapped** names). After the rename it drops the de-selected channels, then calls `raw.load_data()` — so the signal is read from disk **only for the channels actually kept**:
  ```python
  present = [ch for ch in selected_channels if ch in raw.ch_names]   # remapped namespace, post-rename
  raw.drop_channels([ch for ch in raw.ch_names if ch not in present])
  raw.load_data()
  ```
  A `present`-empty guard (e.g. if the rename failed) marks the participant as failed and skips it, so a single bad file never crashes the run.

**Shared utility functions** (defined identically in `quality_overview_voila` and `preprocessing_voila`, used right after `read_raw_edf`):
- `drop_suffix_duplicates(raw)` — MNE ≥ 1.8 appends `-0`/`-1` to duplicate channel names; this keeps only the `-0` variant and drops the rest. Returns `(raw, dropped_list)`.
- `adapt_remap_dict_to_suffixes(raw, remap_dict)` — rewrites the remap dict so a base name (`Fp1`) matches MNE's suffixed name (`Fp1-0`) when duplicates were present, so `rename_channels` still applies.

**Preprocessing steps** (applied in this order, each optional via widget):
1. **Resampling** — `raw.resample(target_freq, npad='auto')`. Target frequency chosen by user; applied before filtering to avoid aliasing. Step is skipped if checkbox is unchecked.
2. **Re-referencing** — applied as specified in JSON config per participant: `'average'` → common average reference; `[list]` → subtract listed channel(s) then drop them; empty → no re-referencing.
3. **Bandpass filter** — FIR zero-double-pass Hamming window, defaults `l_freq=0.1 Hz, h_freq=40 Hz`. Applied via `raw.filter(..., method='fir', phase='zero-double', fir_window='hamming', fir_design='firwin')`.

**Epoching**: 30-second fixed-length epochs created with `mne.make_fixed_length_epochs(raw, duration=30)`. Sleep stage assigned to each epoch from the hypnogram; epochs at the tail beyond the hypnogram length are discarded.

**Epoch rejection — five methods**:

All methods operate on the raw epoch data in µV (`epochs.get_data() * 1e6`, shape `n_epochs × n_channels × n_times`). Rejection masks are boolean arrays of shape `(n_epochs, n_channels)` — a `True` entry means that (epoch, channel) pair was flagged. An epoch is considered **rejected** if *any* channel is flagged by *any* method.

| Method | Signal feature | Default threshold | Notes |
|--------|---------------|-------------------|-------|
| **Amplitude** | Peak-to-peak = `max(epoch) − min(epoch)` | W: 300, N1: 250, N2/N3: 200, REM: 250 µV | Per-stage threshold; W/REM more lenient because muscle and eye-movement artefacts are physiologically common in those stages. Equivalent to MNE's `drop_bad(reject=...)` criterion. |
| **Flat signal** | Peak-to-peak < threshold | 1 µV | Detects disconnected electrodes or amplifier saturation within a single epoch. Logically identical to MNE's `drop_bad(flat=...)` criterion: both compare `ptp` against a low-amplitude threshold. |
| **Gradient** | `max(|diff(epoch)|)` across time | 100 µV/sample | Maximum sample-to-sample absolute difference; sensitive to sudden jumps, electrode pops, and movement artefacts not captured by peak-to-peak. `diff` and `max` both operate on `axis=-1` (time axis) to handle the 3D `(n_epochs, n_channels, n_times)` array correctly. |
| **Manual events** *(placeholder)* | Manually scored events (limb movement, apnea, micro-arousal) | — | Section reserved in code with `# TODO` comment. Format TBD: EDF annotations vs. XML vs. CSV. Not yet implemented. |
| **1/f fit quality** | Specparam aperiodic fit on Welch PSD (4 s windows, 2–30 Hz, `aperiodic_mode='fixed'`, `max_n_peaks=0`) | MAE > 0.15 OR R² < 0.95 | Fit restricted to ≥2 Hz to limit influence of slow-wave non-stationarity. `max_n_peaks=0` skips peak detection for speed (we only need the aperiodic metrics). A failed fit is treated as a double flag (both error and R²). |

**Heatmap — `_rejection_heatmap.png`**: channels (Y-axis) × epochs (X-axis); each cell coloured by the flagging method with priority encoding when multiple methods fire. A hypnogram strip is drawn above the main heatmap. Colour scheme: dark purple = none, red = amplitude, blue = flat, orange = gradient, yellow = 1/f error, green = 1/f R², dark red = multiple. Title includes overall rejection percentage.

**Two-step QC approach** — Phase 2 does **not** drop epochs. It saves ALL epochs (including flagged ones) with an MNE `metadata` DataFrame attached, so downstream Phase 2b can inspect rejected epochs before finalising the rejection.

**Outputs per participant** (in `<derivatives_root>/sub-{file_id}/`):
- `{file_id}_all-epo.fif` — all epochs with `epochs.metadata` DataFrame (columns: `epoch_idx`, `stage`, `reject_flag`, `reject_method`, `flag_amplitude`, `flag_flat`, `flag_gradient`, `flag_1f_error`, `flag_1f_r2`)
- `{file_id}_rejection_mask.tsv` — per-(epoch, channel) rejection table; human-readable and manually editable before Phase 2b (columns: `epoch_idx`, `stage`, `channel`, `reject_flag`, plus one bool column per method)
- `{file_id}_rejection_log.tsv` — rejection counts per stage per method (columns: `file_id`, `stage`, `method`, `n_total`, `n_rejected`, `pct_rejected`)
- `{file_id}_rejection_heatmap.png` — channels × epochs colour-coded heatmap
- `{file_id}_preprocessing_report.html` — MNE HTML report with heatmap and rejection summary table

**Global output** (at `<derivatives_root>/`):
- `preprocessing_phase2_global_rejection.tsv` — concatenation of all `_rejection_log.tsv` files across participants
- `preprocessing_phase2_failed.tsv` — participants that could not be processed (EDF not found, config missing, hypno mismatch, etc.)

**Phase 2b — Interactive QC** *(future)*
Load `_all-epo.fif` + `_rejection_mask.tsv`, display the heatmap for navigation, show the raw signal of flagged epochs for visual inspection, allow manual override of individual entries in the mask, then save a final `{file_id}_clean-epo.fif` with only the validated-clean epochs.

### Test data infrastructure

`tools/generate_test_data.py` produces controlled-defect EDF files from a clean baseline (`tools/test_data/73.edf`). It uses the channel list from `tools/test_data/config_param/remap_reref_persubject.json` to keep only the relevant EEG channels (currently `A2, C3, Fp1, O1`), making each output ~110 MB. Each generated file injects **one** defect on a specified channel; a `combined` file mixes three defects on different channels for integration testing.

Each defect file is given a unique participant-like ID (731–738) so that `select&remap_channels_edf` does not group them all under participant `73`. The ID mapping is defined in `DEFECT_IDS` in `generate_test_data.py`. The script also writes entries for each generated file into `config_param/remap_reref_persubject.json` (copied from the `73` entry).

| Output file (`tools/test_data/`) | Channel | Injected defect | Expected detection |
|---|---|---|---|
| `731_clipping.edf` | Fp1 | Hard clip at ±75 µV | `hist_extreme_pct > 1.0%` |
| `732_combined.edf` | Fp1+C3+O1 | flat (Fp1) + drift (C3) + line noise (O1) | multiple metrics |
| `733_dc_drift.edf` | C3 | Sigmoidal +50 µV baseline shift, width 30 min | `bimodal_distribution` |
| `734_dead_channel.edf` | C3 | Whole channel scaled by 0.01 | `flat_pct > 3.5%` + `std_uV < 5 µV` |
| `735_flat_segment.edf` | Fp1 | 30 min of near-zero signal mid-recording | `flat_pct > 3.5%` |
| `736_line_noise.edf` | O1 | 50 Hz sinusoid, 50 µV peak-to-peak | `spectral_line_50hz` |
| `737_movement_bursts.edf` | Fp1 | 10 bursts (±300 µV, 2–5 s each) | `peak_to_peak_rejection_step2` |
| `738_quantization.edf` | C3 | Coarse quantization to 16 levels over ±200 µV | `multimodal_distribution` |

Each output also has copied companion files (`*_Hypnogram_Export.txt`, `*.edf.XML`, `*_event_xml.csv`). Ground truth is recorded in `tools/test_data/test_data_manifest.tsv`. File `100.edf` is kept as-is as a real-world fixture (naturally multi-peak histogram on Fp2). File `8_N1.edf` serves as a second clean-baseline sanity check.

Run with `python tools/generate_test_data.py` (idempotent — pass `--force` to regenerate). EDF files are written using `edfio.EdfSignal()` and `edfio.Edf()` directly.
