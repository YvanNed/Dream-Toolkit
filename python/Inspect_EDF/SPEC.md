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
1. **Scan** — Select data folder, hypnogram suffix, and output suffix → auto-detects files recursively; reports unique label configurations, flags boundary `?` epochs, and highlights suspicious labels (unrecognized, numeric > 5, length > 3, or containing special characters). Also reports: how many remapped files already exist with the output suffix; a list of `.txt` files matching neither suffix (informational, only when both suffixes are defined); automatically exports `mid_question_epochs_to_verify.tsv` if mid-recording `?` epochs are found. Optional checkbox to exclude participants whose remapped file already exists from all downstream processing.
2. **Mid-recording `?` correction** *(shown only when mid-recording `?` epochs are found)* — One-by-one navigation through all `?` epochs; displays ±N epochs of context (adjustable slider, default 5); individual correction with auto-advance; warns user to cross-check against scoring software. Clicking "Corrections done" updates `mid_question_epochs_to_verify.tsv` with the assigned labels.
3. **Remap labels** — Per-configuration accordion widget with combobox suggestions pre-filled from `DEFAULT_MAPPING`; suspicious labels are highlighted in orange with inline epoch context; warns if the final mapping leaves non-AASM labels; confirmation required before proceeding
4. **Save** — Writes remapped hypnograms next to originals using the output suffix defined in Section 1; end message confirms completion and recalls the suffix used
5. **Verify** — Before/after configuration summary; verdict fails only if non-AASM labels remain (multiple configurations with valid AASM labels are acceptable — e.g. insomnia patients legitimately missing N3)

**Key constants:**
- `DEFAULT_MAPPING`: `0→W`, `1→N1`, `2→N2`, `3→N3`, `4→N3`, `5→R`, `?→W`, `S1→N1`…
- `STANDARD_LABELS`: `{W, N1, N2, N3, R}`
- `ACCEPTABLE_LABELS`: `STANDARD_LABELS | {MT}` — used for AASM compliance warnings (MT = movement time is tolerated)

**Outputs**:
- One `.txt` file per participant with the output suffix (e.g. `_Hypnogram_remapped.txt`), one label per line
- `mid_question_epochs_to_verify.tsv` — written to `<data_folder>/` at scan time whenever mid-recording `?` epochs are found; columns: `participant_id`, `epoch_index`, `epoch_time_sec`, `total_epochs`, `context` (±5 epochs), `corrected_label`; updated with assigned labels after Section 2 corrections

`check_hypno_config.py` is the legacy script that preceded this notebook; kept for reference.

### 4. Spectral Analysis (`SpectralPower_&_AperiodicFit_PSG.py`)

Full PSG spectral pipeline: epoch rejection → PSD (Welch, 4 s windows) → aperiodic fit (SpecParam) → frequency band power extraction (Delta, Theta, Alpha, Sigma, Beta) → group-level statistics. Reads the channel remapping JSON produced by tool #2.

## Planned modules (in development)

### Preprocessing pipeline

A preprocessing module is being developed to follow the inspection and channel-harmonization steps. It is structured in three phases.

**Phase 1 — Quality overview (per participant, per channel)**
Implemented as `tools/quality_overview_voila.ipynb`. Produces one `mne.Report` HTML per participant. For each EEG channel: signal amplitude histogram with Savitzky-Golay smooth + peak detection, time series, metrics table, and a YASA hypnospectrogram (0.1–40 Hz bandpass applied per-channel just before plotting). Flags suspect channels for priority inspection. At the end of each run, generates `dataset_overview.html` — a single-page dataset-level summary with statistics and distribution plots per electrode, consumed by Phase 2 to identify channels to exclude.

**EDF scan**: recursive (`rglob('*.edf')`), so datasets organized in subfolders (e.g. `group1/`, `group2/`) are fully covered without needing to run the tool per subfolder.

**Hypnogram suffix auto-detection**: when the data folder is selected, the tool scans `.txt` files next to each EDF, counts candidate suffixes (all files per EDF are matched — no early break), and auto-fills the `Hypno suffix:` widget with the most common suffix. In case of equal counts, the longest suffix is preferred (more specific = remapped/processed version). All detected suffixes and their counts are shown in a colour-coded info label below the widget (green = all EDFs matched, orange = partial match or no files found).

**Live output warnings**: if a hypnogram is not found, fails to load, or has a length mismatch with the EDF, a plain-text `⚠` warning is printed in the notebook output area immediately after the per-participant result line (in addition to the yellow banner already shown in the HTML report's Overview section).

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

**Phase 2 — Preprocessing + epoch rejection**
Reads the `quality_summary.tsv` produced by Phase 1 (filtering on the `exclude` column) to determine which channels to exclude. Applies filtering (default proposal: FIR `l_freq=0.1, h_freq=40, phase=zero-double, fir_window=hamming`), resampling (with anti-aliasing), and re-referencing (driven by `config_param/remap_reref_persubject.json`). Epochs are 30 s, aligned to the hypnogram. Rejection criteria: peak-to-peak amplitude, point-to-point (gradient) for high-frequency artifacts, flat-signal detection, manually-scored events (limb movement, apnea, micro-arousal — format TBD: EDF annotations vs. XML vs. CSV), and aperiodic 1/f fit quality. Output: one `mne.Report` per participant including a channel × epoch heatmap colored by the rejection method that flagged each epoch, plus a `rejection_log.tsv` and cleaned epochs saved as `.fif`.

**Phase 3 — Interactive review**
Click on the rejection heatmap to inspect any flagged epoch (time series, PSD, etc.). Target user has no programming experience — framework choice will favor the simplest deployable option (Voila / Plotly).

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
