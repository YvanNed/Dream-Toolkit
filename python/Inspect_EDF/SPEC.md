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
│   ├── 1_inspect_edf.ipynb                      # EDF parameter inspector (Jupyter)
│   ├── 1_inspect_edf_voila.ipynb                # EDF parameter inspector (Voila GUI)
│   ├── 1_inspect_edf_perdataset.py              # Batch EDF inspection per dataset
│   ├── 1_inspect_edf_perparticipant.py          # Batch EDF inspection per participant
│   ├── 1bis_anonymize_edf_voila.ipynb              # EDF header anonymizer (Voila GUI)
│   ├── 1bis_anonymize_edf.ipynb                    # EDF header anonymizer (Jupyter)
│   ├── 2_select&remap_channels_edf.ipynb        # Channel selection & harmonization (Jupyter)
│   ├── 2_select&remap_channels_edf_voila.ipynb  # Channel selection & harmonization (Voila)
│   ├── check_hypno_config.py                  # Hypnogram validation (legacy script)
│   ├── 3_remap_hypno.ipynb                      # Hypnogram label remapping (Jupyter)
│   ├── 3_remap_hypno_voila.ipynb                # Hypnogram label remapping (Voila GUI)
│   ├── 4_remap_events_edf.ipynb                 # Event label harmonization (Jupyter)
│   ├── 4_remap_events_edf_voila.ipynb           # Event label harmonization (Voila GUI)
│   ├── 5_quality_overview_voila.ipynb           # Quality overview (Voila GUI)
│   ├── 6_preprocessing_voila.ipynb             # Preprocessing + epoch rejection (Voila GUI)
│   ├── 7_inspect_rejected_epochs_voila.ipynb    # QC of rejected epochs — Phase 2b (Voila GUI)
│   ├── 7_inspect_rejected_epochs_batch.py       # QC of rejected epochs — database-level batch report
│   ├── qc_rejected_epochs_lib.py                # Shared analysis/plotting for tool 7 (notebook + batch)
│   ├── 8_live_explore_1file.ipynb               # Interactive single-file explorer (Jupyter)
│   ├── 8_live_explore_1file_voila.ipynb         # Interactive single-file explorer (Voila GUI)
│   ├── 9_SpectralPower_&_AperiodicFit_PSG.py    # Spectral analysis pipeline
│   ├── generate_test_data.py                  # Inject controlled defects into a clean EDF (test fixtures)
│   ├── test_data/                             # Real EDF fixtures + generated defective files + manifest
│   ├── images/                                # Reference images for quality checks
│   ├── preprocessing_phase1_example_scripts/  # Draft/example scripts used during Phase 1 development
│   └── old/                                   # Versioned development notebooks (archive)
```

**Sibling directory** `../Check_EDF/` contains exploratory notebooks used during development (not production tools).

> **Editing the larger notebooks**: see the *rename-to-`.txt`* rule in CLAUDE.md — `Read`/`Edit`/`NotebookEdit` are blocked or size-capped on big `.ipynb` files.

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
voila tools/1_inspect_edf_voila.ipynb
voila tools/1bis_anonymize_edf_voila.ipynb
voila "tools/2_select&remap_channels_edf_voila.ipynb"
voila tools/3_remap_hypno_voila.ipynb
voila tools/4_remap_events_edf_voila.ipynb
voila tools/5_quality_overview_voila.ipynb
voila tools/6_preprocessing_voila.ipynb
voila tools/7_inspect_rejected_epochs_voila.ipynb
voila tools/8_live_explore_1file_voila.ipynb
```

### Standard Jupyter notebooks
```bash
conda activate inspect_edf
jupyter notebook tools/1_inspect_edf.ipynb
```

### Batch Python scripts
```bash
conda activate inspect_edf
python tools/1_inspect_edf_perdataset.py
python tools/1_inspect_edf_perparticipant.py
```

## Cross-cutting procedures (shared across tools)

Conventions and helpers reused by several tools. The tool sections below reference these by name
instead of restating them; only tool-specific deltas are kept inline.

- **Dual delivery (Jupyter + Voila + batch `.py`)**: every user-facing tool exists as a code-visible
  Jupyter notebook and a code-hidden Voila app (kept in sync); several also ship a batch `.py` twin.
  All forms must be maintained together when a tool changes.
- **Custom EDF header parser + sampling-frequency derivation**: EDF headers are read with a hand-written
  binary parser (robust to encoding/header edge cases; **never** used for signal data).
  `sampling_frequency = samples_per_record / duration_data_record` (the per-channel 8-byte field is the
  *number of samples per data record*, not the rate; the two coincide only when `duration_data_record == 1 s`).
  Kept as a string to preserve `sorted(set(...))` grouping. Shared by `1_inspect_edf*`,
  `2_select&remap_channels_edf*`, and `8_live_explore_1file*`. (Tool 1 carries the worked EDF+ examples.)
- **MNE EDF signal loading pattern**: when actual signal is needed, load with
  `mne.io.read_raw_edf(..., preload=False, include=list(remap.keys()))` — `include=` evaluated **at read
  time** (not a lazy `pick` afterwards) to avoid MNE's partial-read `AssertionError` when the highest-rate
  channel is excluded and to preserve the native EEG rate — then `drop_suffix_duplicates(raw)` and
  `raw.rename_channels(adapt_remap_dict_to_suffixes(raw, remap))`. The two helpers handle MNE ≥ 1.8's
  `-0`/`-1` suffixes on duplicate channel names. Shared by tools 5, 6, 8 (tool 6 documents the full
  rationale and its channel-deselection delta).
- **Path comparison normalization (`os.path.normcase`)**: whenever two filesystem paths, stems, or
  filenames are compared as strings (equality, `in`, `.isin()`, set/dict membership) **and the two sides
  can come from different sources** (one from disk, one from a stored TSV/JSON/widget), both sides are
  wrapped in `os.path.normcase(...)` **at the comparison only** — the stored/displayed value keeps its
  original case. On Windows this removes drive-letter case and `/`↔`\` mismatches (a no-op on POSIX).
  Applied across the skip/merge filters, EDF and hypnogram lookups, participant-id matching, and
  hypnogram-suffix matching of every tool that mixes disk paths with stored configuration.
- **Skip + cumulative-merge workflow**: processing tools decouple folder selection from running, show an
  "N / M already done" info line, offer a **"Skip already processed"** checkbox (on by default), and write
  outputs with **merge/replace** semantics (rows for the items processed this run replace their previous
  rows, all others kept) so the tables stay the full cumulative dataset across runs. Aggregated summary
  files are regenerated from every per-item file present, not just the current run. Shared by tools 1, 2,
  5, 6 (and the hypno/event variants), each keyed on its own identifier (`path`, participant id, `file_id`).
- **Lenient JSON loader**: every read of a config JSON parses strictly first and, on failure, repairs a
  single trailing comma before a closing `}`/`]` (a common hand-edit mistake) before retrying. Shared by
  the config readers of tools 2 and 4.
- **Hypnogram-suffix auto-detection**: on folder selection, the `.txt` files next to each EDF are scanned,
  candidate suffixes counted (all files per EDF matched), and the suffix widget auto-filled, with a
  colour-coded info label (green = all EDFs matched, orange = partial/none). **Selection rule differs by
  intent**: tool 3 prefers the **shortest** suffix on ties (target = raw, unremapped hypnogram); tools 5
  and 6 prefer the **longest** suffix among candidates appearing for ≥ 50% of the maximum count
  (target = the more specific remapped/processed version).
- **Event sourcing (CSV-first / XML-fallback)**: scored events are read via a shared `load_events()` —
  the Compumedics event CSV (`Name, Start, Duration`, default suffix `_event_xml.csv`) first, then the
  `<ScoredEvents>` of the `*.edf.XML` (`CMPStudyConfig`). Shared by tool 4 (harmonization) and tool 8
  (overlay). Tool 4 exposes the CSV suffix as an editable, auto-detected field (its `load_events()`
  takes a `csv_suffix=` argument). Tool 4 returns
  `(list-of-(name, start, duration), source)`; tool 8's `load_events()` returns the same events as a
  `Name/Start/Duration` **DataFrame** plus the `source` tag, because its overlay/navigator code consumes
  a DataFrame.
- **Proactive error handling**: per-item `try/except` with a **fatal** (add to a `failed` list and
  `continue`) vs **non-fatal** (`⚠` warning, continue) distinction; in Voila, every button callback and
  per-item loop is wrapped so a single failure never crashes the run or freezes the UI, and errors are
  always surfaced via a widget or `print()`.
- **Custom (non-AASM) sleep stages**: a project may intentionally keep stage labels outside the AASM set
  (`W/N1/N2/N3/R`), e.g. `N4` or a movement stage. They are declared **once** in a shared flat JSON
  `<data_folder>/config_param/custom_stages.json` (`{"custom_stages": ["N4", …]}`, order = display order),
  **written only by `3_remap_hypno`** and **read** by tools 5/6/7/8. Tools 5/6/7/8 each expose an editable
  `Custom stages` field **auto-filled from the JSON on folder/file selection** — a **volatile per-run
  override** that never rewrites the JSON (management stays in tool 3). Three helpers are **duplicated**
  across the tools (like `get_phys_bounds_uV`): `load_custom_stages(folder)`, `parse_custom_field(text)`,
  and `custom_stage_style(custom_stages) -> (stage_y, stage_colors, ytick_pos, ytick_labels)`. Custom
  stages stack **below N3** on every hypnogram axis (`N3=0 → -1, -2, …` in declaration order); the
  step-line stays gray with **REM in red** (YASA convention) and each custom stage in a fixed non-red
  palette (`#8dd3c7, #ffffb3, #bebada, #80b1d3, #fdb462, #b3de69, #fccde5, #d9d9d9`). Because
  `yasa.plot_spectrogram` / `yasa.Hypnogram` **hard-reject** any non-AASM label, tools 5 & 8 plot the
  hypnospectrogram with a custom **`plot_hypnospectrogram()`** that keeps YASA's stage-agnostic
  spectrogram core (`from yasa.plotting import spectrogram_lspopt`) and draws the hypnogram band itself —
  no new dependency (`spectrogram_lspopt` ships with the already-required `yasa`). Reading the JSON is
  non-fatal (`[]` on absent/corrupt); an unregistered non-AASM label keeps the old behaviour (warning +
  per-stage exclusion), never a crash.
  - **Flat/dead-epoch colour scaling**: YASA's percentile-based colormap range (`np.percentile(Sxx_dB,
    [trimperc, 100-trimperc])` over all pixels) washes the spectrogram out (uniform red, dead-epoch
    stripes) once the fraction of **fully-flat 30 s epochs** — spectrogram columns, **not** the
    channel's sample-level `flat_pct` — exceeds `trimperc` (2.5 %): a flat/disconnected epoch has
    ~zero power → ≈ −400 dB after the display filter, dragging `vmin` to that floor. The shared
    `plot_hypnospectrogram()` therefore excludes near-zero columns from the `vmin/vmax` percentiles
    (a column is valid when its peak dB is within 60 dB of the median epoch peak) and renders the
    excluded columns grey (`#d9d9d9`, "no signal"); clean channels (no dead epoch) are unaffected
    (byte-identical scale and image). The tool-7 *navigator* spectrogram is a separate plot (floored
    at −120 dB, p5–p99) and is left as-is. Diagnostic scripts:
    `tools/simple_hypnospectro_yasa_vs_fix.py` and `tools/compare_flat_spectrogram_fix.{py,ipynb}`.
- **Physical bounds in µV (`get_phys_bounds_uV`)**: MNE 1.9 stores an EDF channel's physical range as
  `physical_max + offset` in `raw._raw_extras` (not explicit `physical_min`/`physical_max`).
  `get_phys_bounds_uV()` reconstructs the µV bounds and **must scale both `physical_max` and `offsets`
  by `extras['units'][ch_idx] * 1e6`** — MNE keeps them in the channel's *native* EDF unit (`1e-6` µV,
  `1e-3` mV, `1.0` V), while `raw.get_data() * 1e6` is always µV. Without the scaling, any channel
  declared in mV (typical for Compumedics EOG/EMG/ECG, `physical_max = 1.0 mV`) is compared against a
  1.0 µV bound — 1000× too small — so `bounds_pct` flags ~98–100 % of a perfectly healthy signal; EEG
  (declared in µV) stays unaffected, which kept the bug latent until non-EEG channels were added.
  Defined in `5_quality_overview_voila`, **duplicated** in `8_live_explore_1file*` — keep in sync.
  (Verified on ICEBERG 117: EOG/EMG/ECG `bounds_pct` 97–98 % → <0.4 %, EEG unchanged.)
- **EOG/EMG/ECG channel-type detection (`detect_channel_types`)**: non-EEG channels are classified by
  **transducer type OR channel name** — EOG = transducer `EOG` / name `eog`; ECG = transducer
  `ECG`/`EKG` / name `ecg`/`ekg`; EMG = transducer `EMG` / name `emg`/`chin`/`menton` (the `chin|menton`
  aliases cover Compumedics chin-EMG labels). EEG uses `KNOWN_EEG_CHANNEL_RE` (full 10-10 + mastoids +
  literal `EEG`) or transducer `EEG`/`AGAGCL ELECTRODE`, excluding anything already matched as
  EOG/ECG/EMG. Defined in `8_live_explore_1file`; reuse it when a tool needs the scoring montage,
  pre-filled as an editable selection so the user can correct misses.
- **ipywidgets `Box` stray per-row scrollbars**: jupyter-widgets ships
  `.widget-box { box-sizing: border-box; overflow: auto; }`, so any bordered + padded `HBox`/`VBox` row
  whose children overflow the content box by even 1–2 px renders a per-row ▲▼ vertical scrollbar the user
  can scroll by accident (shifting the row content). Fix: set `overflow='hidden'` explicitly on such row
  layouts (applied to the Section 3 "Harmonize labels" rows of `4_remap_events_edf*`); keep the intended
  scroll only on the outer list container.

## Tool descriptions

### 1. EDF Inspector (`1_inspect_edf_voila.ipynb`, `1_inspect_edf_perdataset.py`, `1_inspect_edf_perparticipant.py`)

Inspects EDF file parameters across an entire dataset **without loading signal data**. Reads EDF headers using a custom binary parser to handle encoding edge cases robustly (see design decisions in CLAUDE.md).

**Sampling frequency derivation**: the per-channel 8-byte header field is the *number of samples per data record*, **not** the sampling frequency. The parser stores it as `samples_per_record` and computes `sampling_frequency = samples_per_record / duration_data_record`. In classic EDF the data-record duration is 1 s, so the two values coincide; but EDF+ files frequently use a different record duration (e.g. `0.1 s` → `40 / 0.1 = 400 Hz`, or `2 s` → `512 / 2 = 256 Hz`), so dividing by `duration_data_record` is required to match the rate reported by MNE. The computed value is kept as a string (e.g. `'256'`, `'400'`) to preserve the existing `sorted(set(...))` grouping and avoid TSV round-trip type changes. This applies to every tool sharing the custom header parser: `1_inspect_edf_voila.ipynb`, `1_inspect_edf.ipynb`, `1_inspect_edf_perdataset.py`, `1_inspect_edf_perparticipant.py`, and `2_select&remap_channels_edf(_voila).ipynb`.

**Checks performed for EEG, EOG, and ECG channels:**
- Channel configuration and montage consistency across participants
- Sampling frequency consistency
- Filter settings consistency
- Signal units
- Inverted polarity (physical_min > physical_max)
- Signal clipping (dynamic range ≤ 500 µV)
- Poor resolution (dynamic range ≥ 0.1 µV per digital unit)

**Anonymization check** (section 1.3 in Voila / section 1.4 in Jupyter): inspects the EDF+ *Local Patient ID* field (80-byte header) to detect non-anonymized patient names. The EDF+ format encodes this field as `code sex birthdate name` (space-separated); Compumedics writes the name as `LASTNAME_FIRSTNAME` and replaces it with `X_X` on anonymized export. The check isolates the name sub-field (4th token onward), strips placeholder characters (`X`, `x`, `_`, `,`, `;`, whitespace), and flags the file if anything remains. The `,` and `;` separators are stripped so that headers anonymized by *other* systems using a `Lastname,Firstname` placeholder (e.g. `Xxxxxxx,Xxxx`) are still recognized as anonymized, not just the Compumedics `X_X` form. Additionally, each non-placeholder name token (≥ 3 characters, to avoid false positives from short codes or initials) is searched case-insensitively in the file stem to detect PII leaking into the file name. Two warning levels:
- `PII in header AND file name` — real name found in both header and file name.
- `header NOT anonymized (file name looks clean)` — real name in header but file name appears clean; the important case where the file was renamed but the header was forgotten.
Files where the check cannot be performed (read failures) are already captured in `failed_edf_read.tsv`. The `patient_name` field (raw name sub-field, before cleaning) is also stored as a column in `FULL_summary_table_edf.tsv` for quick cross-reference.

**Outputs** (all written to `<study_folder>/summary_inspection/`, created on first run with a `README.md`):
- `FULL_summary_table_edf.tsv` — all parameters for all channels/files, including `patient_name` column (raw name sub-field from the EDF+ `patient_id` field)
- `anonymization_check_edf.tsv` — per-file anonymization status; columns: `subject`, `path`, `patient_id`, `name_subfield`, `patient_id_format`, `header_anonymized`, `name_in_filename`, `anon_warning`
- `EEG_summary_table.tsv`, `EOG_summary_table.tsv`, `ECG_summary_table.tsv`
- `EEG_inverted_polarity_edf.tsv`, `EEG_bad_dynamic_range_edf.tsv`, `EEG_bad_resolution_edf.tsv`
- `EDF_inspection_report.html`, `EDF_perParticipant_report.html`
- `failed_edf_read.tsv` — files that could not be read
- `README.md` — describes each output file and which tool generates it

**Skip + incremental workflow (`1_inspect_edf_voila.ipynb` and `1_inspect_edf.ipynb`)**: the Voila inspector decouples folder selection from running. Selecting the study folder only refreshes an info line ("N / M EDF file(s) already inspected" — counted against `FULL_summary_table_edf.tsv`); the scan runs on an explicit **Run inspection** button. A **"Skip files already inspected"** checkbox (checked by default) limits the run to EDF files not already in `FULL_summary_table_edf.tsv` (matched on the `path` column, normalized with `os.path.normcase`). All output tables use **merge/replace** semantics instead of being overwritten: rows for files processed this run replace their previous rows and rows for other files are kept, so the tables stay the full cumulative dataset across runs (`FULL_summary_table_edf.tsv` and `failed_edf_read.tsv` merge on file path; the per-channel summary, polarity, dynamic-range, resolution and missing/suspect tables merge on `subject`). The in-notebook **display** (file/EEG counts, channel/sampling-frequency/unit configs, polarity/range/resolution sections) reflects only the files read this run, so a full re-run (skip off) shows everything while an incremental run shows only the new files. TSV outputs are written with `index=False`. The batch `.py` scripts are unchanged.

The Jupyter notebook `1_inspect_edf.ipynb` offers the same skip + cumulative-merge behaviour, adapted to its sequential cell-by-cell model: the **"Skip files already inspected"** checkbox and the "N / M already inspected" info line live in the folder-selection cell — running the scan cell is the "Run" action, there is no separate Run button. All output tables use the same merge/replace semantics through a shared `_merge_save()` helper (`subject` key for the per-channel tables; `path` key, normcase-matched, for `FULL_summary_table_edf.tsv` and `failed_edf_read.tsv`), written with `index=False`. The `EOG_suspect_edf.tsv` / `ECG_suspect_edf.tsv` tables were renamed from `.csv` for consistency (in both the Jupyter and Voila versions). When skip is ON, the in-notebook displays and the group/session inference reflect only the files read this run, while the saved tables stay cumulative.

### 1bis. EDF Anonymizer (`1bis_anonymize_edf_voila.ipynb`, `1bis_anonymize_edf.ipynb`)

Writes **header-anonymized copies** of EDF files in batch, so a non-anonymized dataset can be cleaned **without re-exporting** from the acquisition software (which is slow). It is the write-side companion to the EDF Inspector's anonymization check (section 1): the Inspector *detects* non-anonymized headers, this tool *fixes* them.

**Safety — signal integrity guaranteed**: the tool never modifies an original file. For each file it does `shutil.copy2(src, dst)` then overwrites **only** the fixed-width identity fields of the 256-byte EDF general header. Every byte from offset 256 onward (per-channel headers + all signal data records) is left untouched. This is verified per file by comparing `sha256(file[256:])` of source vs. output **and** the total file size; both are recorded in the log (`signal_identical`, `size_identical`). A re-read of the output confirms the header now passes the anonymization check (`verified_anonymized`).

**Fields rewritten** (matching the Compumedics anonymized-export format exactly, verified against real anonymized files in the dataset):
- `patient_id` (bytes 8–88) → `X X 30-DEC-1899 X_X` — removes both the name (`X_X`) and the birthdate (replaced by the `30-DEC-1899` placeholder); code and sex become `X`.
- `recording_id` (bytes 88–168) → the real `Startdate <dd-MMM-yyyy>` token is **kept**, only the trailing admin-code / technician / equipment fields are blanked to `X X X`. If no `Startdate` token is present, the date is derived from the `start_date` header field. Toggled by a checkbox (on by default).
- `start_date` / `start_time` (bytes 168–184) → **left untouched**. Verified finding: Compumedics anonymization keeps the real recording night (it is not direct PII, and downstream sleep tools rely on it); only the birthdate inside `patient_id` is anonymized.

**Filename anonymization**: when the patient name leaks into the file name (the `name_in_filename` case from the Inspector check), the tool proposes a new file name with the name token(s) (≥ 3 chars) stripped (e.g. `01016_DUPONT_N2` → `01016_N2`). The suggestion is editable per file in the review table; files whose name is already clean keep their name.

**Companion files**: files sharing the EDF stem at a separator boundary (next char after the stem is `.`/`_`/`-`/space — so `73` does not match `731_...`) are copied and renamed to the new stem. Their **content is copied unchanged**: verified that Compumedics companions carry no patient name (`.edf.XML` is a `CMPStudyConfig` of scoring/montage settings; `*_event_xml.csv` and `_Hypnogram_*.txt` are event/stage data only). A warning reminds the user to re-check their own export profile.

**Workflow (Voila)**:
1. Select the data folder → every EDF is scanned (headers only) and classified; an info line reports `N total, M not anonymized, K already anonymized`.
2. Review table — one row per file to process: an include checkbox (pre-ticked for non-anonymized files), an editable new-name field, and a colour-coded badge (`name in header AND file name` / `name in header only` / `already anonymized`). Options: anonymize `recording_id` trailing fields (on), skip files already in the output folder / recompute everything (on), also list already-anonymized files (off).
3. Run → anonymized copies are written and the log is updated.

**Outputs** (under `<study_folder>/anonymized/`, mirroring the EDF sub-folder tree):
- `<subtree>/<new_stem>.edf` — header-anonymized copy (signal bytes identical to the original)
- `<subtree>/<companions>` — hypnogram / `.edf.XML` / event companions, copied + renamed, content unchanged
- `anonymization_log.tsv` — one row per processed file; columns: `original_path`, `anonymized_path`, `renamed`, `patient_id_before`, `patient_id_after`, `recording_id_before`, `recording_id_after`, `name_subfield`, `companions_copied`, `signal_identical`, `size_identical`, `verified_anonymized`, `status`. Merged across runs on the normalized `original_path`.
- `anonymization_failed.tsv` — files that could not be anonymized (written only if any failed)
- `FILES_DESCRIPTION.md` — describes the output files

**Important**: originals are never modified. The tool explicitly instructs the user to **delete the original non-anonymized files themselves** after verifying the `anonymized/` folder.

### 2. Channel Selection & Remapping (`2_select&remap_channels_edf_voila.ipynb`)

Interactive tool to select channels of interest and harmonize their labels across a heterogeneous dataset. Produces a JSON remapping configuration consumed by downstream analysis tools.

**Channel detection robustness**: The EDF scan section uses a 3-condition mask to identify EEG/EOG channels:
1. `transducer_type` contains `EEG`, `AGAGCL ELECTRODE`, or `EOG` (standard acquisition systems)
2. Channel name contains `EOG` (fallback for EOG channels)
3. Channel name matches `KNOWN_EEG_CHANNEL_RE` — an anchored regex covering the full 10-10 system, mastoids (M1/M2, A1/A2), and common EOG labels (LOC, ROC, E1/E2)

Condition 3 is required for EDFs exported by `mne.export.export_raw()`, which writes an empty `transducer_type` field.

**Section 4 — Define re-reference method**: Section 4 groups configurations by their post-remap canonical channel set (the channels resulting from the Section 3 harmonization), so a re-reference method is defined **once per unique harmonised montage** rather than once per raw configuration. Original configurations that become identical after remapping share a single panel (its title lists the original configs it covers, e.g. `config. 1 (n=26) + config. 3 (n=4)`). On save, the chosen method is fanned out to every original configuration in the group, so `reref_plan_by_config` stays keyed by the original config label and Sections 5/6 are unaffected.

**Section 5 — Preview & save JSON**: section 5 exposes a single "Preview & Save" button. Clicking it builds the per-participant dict for the participants configured this session and **merges** it into any existing `<data_folder>/config_param/remap_reref_persubject.json` (entries for re-configured participants are replaced, all others kept), so the file stays the full cumulative configuration. The saved file is sorted by participant id; the on-screen preview lists this session's participants first, then the previous ones, with a note of how many were added/updated and the new total. There is no separate save step. `mne_reref_plan.json` is no longer generated (it was redundant with `remap_reref_persubject.json` which already carries `ref_channels` per participant).

**Skip + incremental workflow**: like the EDF inspector, folder selection is decoupled from running. Selecting the folder refreshes an info line ("N / M participant(s) already configured" — counted against the existing `remap_reref_persubject.json`); the scan runs on an explicit **Run scan** button. A **"Skip participants already configured"** checkbox (checked by default) excludes participants already in the JSON from the scan, so the configurations and every downstream section involve only the new participants, whose entries are merged into the JSON on save (see Section 5). Participant ids and file paths are compared with `os.path.normcase` (case/separator-insensitive). `failed_edf_read.tsv` is merged on file path the same way, so it reflects the current config state rather than only the last scan.

**Robust JSON loading** (see *Cross-cutting procedures*): every read of `remap_reref_persubject.json` (info line, scan filter, save-merge, section 6 test) goes through the shared lenient loader.

**Section 6 — Test the JSON**: applies each participant's remap + re-reference and reports the resulting channel configurations (harmonization succeeds when a single configuration remains). A **scope** toggle selects what to test — **Whole database** (default; every participant in the JSON that has an EDF in the folder, normcase-matched) or **New files (this session)**. The toggle and run button persist in their own area with results rendered below, so the test can be re-run with a changed scope (e.g. verify the just-modified participants, then the whole database) without redoing the workflow.

**Jupyter twin parity (`2_select&remap_channels_edf.ipynb`)**: the Jupyter version implements the same skip + cumulative-merge workflow as the Voila, adapted to its sequential model. The **"Skip participants already configured"** checkbox and the "N / M already configured" info line live in the folder-selection cell (running the scan cell is the "Run" action). Save (Section 5) now **merges** this session's entries into the existing `remap_reref_persubject.json` instead of overwriting it, and the preview shows the merged result. Section 6 gains a **scope toggle** (`New files (this session)` / `Whole database`) placed in its own cell just above the test cell — change the toggle then re-run the test cell. All JSON reads use the same lenient loader as the Voila.

### 3. Hypnogram Label Remapping (`3_remap_hypno_voila.ipynb`, `3_remap_hypno.ipynb`)

Interactive tool to harmonize sleep stage labels across a heterogeneous database, converting different scoring conventions (e.g. `0,1,2,3,4` or `W,S1,S2,S3,S4`) to the standard AASM format (`W`, `N1`, `N2`, `N3`, `R`).

**Workflow (5 sections):**
1. **Scan** — Select data folder, hypnogram suffix, and output suffix → auto-detects files recursively; reports unique label configurations and flags problematic epochs in a **single combined message** — mid-recording `?` epochs and suspicious labels (any label not in `DEFAULT_MAPPING`, e.g. `U`/`M`) are listed together, one aligned row per label (`?` treated as a label), each showing the file count, epoch count, and the **affected file names** (the quoted label is left-padded so the columns align for multi-character labels). Boundary `?` epochs (first/last ~10) are not flagged here — they are mapped in Section 3. Also reports: how many remapped files already exist with the output suffix; a list of `.txt` files matching neither suffix (informational, only when both suffixes are defined); automatically exports `mid_uncertain_epochs_to_verify.tsv` if mid-recording `?` epochs or unexpected labels are found. If the TSV already exists (e.g. edited manually between sessions), it is loaded automatically and its corrections are pre-applied in memory; the flags are computed on the **original** labels (before applying the TSV), so the pre-corrected epochs still appear in Section 2 (pre-filled) instead of disappearing — the TSV file is never overwritten at scan time. Optional checkbox to exclude participants whose remapped file already exists from all downstream processing.
2. **Uncertain and unexpected epoch review** — A **single review widget** (shown whenever any mid-recording `?` or unexpected label is found) navigates all flagged epochs in one flat list, handling both kinds together: mid-recording `?` epochs and epochs carrying a label not in `DEFAULT_MAPPING` (e.g. `U`, `M`).
   - **Show** filter dropdown: *All issues* / *Mid-recording `?`* / one entry per unexpected label, to focus the navigation.
   - ±N context epochs (adjustable slider, default 5); the current epoch in red, already-corrected epochs in green ✓, other still-pending flagged epochs in orange.
   - **Apply to this epoch** (assign + auto-advance) and **Apply to all in this file** (batch-replace every flagged occurrence of the current label within the current participant — now works for mid-`?` too).
   - Combobox pre-filled from `DEFAULT_MAPPING` **only for unexpected labels**; mid-recording `?` are left blank (a mid-night `?` is unlikely to be Wake, so the user picks explicitly from the full list).
   - The confirmation line shows `old → new` using the original flagged label (correct even when the epoch was pre-loaded from an existing TSV, where `hypno_data` already holds the new value).
   - Corrections update `hypno_data` in memory and a single internal `STATE.corrections` dict; fully corrected labels disappear from Section 3 configs on next run. Clicking "Corrections done" persists all corrections into `mid_uncertain_epochs_to_verify.tsv`.
3. **Remap labels** — Per-configuration accordion widget with combobox suggestions pre-filled from `DEFAULT_MAPPING`; suspicious labels are highlighted in red with inline epoch context; warns if the final mapping leaves non-AASM labels; confirmation required before proceeding
4. **Save** — Writes remapped hypnograms next to originals using the output suffix defined in Section 1; end message confirms completion and recalls the suffix used
5. **Verify** — Before/after configuration summary; verdict fails only if non-AASM labels remain (multiple configurations with valid AASM labels are acceptable — e.g. insomnia patients legitimately missing N3)

**Custom (non-AASM) stages** (see *Cross-cutting procedures*): tool 3 is the **only** writer of `config_param/custom_stages.json`. A `Custom stages` field (Section 1, comma-separated, auto-filled from any existing JSON) lists labels deliberately kept outside the AASM set; `current_acceptable()` = `STANDARD_LABELS | {MT} | <field>`, so those labels no longer trip the Section 5 verdict. Section 3's **Save remapping** detects non-AASM *target* labels and offers a **➕ Register** button that appends them to the field; Section 4's **Save files** then **merges** the declared stages that actually survive in the remapped output into `custom_stages.json`. Downstream, tools 5/6/7/8 auto-load this file so the kept labels are recognised (hypnospectrogram, per-stage tables, rejection) instead of being flagged as unrecognised.

**Declared custom stages are first-class, not errors** (UX): a label listed in the Section 1 `Custom stages` field is no longer treated as suspect/unexpected. **Section 1** reports it on its own info line (`'M' : N file(s), K epoch(s) — ids…`) and excludes it from the **Section 2** review widget and the `mid_uncertain_epochs_to_verify.tsv`. **Section 3** suggests the *raw* custom label as its own target (identity, "keep" — the user may still rename it), shown in blue with a "rename if needed" note rather than suspect-red, and offers it in the combobox options. **All post-remap reporting keys off the labels actually present in the OUTPUT, not the raw field**, so renaming e.g. `M→SD` reports and saves only `SD` (never the now-unused `M`): Section 3's save splits output non-AASM targets into *registered* (already in the field → green "kept" line) vs *unregistered* (→ the ➕ Register warning); **Section 4** shows declared custom stages as a green **info box** and only **genuine** non-AASM, non-custom labels block saving / require the confirm-checkbox (MT is treated as acceptable); the **conclusion** notes that tools 5/6/7/8 will recognise the kept stages; **Section 5**'s success line lists the custom stages actually present in the reloaded files. The Section 1 field is never auto-pruned of renamed-away labels (they may still be in use by another configuration), and `save_custom_stages` filters to output-present stages, so `custom_stages.json` stays correct regardless.

**Hypnogram suffix auto-detection** (see *Cross-cutting procedures*): tool 3 auto-fills the `Hypnogram suffix:` widget with the **shortest** candidate suffix on ties — the goal is the raw (unremapped) hypnogram, not an already-processed one (the reverse of tools 5/6, which prefer the longest).

**Key constants:**
- `DEFAULT_MAPPING`: `0→W`, `1→N1`, `2→N2`, `3→N3`, `4→N3`, `5→R`, `?→W`, `S1→N1`…
- `STANDARD_LABELS`: `{W, N1, N2, N3, R}`
- `ACCEPTABLE_LABELS`: `STANDARD_LABELS | {MT}` — used for AASM compliance warnings (MT = movement time is tolerated)

**Outputs**:
- One `.txt` file per participant with the output suffix (e.g. `_Hypnogram_remapped.txt`), one label per line
- `mid_uncertain_epochs_to_verify.tsv` — written to `<data_folder>/` at scan time when mid-recording `?` epochs or unexpected labels are found; columns: `participant_id`, `epoch_index`, `epoch_time_sec`, `total_epochs`, `original_label` (`?` for mid-recording unscored epochs, or the raw unexpected label e.g. `U`, `M`), `context` (±5 epochs), `corrected_label`; updated with all corrections after "Corrections done" is clicked. If the file already exists at scan time it is loaded and applied in memory instead of being overwritten; the pre-loaded corrections are shown in a summary panel (as `ep.N: old→new`) and also appear **pre-filled** in the Section 2 review widget (flags are computed on the original labels, so they are not hidden). Corrections whose `original_label` no longer matches the current hypnogram value (re-scored since the TSV was written) are flagged as conflicts and shown in a warning panel in Section 2.

- `config_param/custom_stages.json` — written/merged when the user keeps non-AASM labels as custom stages (see *Custom (non-AASM) stages* above); a flat `{"custom_stages": [...]}` list consumed by tools 5/6/7/8.

`check_hypno_config.py` is the legacy script that preceded this notebook; kept for reference.

### 4. Event Label Harmonization (`4_remap_events_edf_voila.ipynb`, `4_remap_events_edf.ipynb`)

Interactive tool to **visualize** the scored-event configurations present across a heterogeneous database and **harmonize** their raw labels to a single canonical vocabulary — the event analogue of tool #2 (channel selection & remapping). Scored events are annotated during sleep scoring and exported by Profusion/Compumedics (apnea, hypopnea, arousals, limb movements, PLM, SpO2 desaturation…).

**Event sourcing (CSV-first / XML-fallback)**: a shared `load_events(edf_path, csv_suffix='_event_xml.csv')` helper reads the Compumedics event CSV (`Name, Start, Duration` in seconds) first and, when it is absent, parses the `<ScoredEvents>` of the `*.edf.XML` (Profusion `CMPStudyConfig`; `<Input>` is ignored). Both sources were verified equivalent on ICEBERG. The `<ScoredEventSettings>` catalogue (the profile's possible event types) is **not** used for grouping.

**Configurable event-CSV suffix**: the Compumedics event-CSV name is no longer hardcoded. A **`CSV suffix:`** text field (default `_event_xml.csv`) drives `event_companion_paths(edf_path, csv_suffix)` / `load_events(edf_path, csv_suffix)`, so datasets exported with a different CSV suffix are supported. On folder selection the suffix is **auto-detected** (mirroring the hypnogram-suffix detection of `5_quality_overview`, see *Cross-cutting procedures*): the `.csv` files next to each EDF are scanned, candidate suffixes counted, and the field auto-filled with the **most frequent** suffix (shortest on ties — events have no "more specific remapped" variant, unlike hypnograms), with a colour-coded info line (green = all EDFs matched, orange = partial/none). The chosen suffix feeds both the Section 1 scan and the Section 1bis CSV-vs-XML check; the XML fallback (`.edf.XML`) is unchanged.

**Configuration grouping**: files are grouped by their `frozenset` of **unique event labels actually present** — two files with the same unique labels share one configuration even if their event counts/timing differ.

**Workflow (sections):**
1. **Scan** — select the data folder (recursive `rglob('*.edf')`); selecting the folder only refreshes an info line **and auto-detects the event-CSV suffix** (editable `CSV suffix:` field, colour-coded detection line — see *Configurable event-CSV suffix* above), the scan runs on an explicit **Run scan** button. A **"Skip labels already mapped"** checkbox (on by default) hides labels already present in an existing `event_remap.json` (incremental harmonization when a new cohort is added).
1bis. **(Optional) CSV vs XML consistency check** — opt-in button; for every file having both companions, compares the CSV and the XML `<ScoredEvents>` (name + start + duration, tolerance 1e-3 s) and writes `event_source_mismatch.tsv`. Opt-in because it forces reading both files for every EDF.
2. **Visualize configurations** — because two files sharing the same label names can still form distinct configs (a config = the exact set of labels *present*, so a missing label splits it off), the configs are not stacked: a **dropdown** ("Show config:") selects one configuration to detail (its sorted unique labels + file/label counts), and a **"Show file ids" toggle button** (replacing the old `<details>` arrow) shows/hides that config's file-id list in a scrollable box. The global table of every raw label (file count + total occurrences + suggested canonical) is kept below.
3. **Harmonize labels** — one editable row per unique raw label (combobox pre-filled from `DEFAULT_EVENT_MAPPING`, free text allowed), with an **ignore** toggle (stored as `null`); filtered by the skip checkbox. Each row is a bordered, column-aligned line for readability (the row layout forces `overflow='hidden'` so the jupyter-widgets default `.widget-box { overflow:auto }` does not raise a stray per-row scrollbar — see CLAUDE.md). A **"Validate mapping & ignores"** button summarizes the choices (N mapped / N ignored / N left empty, warning on empties) and **unlocks** the Section 4 save button (which starts disabled). Editing any row after validating (or re-running the scan / toggling the skip checkbox) re-locks Section 4 and clears the previous save preview, so the saved JSON always reflects the latest Section 3 selection.
4. **Preview & save** — enabled only after Section 3 validation; builds a flat `{raw_label: canonical_label}` mapping and **merges** it into `config_param/event_remap.json` via the lenient JSON loader (labels mapped this session replace their old value, all others kept; keys sorted). Unmapped non-ignored labels are reported and not saved.
5. **Verify** — applies the saved mapping to every configuration, reports the resulting harmonized labels, and passes when no raw label is left unmapped (ignored labels count as handled). A scope dropdown can restrict the view to configs with unmapped labels.

**`DEFAULT_EVENT_MAPPING`** (editable suggestions, snake_case canonical vocabulary): apnea subtypes kept (`apnea_obstructive` / `apnea_central` / `apnea_mixed`), `hypopnea`, `spo2_desaturation`, arousal subtypes kept (`arousal_respiratory` / `arousal_spontaneous` / `arousal_limb` / `arousal`), limb laterality collapsed (`limb_movement`, `plm`). `suggest_canonical()` also tolerates a trailing `(Left)`/`(Right)` marker.

**Outputs** (under `<data_folder>/config_param/`):
- `event_remap.json` — global flat `{raw_label: canonical_label}` (`null` = ignore), merged across runs
- `event_source_mismatch.tsv` — only if the CSV-vs-XML check is run; columns `file_id, n_csv, n_xml, n_only_in_csv, n_only_in_xml, labels_only_in_csv, labels_only_in_xml, status`
- `failed_event_read.tsv` — files with no readable event companion (only if any failed)

### 5. Quality overview (`5_quality_overview_voila.ipynb`)

*(Stable tool — formerly “Phase 1” of the preprocessing pipeline.)*
Implemented as `tools/5_quality_overview_voila.ipynb`. Produces one `mne.Report` HTML per participant. For each EEG channel: signal amplitude histogram with Savitzky-Golay smooth + peak detection, time series, metrics table, and a YASA hypnospectrogram (0.1–40 Hz bandpass applied per-channel just before plotting). Flags suspect channels for priority inspection. At the end of each run, generates `dataset_overview.html` — a single-page dataset-level summary with statistics and distribution plots per electrode, consumed by Phase 2 to identify channels to exclude.

**EDF scan**: recursive (`rglob('*.edf')`), so datasets organized in subfolders (e.g. `group1/`, `group2/`) are fully covered without needing to run the tool per subfolder.

**Hypnogram suffix auto-detection** (see *Cross-cutting procedures*): tools 5 and 6 prefer the **longest** suffix among candidates appearing for ≥ 50% of the maximum count (more specific = remapped/processed version; count breaks ties on equal length). The 50% threshold prevents a rare accidental long suffix from winning when remapping is far from complete.

**Live output warnings**: if a hypnogram is not found, fails to load, has a length mismatch with the EDF, or contains unrecognised stage labels, a plain-text `⚠` warning is printed in the notebook output area immediately after the per-participant result line (in addition to the yellow banner already shown in the HTML report's Overview section).

**Hypnogram label validation**: after mapping labels to YASA integers (`W→0`, `N1→1`, `N2→2`, `N3→3`, `R→4`), the tool checks for unrecognised labels. `MT` (movement time) is silently tolerated — YASA treats it as an artifact epoch (NaN). Any other unrecognised label triggers a warning. Two severity levels: if > 10% of epochs are unrecognised, `hypno_vec` is set to `None` (spectrogram skipped, error-level warning pointing to `3_remap_hypno_voila`); if ≤ 10%, `hypno_vec` is kept but a warning lists the unrecognised labels and their count. This catches hypnograms that were not yet remapped to AASM convention (e.g. `S1/S2/S3/S4` or raw numeric labels).

**Custom (non-AASM) stages** (see *Cross-cutting procedures*): an editable `Custom stages` field (auto-filled from `config_param/custom_stages.json`) registers project-specific labels. Registered labels are treated as **recognised** — they no longer count toward the >10% skip and are not set to `None` — the YASA hypnospectrogram is replaced by the shared `plot_hypnospectrogram()` (custom labels stacked below N3, their own colours), and `quality_summary_by_stage.tsv` + the `dataset_overview.html` by-stage tables/boxplot insets iterate `['W','N1','N2','N3','R'] + custom_stages` (`stage_map` extends the YASA codes with `5+i`).

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

`flat_pct` = fraction of consecutive sample pairs with \|diff\| < `max(2×ADC_step, 0.06 µV)`. `bounds_pct` = fraction of samples within 0.5 µV of the EDF physical_min/max header limits, reconstructed in µV by `get_phys_bounds_uV()`. **Unit handling**: MNE keeps `_raw_extras['physical_max']`/`['offsets']` in each channel's *native* EDF physical unit (µV, mV, V…), not in µV, whereas the signal is read as `raw.get_data() * 1e6` (always µV). `get_phys_bounds_uV()` therefore multiplies both bounds by `extras['units'][ch_idx] * 1e6` (MNE's native-unit→volts factor: `1e-6` µV, `1e-3` mV, `1.0` V) so the comparison is unit-consistent. This is essential for non-EEG channels: Compumedics EOG/EMG/ECG are declared in **mV** (e.g. `physical_max = 1.0 mV`), so without the conversion they are compared against a 1.0 µV bound — 1000× too small — and `bounds_pct` flags ~98–100 % of a perfectly healthy signal. EEG channels (declared in µV) are unaffected. `hist_extreme_pct` = fraction of samples in the outermost histogram bins. **Kurtosis is intentionally NOT used as a flagging criterion** — normal PSG EEG has physiologically high kurtosis (spindles, K-complexes produce values of 100–500), making it unreliable without per-subject normalization. `p99_abs_uV` and `p999_abs_uV` (99th and 99.9th percentile of |amplitude|) are recorded as informational metrics in `quality_summary.tsv` but are not currently used for flagging; they are useful for cross-channel and cross-dataset amplitude comparison.

**`dataset_overview.html` — dataset-level summary**: generated at the end of every run from the full cumulative `quality_summary.tsv` (reflects all participants processed to date, not just the current run). Contains two levels:
- **Global section (all electrodes pooled)**: stats table (mean / median / p5 / p25 / p75 / p95 per metric), **mean (median) by sleep stage table** (metrics as rows — `mean_uV`, `std_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `p99_abs_uV`, `p999_abs_uV` — stages W/N1/N2/N3/R as columns; shown only when `quality_summary_by_stage.tsv` is present), n_peaks frequency table by electrode (flags DC drift and quantization cases), pooled boxplots (one subplot per key metric; each subplot contains a **stage inset** in its upper-right corner showing per-stage boxplots with stage colours W/N1/N2/N3/R — inset shown only when stage data is available; stage colours: W `#969696`, N1 `#9e9ac8`, N2 `#807dba`, N3 `#6a51a3`, R `#c994c7`), grouped boxplots (one subplot per key metric, x-axis = electrode — compares electrodes side by side).
- **Per-electrode sections**: stats table for that electrode only, boxplots of each metric's distribution across participants with individual data points overlaid (reveals outlier participants for that channel).

Key metrics shown in plots: `std_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `p99_abs_uV`, `p999_abs_uV`. All numeric metrics (`mean_uV`, `kurtosis`, `skewness` included) appear in the stats tables. HTML is standalone (figures embedded as base64 PNG, no external dependencies).

**Outputs per run**:
- `<data_folder>/reports_quality_overview/<relative_subfolder>/<file_id>_quality_overview.html` — HTML reports mirror the EDF subfolder structure under `reports_quality_overview/`
- `<data_folder>/reports_quality_overview/quality_summary.tsv` — all numeric metrics for all channels, cumulative across runs (always at root); columns: `file_id`, `channel`, `mean_uV`, `std_uV`, `kurtosis`, `skewness`, `p99_abs_uV`, `p999_abs_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `n_peaks`, `suspect_reason`, `exclude` (last two columns)
- `<data_folder>/reports_quality_overview/dataset_overview.html` — dataset-level statistics and plots (always at root, regenerated each run)
- `<data_folder>/reports_quality_overview/quality_summary_by_stage.tsv` — key metrics split by sleep stage, one row per `file_id × channel × stage`; columns: `file_id`, `channel`, `stage`, `mean_uV`, `std_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `p99_abs_uV`, `p999_abs_uV`. Populated only for participants with a valid hypnogram. Cumulative merge: rows for `file_id ∈ attempted_ids` are replaced, all others kept (same as `quality_summary.tsv`).
- `<data_folder>/reports_quality_overview/failed_files.tsv` — files that could not be read (at root)

**End-of-run summary** (printed in the notebook output): participants processed, participants with ≥1 flagged channel, total flagged channels, files failed to load, path to `dataset_overview.html`.

### 6. Preprocessing + epoch rejection (`6_preprocessing_voila.ipynb`)

*(Stable tool — formerly “Phase 2” of the preprocessing pipeline.)*

Implemented as a Voila notebook with four sections: (1) path configuration, (2) preprocessing and rejection parameters, (3) participant selection, (4) processing loop. **Section 2 is organised into two headed sub-sections**: **Preprocessing** (resampling + bandpass filter) and **Epoch rejection** (peak-to-peak amplitude per stage, flat signal & gradient, 1/f fit quality, and the optional event-based rejection — see below).

**Inputs**:
- `quality_summary.tsv` from Phase 1 — `exclude` column identifies channels to drop before preprocessing
- `remap_reref_persubject.json` from `2_select&remap_channels_edf` — drives channel remapping and re-referencing per participant
- Raw EDF files and remapped hypnograms (default suffix `_Hypnogram_remapped.txt`)
- *(optional, for event-based rejection)* `config_param/event_remap.json` from `4_remap_events_edf` and the per-EDF scored-event companions (`*_event_xml.csv` / `*.edf.XML`). Section 1 has an explicit `event_remap.json` `FileChooser` (auto-pointed at `<edf_folder>/config_param/` when present) and an editable **`Event CSV suffix:`** field auto-detected from the `.csv` companions next to the EDFs (most frequent suffix, shortest on ties; colour-coded info line), mirroring the suffix auto-detection of tools 4 / 5.

**Channel-name handling when loading EDF data from notebook outputs** (critical):

The EDF files on disk carry their **original** acquisition channel names (e.g. `Fp1, C3, O1, A2`), usually alongside non-montage channels (ECG, respiration, SpO2, position) that may be sampled at a **higher rate** than the EEG (e.g. ECG at 512 Hz while EEG is at 256 Hz). The config/report files refer to channels by their **harmonized (remapped)** names:
- `remap_reref_persubject.json` stores the mapping original → remapped, e.g. `{"Fp1":"F","C3":"C","O1":"O","A2":"M"}`. Its **keys are the original EDF names**, its values the harmonized names.
- `quality_summary.tsv` lists each channel under its **remapped** name, because `5_quality_overview_voila` renames the channels (`raw.rename_channels`) *before* computing per-channel metrics. So the `exclude` flags and the channel checkboxes derived from it are expressed in remapped names.

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
- **No inverse-remap dictionary**: the include list is `list(remap.keys())` directly, identical to `5_quality_overview_voila` — the two tools stay consistent.

The two tools then differ only in what follows:
- `5_quality_overview_voila` analyses **every** remap channel, so it keeps them all (uses `preload=True`, no further selection).
- `6_preprocessing_voila` lets the user deselect channels in the UI (selection in **remapped** names). After the rename it drops the de-selected channels, then calls `raw.load_data()` — so the signal is read from disk **only for the channels actually kept**:
  ```python
  present = [ch for ch in selected_channels if ch in raw.ch_names]   # remapped namespace, post-rename
  raw.drop_channels([ch for ch in raw.ch_names if ch not in present])
  raw.load_data()
  ```
  A `present`-empty guard (e.g. if the rename failed) marks the participant as failed and skips it, so a single bad file never crashes the run.

**Shared utility functions** (`drop_suffix_duplicates(raw)` → `(raw, dropped_list)`, and `adapt_remap_dict_to_suffixes(raw, remap_dict)`): defined identically in `5_quality_overview_voila` and `6_preprocessing_voila`, used right after `read_raw_edf` — see *Cross-cutting procedures → MNE EDF signal loading pattern*.

**Preprocessing steps** (applied in this order, each optional via widget):
1. **Resampling** — `raw.resample(target_freq, npad='auto')`. Target frequency chosen by user; applied before filtering to avoid aliasing. Step is skipped if checkbox is unchecked.
2. **Re-referencing** — applied as specified in JSON config per participant: `'average'` → common average reference; `[list]` → subtract listed channel(s) then drop them; empty → no re-referencing.
3. **Bandpass filter** — FIR zero-double-pass Hamming window, defaults `l_freq=0.1 Hz, h_freq=40 Hz`. Applied via `raw.filter(..., method='fir', phase='zero-double', fir_window='hamming', fir_design='firwin')`.

**Epoching**: 30-second fixed-length epochs created with `mne.make_fixed_length_epochs(raw, duration=30)`. Sleep stage assigned to each epoch from the hypnogram; epochs at the tail beyond the hypnogram length are discarded.

**Epoch rejection — six methods** (the sixth, *Event overlap*, is optional and off by default):

All methods operate on the raw epoch data in µV (`epochs.get_data() * 1e6`, shape `n_epochs × n_channels × n_times`). Rejection masks are boolean arrays of shape `(n_epochs, n_channels)` — a `True` entry means that (epoch, channel) pair was flagged. An epoch is considered **rejected** if *any* channel is flagged by *any* method. A single ordered registry — `METHOD_ORDER = ['amplitude', 'flat', 'gradient', '1f_error', '1f_r2', 'event']` with `METHOD_CODE`/`MULTIPLE_CODE` — is the one source of truth shared by the mask builder, heatmap, per-epoch log, per-stage summary and global table, so `event` appears in every output **only when event rejection actually ran** and older event-free outputs keep their original column set.

| Method | Signal feature | Default threshold | Notes |
|--------|---------------|-------------------|-------|
| **Amplitude** | Peak-to-peak = `max(epoch) − min(epoch)` | W: 300, N1: 250, N2/N3: 200, REM: 250 µV | Per-stage threshold; W/REM more lenient because muscle and eye-movement artefacts are physiologically common in those stages. Equivalent to MNE's `drop_bad(reject=...)` criterion. |
| **Flat signal** | Peak-to-peak < threshold | 1 µV | Detects disconnected electrodes or amplifier saturation within a single epoch. Logically identical to MNE's `drop_bad(flat=...)` criterion: both compare `ptp` against a low-amplitude threshold. |
| **Gradient** | `max(|diff(epoch)|)` across time | 100 µV/sample | Maximum sample-to-sample absolute difference; sensitive to sudden jumps, electrode pops, and movement artefacts not captured by peak-to-peak. `diff` and `max` both operate on `axis=-1` (time axis) to handle the 3D `(n_epochs, n_channels, n_times)` array correctly. |
| **1/f fit quality** | Specparam aperiodic fit on Welch PSD (4 s windows, 2–30 Hz, `aperiodic_mode='fixed'`, `peak_width_limits=[0.5, 20]`, `min_peak_height=0.3`) | MAE > 0.15 OR R² < 0.95 | Fit restricted to ≥2 Hz to limit slow-wave influence. A **full peak model is used deliberately** — periodic components (spindles, alpha…) are modelled and removed *before* assessing the aperiodic fit quality. Forcing `max_n_peaks=0` would push all peak power into the aperiodic component, degrading R² and over-rejecting nearly every N2/REM epoch. Metrics read via `get_metrics('error','mae')` / `get_metrics('gof','squared')` (specparam 2.x). A failed fit is treated as a double flag (both error and R²). |
| **Event containment** *(optional, off by default)* | 30 s epoch **containing the onset** of any **selected** canonical scored-event type (arousal, apnea, hypopnea, limb movement, SpO2 desaturation…) | **onset-only** — the epoch holding the event `Start`; the annotated `Duration` is **intentionally ignored** (clinicians often score only the onset without a reliable duration), so each event flags exactly one epoch | **Epoch-level** flag, replicated across all channels → single `flag_event` column. Events read with the shared CSV-first / XML-fallback `load_events(edf, csv_suffix)`; raw labels mapped to canonical via `event_remap.json` (tool 4). UI: a checkbox to activate, **a wrapping row of checkboxes for the canonical types** (all shown at once, populated from the chosen `event_remap.json`), and an inline note explaining the onset-only rule so the choice is informed. A **"Count affected epochs"** button reports, over the participants currently checked in Section 3, how many epochs each selected type would flag (overall + per stage) using only the hypnogram length/stages and events — no signal is read. Missing/unreadable event companions are non-fatal (the file keeps the other 5 methods, never added to `failed`). |

**Custom (non-AASM) stages** (see *Cross-cutting procedures*): an editable `Custom stages` field (Section 1, auto-filled from `config_param/custom_stages.json`) extends the per-stage logic. Each custom stage gets its **own amplitude-threshold widget** (default 250 µV, generated dynamically when the field changes) feeding `ptp_thresholds`; the per-participant summary, `global_rejection_by_stage.tsv`, the heatmap hypnogram strip and the **"Count affected epochs"** estimate all iterate `['W','N1','N2','N3','R'] + custom_stages`. Custom-stage epochs are still rejected by the other (stage-independent) methods regardless.

**Heatmap** (rendered into `{file_id}_preprocessing_report.html`, **not** saved as a standalone PNG): channels (Y-axis) × epochs (X-axis); each cell coloured by the flagging method with priority encoding when multiple methods fire. A hypnogram strip is drawn above the main heatmap. Colour scheme: dark purple = none, red = amplitude, blue = flat, orange = gradient, yellow = 1/f error, green = 1/f R², **magenta = event**, dark red = multiple. Title includes overall rejection percentage.

**Two-step QC approach** — Phase 2 does **not** drop epochs. It saves ALL epochs (including flagged ones) with an MNE `metadata` DataFrame attached, so downstream Phase 2b can inspect rejected epochs before finalising the rejection.

**Outputs per participant** — the `.fif` + its params sidecar under `<output_folder>/derivatives/<edf_subtree>/`; the TSV/HTML reports under `<output_folder>/reports_preprocessing/`:
- `{file_id}_all-epo.fif` — all epochs with `epochs.metadata` DataFrame (columns: `epoch_idx`, `stage`, `reject_flag`, `reject_method`, `flag_amplitude`, `flag_flat`, `flag_gradient`, `flag_1f_error`, `flag_1f_r2`, plus `flag_event` **when event rejection ran**). The `flag_<method>` columns are **per-epoch "any channel" booleans** — the per-(epoch, channel) mask is **not** persisted. The per-epoch/per-stage TSVs and `global_rejection_by_stage.tsv` gain the matching `flag_event` / `event` entries the same way — additively, so event-free runs stay byte-compatible with earlier outputs.
- `{file_id}_preprocessing_params.json` — the resampling/filter settings + the per-stage rejection thresholds actually used (amplitude p-p per stage, flat, gradient, 1/f MAE/R²) + `methods_run`. Read back by **tool 7** (QC of rejected epochs) to draw threshold reference lines and recompute per-channel margins. Written non-fatally.
- `{file_id}_epoch_rejection.tsv` — per-epoch rejection table (columns: `file_id`, `epoch_idx`, `stage`, `reject_flag`, one `flag_<method>` bool per method, incl. `flag_event` when event rejection ran).
- `{file_id}_rejection_summary.tsv` — rejection counts per stage per method (columns: `file_id`, `stage`, `method`, `n_total`, `n_rejected`, `pct_rejected`).
- `{file_id}_preprocessing_report.html` — MNE HTML report with the heatmap and two rejection tables. The **per-stage rejection table** has one row per stage (W/N1/N2/N3/R + custom) and one column per method (Amplitude, Flat, Gradient, 1/f error, 1/f R², and Event when event rejection ran), each cell showing the **% in front and the raw count `(n)` in parentheses**; % is relative to the stage total (same denominator as `global_rejection_by_stage.tsv`). When event rejection is enabled, a second **"Event-based rejection by type"** table lists one row per selected canonical event type (+ a bold `(any selected)` union row) with the flagged epoch count, % of all epochs, then `%(n)` per stage — so event types that reject too many epochs can be spotted and de-selected. Both tables are report-only (HTML); no TSV schema changes. Built by `build_stage_method_html()` / `build_event_type_html()`, with per-type epoch masks computed in the run loop via `compute_event_epoch_mask(..., [t], ...)`.

**Global output** (in `<output_folder>/reports_preprocessing/`):
- `global_epoch_rejection.tsv` — concatenation of all `{file_id}_epoch_rejection.tsv` across participants
- `global_rejection_by_stage.tsv` — concatenation of all `{file_id}_rejection_summary.tsv` across participants
- `preprocessing_failed.tsv` — participants that could not be processed (EDF not found, config missing, hypno mismatch, etc.)

### 7. Interactive QC of rejected epochs — Phase 2b (`7_inspect_rejected_epochs_voila.ipynb`, `7_inspect_rejected_epochs_batch.py`, `qc_rejected_epochs_lib.py`)

Manual quality control of the epochs that `6_preprocessing_voila` flagged: inspect each rejected epoch, override the keep/reject decision, and export a validated `{file_id}_clean-epo.fif`. Reads one participant's tool-6 outputs (`{file_id}_all-epo.fif` + optional `{file_id}_preprocessing_params.json`) — **the raw EDF is never reloaded**; the analysis channel set is exactly what the `.fif` holds (EEG-only in practice). Tool-6 outputs are never modified. The Voila app is the primary delivery (a code-visible Jupyter twin is planned); the `.py` batch twin produces the Section-2 report over a whole database.

**Shared library (`qc_rejected_epochs_lib.py`)**: the analysis + plotting used by both the notebook and the batch live in one module (imported by both) to avoid drift — `METHOD_ORDER`, the heatmap/method colour palette, the custom-stage helpers, the Welch-PSD config and the specparam 1/f fit are copied from tool 6 and kept in sync. Because tool 6 does not persist a per-(epoch, channel) mask, the **per-epoch reject decision is authoritative from `epochs.metadata`**, while the **per-channel attribution** shown here (which channel drove a flag, margins to threshold) is **recomputed** from the signal with the same formulas + the persisted thresholds (fallback: tool-6 defaults when no params JSON is present).

**Section 1 — Load**: a `FileChooser` for the derivatives root lists every `*_all-epo.fif` found recursively; a participant dropdown loads one file (`mne.read_epochs`), reads the params JSON (thresholds) and auto-fills the editable `Custom stages` field. A summary reports epoch/rejection counts, channels, sfreq, and whether thresholds came from the params JSON or defaults.

**Section 2 — Global per-stage report** (Run): per sleep stage (`W/N1/N2/N3/R` + custom) — a **PSD overlay** (per-epoch mean-across-channel PSD: clean epochs grey + median/IQR band, rejected epochs coloured by their reject method), **metric distributions** (p-p, gradient, 1/f MAE, 1/f R² — clean vs rejected, with threshold lines), a **p-p vs gradient scatter** (method-coloured), and a **stage × method rejection table** (from the metadata flags). 1/f is fitted per epoch (~1 min for a full night); the editable thresholds drive only the reference lines. Assembled into an `mne.Report`; **Save** writes `{file_id}_qc2b_report.html`.

**Section 3 — Per-epoch navigator** (Run): walks the rejected epochs (filterable by method or stage). For the current epoch: a **stacked montage** (± context) with the current-epoch trace of each channel coloured by its recomputed flagging method, the gradient-max sample marked, and method-coloured channel labels; plus a **detail panel** — PSD + aperiodic fit (worst-R² channel), mean band power (δ/θ/α/σ/β) + 50 Hz ratio, an epoch spectrogram, and a per-channel metric table (value vs threshold). A **keep / reject** toggle overrides the decision in both directions — confirm a rejection or *rescue* a clean-looking flagged epoch.

**Section 4 — Manual override & save**: a review strip (hypnogram + final keep/reject per epoch, overridden epochs marked) and counts (kept / rejected / rescued / newly-rejected). **Save** writes, next to the `.fif`: `{file_id}_clean-epo.fif` (kept epochs only; metadata carries `manual_override` + `final_reject`), `{file_id}_epoch_rejection_reviewed.tsv` (per-epoch metadata + the two override columns), and `{file_id}_qc2b_review_log.tsv` (one row per overridden epoch: `epoch_idx`, `stage`, `orig_reject`, `final_reject`, `action` ∈ rescued/added). Override granularity is **whole-epoch** (only 3–4 EEG channels, and `clean-epo.fif` drops whole epochs anyway).

**Batch twin (`7_inspect_rejected_epochs_batch.py`)**: `python 7_inspect_rejected_epochs_batch.py <derivatives_root> [--no-1f] [--limit N] [--out DIR]`. Runs Section 2 for every participant (writing each `{file_id}_qc2b_report.html`) plus a **database-level aggregate** — rejection rate per participant, rejection rate by stage across participants, pooled metric distributions + scatter — in `qc2b_database_report.html`, alongside `qc2b_database_rejection_summary.tsv` (one row per `file_id × stage`: `n_total`, `n_rejected`, per-method counts, `pct_rejected`). `--no-1f` skips the slow 1/f fitting. The batch is report-only (no `.fif`/override written). Outputs default to `<derivatives_root>/qc2b_reports/`.

### 8. Live single-file explorer (`8_live_explore_1file.ipynb`, `8_live_explore_1file_voila.ipynb`)

Interactive inspection of **one EDF file at a time** — load it once, then explore it live (inspired by ScoringHero). Unlike the batch tools it preloads the signal and stays interactive. It is a **QC + scoring-review companion**: it never modifies the EDF and writes outputs only on explicit button presses. Reuses the quality plots of `quality_overview` and the rejection logic of `6_preprocessing_voila`, applied to a single recording. Both delivery forms are kept in sync; the Voila version hides code (Voila strips sources by default).

**Section 1 — Load**: free file pickers for the EDF (required), the hypnogram `.txt` (optional), and a `remap_reref_persubject.json` (optional, matched to the participant by EDF filename stem, normcase-insensitive). The EDF header is parsed with the **custom binary parser** (no signal loaded) to list channels and auto-detect EEG/EOG/EMG; each channel set is offered as **checkboxes** with the detected channels pre-checked (correct as needed). EEG analysis channels come from the matched JSON remap keys when available, otherwise from the EEG checkboxes. Loading uses the established `include=`-at-read-time + `drop_suffix_duplicates` + `adapt_remap_dict_to_suffixes` pattern; EOG/EMG keep their original names. **Sampling rate**: the EEG native rate is the reference — the assembled montage is resampled to it so the whole scoring montage shares one time base. A global **Signal** toggle (Raw ↔ Preprocessed) governs every section; the preprocessed copy applies the per-participant re-reference (JSON `ref_channels`, or a manual average/none when no JSON) plus an optional bandpass (default 0.1–40 Hz). Reference channels are not dropped in the preprocessed copy so the channel set stays identical to the raw copy.

**Section 2 — Whole-recording overview (EEG)**: per-electrode quality view (amplitude histogram with Savitzky-Golay smoothing + peak detection, full time series, metrics table with the same flags/thresholds as `quality_overview`, YASA hypnospectrogram, whole-night PSD, and **mean PSD per sleep stage**). A **Run** button precomputes the figures for **all** electrodes for **both** the raw and preprocessed signal and caches them; switching the electrode dropdown or the Raw/Preprocessed toggle then displays the cached plots instantly, so the impact of preprocessing is immediately visible without recomputation.

**Section 3 — Epoch explorer (EEG + EOG + EMG)**: gated behind a **Run** button. The **navigator** shows a spectrogram of a reference channel (viridis, colour-scaled to the p5–p99 dB range so it isn't flattened by the silent floor / artefact spikes) with the hypnogram line overlaid (twin axes: frequency left, stage right; white line with a dark halo for contrast), a bold cursor at the current epoch and enlarged event / modified-epoch markers. Run precomputes the navigator spectrogram for **both** signal sources so the Raw/Preprocessed toggle switches instantly (the epoch view is recomputed live); changing the reference channel needs a new Run. The **epoch view** stacks the scoring montage (EEG, then EOG, then EMG) for the current 30 s epoch (± context), with `*_event_xml.csv` annotations drawn as shaded spans, plus the per-epoch PSD and a p-p/gradient readout. **Navigation**: prev/next, jump-to-epoch slider, next-stage-change, and (when `ipyevents` is installed) **keyboard shortcuts** — ←/→ to step epochs, and `w/1/2/3/r/m/0` to score the current epoch. **Rescoring**: stage buttons reassign the current epoch (held in memory, modified epochs ticked on the navigator); **Save** writes `<hypno_stem>_rescored.txt` + `<hypno_stem>_rescore_log.tsv` next to the input hypnogram (AASM labels, one per line, original never overwritten).

**Section 4 — Quick epoch rejection (EEG, standalone)**: a **Compute** button runs the five rejection methods (per-stage amplitude, flat, gradient, 1/f error, 1/f R² — same logic and defaults as `6_preprocessing_voila`; 1/f off by default for speed and only if `specparam` is available) on the EEG channels of the active signal source. Shows the channels × epochs rejection heatmap with a hypnogram strip, a per-stage rejection summary, the **mean PSD per stage over clean vs rejected epochs** (gold-standard sanity check), and a **rejected-epoch inspector** (a dropdown of flagged epochs renders that epoch's montage + PSD inline and jumps the Section 3 explorer to it). **Save** exports `<edf_stem>_live_rejection_mask.tsv` (same per-(epoch, channel) schema as `6_preprocessing_voila`'s mask) next to the EDF — standalone QC, no `.fif` is written.

**Rendering**: static matplotlib PNG for v1 (fast, no extra dependency, embeds into reports). A future migration of Section 3 to an interactive canvas (`ipympl`/Plotly) for click-to-seek and direct mouse annotation is tracked in `tools/TODO_live_explore_interactive_migration.md`.

**Event annotations**: scored events are loaded with the shared **CSV-first / XML-fallback** `load_events(edf_path)` — the Compumedics `*_event_xml.csv` (`Name, Start, Duration` in seconds) is read first and, when absent, the `<ScoredEvents>` of the `*.edf.XML` (Profusion `CMPStudyConfig`) is parsed as a fallback (returns a `Name/Start/Duration` DataFrame plus a `source` tag, surfaced in the load summary as `events (from csv|xml)`). Its events (arousals, apnea/hypopnea, limb movement, SpO2 desaturation…) are overlaid colour-coded on the Section 3 epoch view and ticked on the navigator.

**Custom (non-AASM) stages** (see *Cross-cutting procedures*): an editable `Custom stages` field (Section 1, auto-filled from the EDF folder's `config_param/custom_stages.json`) sets a global `CUSTOM_STAGES` at load. The Section 2 hypnospectrogram uses the shared `plot_hypnospectrogram()` (custom labels stacked below N3); `validate_hypno` treats registered labels as known (no warning); and the per-stage mean PSD, the Section 4 clean-vs-rejected PSD, the rejection heatmap strip, the per-stage summary and the navigator hypnogram axis all iterate `STAGES + CUSTOM_STAGES`. Section 4's per-stage amplitude rejection uses the existing `ptp_thresholds.get(stage, 250.0)` fallback for custom stages (quick single-file preview — the authoritative per-stage thresholds live in tool 6).

**Outputs** (written only on explicit Save): `<hypno_stem>_rescored.txt` + `<hypno_stem>_rescore_log.tsv` (next to the input hypnogram), `<edf_stem>_live_rejection_mask.tsv` (next to the EDF).

### 9. Spectral Analysis (`9_SpectralPower_&_AperiodicFit_PSG.py`)

Full PSG spectral pipeline: epoch rejection → PSD (Welch, 4 s windows) → aperiodic fit (SpecParam) → frequency band power extraction (Delta, Theta, Alpha, Sigma, Beta) → group-level statistics. Reads the channel remapping JSON produced by tool #2.

**Planned**: adapt this batch script into a Voila/Jupyter notebook (keeping a `.py` batch twin) so it integrates with the rest of the toolbox like the other tools.

## Planned modules (in development)

Modules still in development. The quality-overview (5), preprocessing (6) and QC-of-rejected-epochs (7) tools above are stable and were promoted out of this section.

### Event-based epoch rejection (Phase A) — **implemented in tool 6**

Promoted out of this section: the **"Event overlap"** rejection method is now part of `6_preprocessing_voila.ipynb` (Section 2 → *Epoch rejection*). See the tool-6 description above for the UI, the `event_remap.json` chooser + auto-detected `Event CSV suffix`, the "Count affected epochs" button, and the additive `flag_event` outputs. (The mirror in `8_live_explore_1file`'s quick rejection remains a possible follow-up.)

### Event-epoch visualizer (Phase B)

Helps decide whether a given event *type* is worth feeding into Phase A. Extends `8_live_explore_1file`'s Section 3 (which already overlays scored-event spans via the shared CSV-first / XML-fallback `load_events()`):
- an **event-type filter** on the navigator: "jump to next/previous epoch overlapping event type X", with a per-type count;
- a per-type **accept/reject decision** widget that writes a small `event_type_decisions.tsv` feeding Phase A's default selection.

### Interactive QC of rejected epochs (Phase 2b) — **implemented as tool 7**

Promoted out of this section: implemented as `7_inspect_rejected_epochs_voila.ipynb` (+ `7_inspect_rejected_epochs_batch.py`, `qc_rejected_epochs_lib.py`) — see the **tool-7 section above**. It loads `{file_id}_all-epo.fif` (+ optional `{file_id}_preprocessing_params.json`), lets the user navigate flagged epochs and override the whole-epoch keep/reject decision, and saves a validated `{file_id}_clean-epo.fif`. **Possible follow-ups**: the code-visible Jupyter twin; optional EOG/EMG context in the per-epoch view (needs the context-channel declaration in tool 2 — see `tools/plan_tool2_eog_emg_context_channels.md`); per-(epoch, channel) override if a future tool-6 run persists the channel-level mask.

### Test data infrastructure

`tools/generate_test_data.py` produces controlled-defect EDF files from a clean baseline (`tools/test_data/73.edf`). It uses the channel list from `tools/test_data/config_param/remap_reref_persubject.json` to keep only the relevant EEG channels (currently `A2, C3, Fp1, O1`), making each output ~110 MB. Each generated file injects **one** defect on a specified channel; a `combined` file mixes three defects on different channels for integration testing.

Each defect file is given a unique participant-like ID (731–738) so that `2_select&remap_channels_edf` does not group them all under participant `73`. The ID mapping is defined in `DEFECT_IDS` in `generate_test_data.py`. The script also writes entries for each generated file into `config_param/remap_reref_persubject.json` (copied from the `73` entry).

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
