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
│   ├── 7_reject_manually_voila.ipynb           # Manually reject flagged epochs — Phase 2b (Voila GUI)
│   ├── 7_reject_manually_batch.py              # Manually reject flagged epochs — database-level batch report
│   ├── qc_rejected_epochs_lib.py                # Shared analysis/plotting for tool 7 (notebook + batch)
│   ├── 7bis_reject_automatically_voila.ipynb   # Automatic epoch rejection (channel-first then epoch, Voila GUI)
│   ├── 8_live_explore_1file.ipynb               # Interactive single-file explorer (Jupyter)
│   ├── 8_live_explore_1file_voila.ipynb         # Interactive single-file explorer (Voila GUI)
│   ├── 9_SpectralPower_&_AperiodicFit_PSG.py    # Spectral analysis pipeline
│   ├── generate_test_data.py                  # Inject controlled defects into a clean EDF (test fixtures)
│   ├── test_data/                             # Real EDF fixtures + generated defective files + manifest
│   ├── images/                                # Reference images for quality checks
│   ├── preprocessing_phase1_example_scripts/  # Draft/example scripts used during Phase 1 development
│   └── old/                                   # Versioned development notebooks (archive)
└── tools_curry/                             # Experimental Curry 9 (.cdt) port — see "Curry 9 support" below
    ├── curry_header.py                        # Header-only .cdt.dpo parser (Curry analogue of the EDF header parser)
    ├── curry_io.py                            # Curry signal / hypnogram / events (text export) loaders
    ├── {1,2,3,4,5,6}_*_curry_voila.ipynb      # Curry twins of tools 1–6 (Voila)
    ├── 7bis_reject_automatically_curry_voila.ipynb  # Curry twin of tool 7bis (verbatim copy — 7bis is format-agnostic)
    └── _make_tool{2,3,4,5,6,7bis}_curry.py    # Re-runnable generators (regenerate a twin from its EDF source)
```

**Sibling directory** `../Check_EDF/` contains exploratory notebooks used during development (not production tools).

### Editing the larger notebooks (agent procedure)

`Read` ignores `offset`/`limit` on `.ipynb` and fails once total size passes ~25 k tokens; `Edit` is
blocked on the `.ipynb` extension and `NotebookEdit` needs a prior whole-file `Read`. So for the large
notebooks, rename around the extension blocks (**rename-to-`.txt`** method):

1. `mv "<nb>.ipynb" "<nb>.ipynb.txt"` (Bash tool).
2. Edit the **raw notebook JSON** with Grep / Read (`offset`/`limit`) / Edit — each source line is a
   `"…\n",` array element; preserve escaping (`\"`, `\\`, `\n`) and array commas (the last element of a
   `source` array has no trailing comma).
3. Validate:
   `& "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" -c "import json; json.load(open(r'<nb>.ipynb.txt', encoding='utf-8'))"`.
4. `mv "<nb>.ipynb.txt" "<nb>.ipynb"` — keep the whole rename in one task so the git diff stays byte-exact.

- **Small notebooks** (`1bis_anonymize_edf*`, `4_remap_events_edf*`): normal `Read` + `NotebookEdit`, no
  rename needed.
- **JSON-surgery fallback** (many repetitive / escaping-heavy edits): a one-shot Python script that
  `json.load`s, string-replaces inside the parsed cell `source` (assert each pattern matches exactly once),
  and `json.dump`s back with `indent=1, ensure_ascii=False` + trailing newline. **Preserve the original
  type of `cell["source"]`** — a single string must stay a single string (assigning a list back
  re-serializes one physical line into ~1500, exploding the diff); if you must write a list, split with
  `splitlines(keepends=True)` (never `split("\n")`) and re-`compile()` the joined source as a syntax gate.
  An open IDE may re-serialize the notebook between calls — re-read before each pass. Prefer rename-to-`.txt`
  when diff minimality matters.

## Conda environment

Defined in `environment.yml`. Key packages:
- **Python 3.12.10**
- **MNE 1.12** — EDF reading, epoching, signal processing (bumped from 1.9 for Curry `.cdt` support; the EDF tools run unchanged on it)
- **YASA 0.6** — sleep staging, hypnogram handling, spectral helpers
- **pandas 2.2, numpy 2.2** — data manipulation
- **voila 0.5, ipywidgets 8.1, ipyfilechooser 0.6** — interactive GUI layer
- **chardet 5.2** — encoding detection for EDF headers
- **edfio** — EDF read/write (used directly for export with per-channel physical range control)
- **specparam** — aperiodic/periodic spectral decomposition (1/f fitting)

## How to run the tools

**Canonical working directory = the repo root (`Inspect_EDF/`).** Launch every tool from there (paths below
are relative to it). This keeps one consistent cwd for the shared libraries, `config_param/*.json`, hypnogram
lookups and outputs across the whole suite. As a safety net the tools also resolve their shared modules
independent of cwd — the EDF tools' `tools/qc_rejected_epochs_lib.py` and the Curry tools'
`tools_curry/curry_header.py` / `curry_io.py` are found whether Voila is launched from the repo root **or**
from the tool's own folder — but the repo root is the documented, recommended launch directory.

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
voila tools/7_reject_manually_voila.ipynb
voila tools/7bis_reject_automatically_voila.ipynb
voila tools/8_live_explore_1file_voila.ipynb
```

### Curry 9 (`.cdt`) Voila twins (run from the repo root, same as the EDF tools)
```bash
conda activate inspect_edf
voila tools_curry/1_inspect_curry_voila.ipynb
voila "tools_curry/2_select&remap_channels_curry_voila.ipynb"
voila tools_curry/3_remap_hypno_curry_voila.ipynb
voila tools_curry/4_remap_events_curry_voila.ipynb
voila tools_curry/5_quality_overview_curry_voila.ipynb
voila tools_curry/6_preprocessing_curry_voila.ipynb
voila tools_curry/7bis_reject_automatically_curry_voila.ipynb
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
  - **Interruption-safety (the skip gate requires data, not just a report)**: an item is skipped only when
    **both** its human report **and** its durable per-item data are on disk; if exactly one is present
    (run interrupted between the two, or a manual delete) the item is **reprocessed** and a `⚠` mismatch
    warning is surfaced (tool 5: in the loop output + the folder-selection info line; tool 6: in the
    Section 3 participant-loading `load_info`). Per-item data is written **before** the report so
    "report ⇒ data" always holds going forward, and every cumulative/aggregated table is rebuilt by
    **globbing the per-item files on disk** — never from an in-memory list keyed on the current run's
    `attempted_ids`. This closes the failure where an interrupted-then-skipped item's rows were lost
    forever. Tool 5 markers: `{file_id}_quality_metrics.tsv`; tool 6 skip marker: `{file_id}_epoch_channel_rejection.tsv`
    (the per-(epoch, channel) TSV — gated so files processed by an older version, which lack it, are
    reprocessed and back-filled; its `global_epoch_rejection.tsv` is globbed from disk like
    `global_rejection_by_stage.tsv`, excluding the global file itself since it shares the
    `_epoch_rejection.tsv` suffix, and the per-file `_epoch_channel_rejection.tsv` is *not* concatenated
    globally — per-file only).
- **Lenient JSON loader**: every read of a config JSON parses strictly first and, on failure, repairs a
  single trailing comma before a closing `}`/`]` (a common hand-edit mistake) before retrying. Shared by
  the config readers of tools 2 and 4.
- **Hypnogram-suffix auto-detection**: on folder selection, the `.txt` files next to each EDF are scanned,
  candidate suffixes counted (all files per EDF matched), and the suffix widget auto-filled, with a
  colour-coded info label (green = all EDFs matched, orange = partial/none). **Selection rule differs by
  intent**: tool 3 prefers the **shortest** suffix on ties (target = raw, unremapped hypnogram); tools 5
  and 6 prefer the **longest** suffix among candidates appearing for ≥ 50% of the maximum count
  (target = the more specific remapped/processed version).
  **Event `.txt` exclusion + display/selection split** (tools 3, 5, 6): the auto-*selection* skips
  candidate suffixes containing `event` (case-insensitive, `sel_counts`), while the *displayed* candidate
  list keeps **every** `.txt` suffix found, with `← selected` marking the one auto-filled into the widget.
  Curry exports scored events as `*_ScoredEvents_Export.txt` — a `.txt` living next to the hypnograms whose
  suffix (24 chars) is **longer** than `_Hypnogram_remapped.txt` (23) — so the "prefer the longest" rule of
  tools 5/6 would otherwise auto-select the event export and read event lines as sleep stages (`np.loadtxt`
  then fails with a *column count changed* error whose numbers vary per file, since it splits event labels
  on whitespace). Excluding `event` rather than requiring `hypno` leaves hypnogram naming unconstrained;
  displaying all candidates keeps a mis-detection visible and hand-correctable. The ≥ 50% threshold is
  computed **within** the selectable (non-event) set, so a numerous event export cannot raise the bar high
  enough to disqualify a partially-remapped hypnogram; when *only* event suffixes exist the selection falls
  back to them rather than failing. Same `event` convention, opposite polarity, in tool 4 (`"event" in name`
  — it *wants* the export). **Live in the EDF `4_remap_events`** (auto-detects the `*_ScoredEvents_Export.txt`
  export, see below) as well as the Curry twins; still inert in tool 8 (its `load_events()` remains
  `.csv`/`.XML` only — extending it to the `.txt` is a possible follow-up). Tool 3 was already immune via
  "shortest wins", the exclusion only makes it explicit.
- **Event sourcing (TXT-first / CSV / XML-fallback)**: scored events are read via a shared `load_events()`.
  **Tool 4 (`4_remap_events_edf*`)** reads three Compumedics companions in priority order: the
  `*_ScoredEvents_Export.txt` **text export** (the "classic" Profusion/Curry French export — UTF-16-or-UTF-8
  comma-separated, **no header**, clock-time + `M:SS[.s]` duration; BOM-sniffed and parsed inline by
  `load_events_from_txt`, mirroring `curry_io._parse_events_txt`), then the event CSV
  (`Name, Start, Duration`, default suffix `_event_xml.csv`), then the `<ScoredEvents>` of the `*.edf.XML`
  (`CMPStudyConfig`). Two editable, auto-detected suffix fields (**TXT suffix** + **CSV suffix**);
  `load_events(edf, txt_suffix, csv_suffix)` returns `(list-of-(name, start, duration), source∈{txt,csv,xml})`.
  For the `.txt` the harmonization scan needs **only the names**, so Start/Duration are left `NaN` (no EDF
  header read); the recording-start datetime (`read_edf_start_datetime`, EDF header offsets 168/176) is read
  **only** by the 1bis consistency check to convert the text-export clock times to seconds. The **Curry twin**
  keeps the single `.txt` source via `curry_io.load_events_curry` (no 1bis).
  **Tool 6 (`6_preprocessing*`)** uses the **same TXT-first / CSV / XML-fallback** `load_events()`, but
  returning a `Name/Start/Duration` **DataFrame** (its `compute_event_epoch_mask` consumes a DataFrame).
  Crucially — unlike tool 4's harmonization scan, which needs only the event *names* — tool 6 flags an
  epoch by the event **onset (`Start`)**, so its `.txt` branch **reads the recording-start datetime
  (`read_edf_start_datetime`, EDF header offsets 168/176) and converts the export's clock times to real
  seconds** (`_events_df_from_txt`). A `.txt` present without a readable start datetime is skipped in favour
  of the CSV/XML companions. Two editable, auto-detected suffix fields (**Event TXT suffix** +
  **Event CSV suffix**); `load_events(edf, txt_suffix, csv_suffix)` returns `(DataFrame, source∈{txt,csv,xml})`.
  The **Curry twin** mirrors the same chain (`_make_tool6_curry.py` keeps `_events_df_from_txt/csv/xml` +
  the dispatcher unchanged and swaps only `read_edf_start_datetime` for a `.cdt`-header read via
  `curry_header`/`rec_start_from_header`), so a Curry dataset shipping a `_event_xml.csv` is picked up too.
  **Tool 8** is unchanged: CSV-first / XML-fallback, returning a `Name/Start/Duration` **DataFrame** plus
  the `source` tag (its overlay/navigator code consumes a DataFrame).
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
    (byte-identical scale and image). **Dynamic-range cap (`vmin = max(vmin, vmax − 45)`)**: the
    60 dB column gate only removes *fully* dead epochs; a recording with a **continuum of partly-flat /
    clipped low-power epochs** (e.g. `1WIBE0543_N1`: ~20 % flat + ±500 µV clipping on every EEG channel)
    leaves degraded-but-valid columns whose peak dB sits ~20–60 dB below the median. Those survive the
    gate but still drag the pixel-percentile `vmin` to ≈ −59 dB, pushing RdBu_r's white midpoint down to
    ≈ −22 dB so all real EEG (> −20 dB) washes to red again. Capping the colour **span** at 45 dB after
    the percentile fixes this: a healthy hypnospectrogram's real structure fits within ~35 dB (measured
    max across all clean `test_data` channels = 34.6 dB), so the cap **never touches** clean channels
    (their span < 45 → `max()` is a no-op → still byte-identical) yet restores contrast on degraded ones;
    the degraded epochs stay **visible in blue** (low power), only the dead ones are greyed. The quantitative
    QC flags (`flat_pct`/`bounds_pct`/…/`exclude`) are independent of this colour scale, so detection is
    unaffected. **Absolute-dB colorbar** (`render_hypnospectrogram`, tool 5 per-channel + median):
    because the colour scale is *relative* (auto-scaled per channel), a **globally attenuated / low-gain**
    channel looks structurally normal — the within-channel greying/blue logic can't reveal a uniform
    problem. A `Power [dB]` colorbar (attached to **both** the hypnogram strip and spectrogram axes so
    they stay x-aligned) exposes it: e.g. the ×0.01 dead channel `734 C3` sits at −55…−23 dB vs a healthy
    channel's −15…+18 dB, an ~40 dB downward shift visible at a glance. This is a display-only change
    (it re-lays-out every per-channel figure, clean ones included, on purpose); the numeric `std_uV` /
    "All electrodes" figures remain the primary catch for globally-bad channels. Present in **tools 5 & 8**
    (tool 8's monolithic `plot_hypnospectrogram` carries the same colorbar block — keep in sync). Applied to the **across-channel median** spectrogram of tool 5's
    Overview (see §5), the very same criterion yields *majority-of-channels* semantics for free: a
    column is greyed only when most channels are dead at that epoch, since that is what makes the
    median drop — no separate majority-vote code. The tool-7 *navigator* spectrogram is a separate plot
    (floored at −120 dB, p5–p99) and is left as-is. Diagnostic scripts:
    `tools/simple_hypnospectro_yasa_vs_fix.py` and `tools/compare_flat_spectrogram_fix.{py,ipynb}`.
- **Time-series display cap for DC-coupled data (±500 µV physiological ceiling)**: DC-coupled recordings
  (Curry `.cdt`, and any acquisition exported in DC with no clipping) carry no export clipping, so slow
  drift or artefacts can push the p99.9-based autoscale far past physiological range and crush the real
  EEG in a time-series / butterfly plot. For such tools the shared amplitude limit is **capped** at a wide
  physiological ceiling — `y_lim = min(max_p999, 500.0)` (constant `DISPLAY_YLIM_UV = 500.0`) — never a
  hard fixed window, so clean low-amplitude channels still auto-zoom below the cap and the shared
  cross-channel scale is kept. Applied to the per-channel + butterfly time series **and** the histogram
  X-axis (`x_lim_hist` follows `y_lim_ts`). Currently **Curry-only** (injected by
  `tools_curry/_make_tool5_curry.py`, block "cap time-series y-limit…" — re-run the generator after
  editing); the EDF tools keep the uncapped autoscale on purpose (full range helps spot export clipping).
  Extend the same cap to any future DC-source tool.
- **Physical bounds in µV (`get_phys_bounds_uV`)**: MNE stores an EDF channel's physical range as
  `physical_max + offset` in `raw._raw_extras` (not explicit `physical_min`/`physical_max`); the
  `units`/`physical_max`/`offsets` keys are present and identical on MNE 1.9 and 1.12.
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
- **Context-channel detection (`detect_context_channels`, tool 2)**: a **header-scan** variant of the same
  *transducer-type OR name* convention, used by `2_select&remap_channels_edf*` Section 2bis to declare
  per-configuration EOG-Left/EOG-Right/EMG/ECG "context" channels. Unlike `detect_channel_types` (which
  reads an MNE `Raw`), this one operates on the scanned header dataframe (`df_full`) because tool 2 never
  loads signal, and it **tolerates a missing `transducer_type` column** so the generated Curry twin still
  works by channel name alone. EOG-side heuristics: `E1|LOC|left|gauche` → left, `E2|ROC|right|droit` →
  right (any unsided EOG fills the first empty L/R slot; user can swap).
- **ipywidgets `Box` stray per-row scrollbars**: jupyter-widgets ships
  `.widget-box { box-sizing: border-box; overflow: auto; }`, so any bordered + padded `HBox`/`VBox` row
  whose children overflow the content box by even 1–2 px renders a per-row ▲▼ vertical scrollbar the user
  can scroll by accident (shifting the row content). Fix: set `overflow='hidden'` explicitly on such row
  layouts (applied to the Section 3 "Harmonize labels" rows of `4_remap_events_edf*`); keep the intended
  scroll only on the outer list container.
- **Checkbox-revealed parameter widgets (initial `display` follows the checkbox)**: the optional-parameter
  boxes (resampling target, high-pass corner, notch frequency) are hidden until their checkbox is ticked,
  via an `observe(..., names='value')` handler setting `w.layout.display = '' if change['new'] else 'none'`.
  An observer only fires on a **change**, so a box whose layout hard-codes `display='none'` stays invisible
  under a checkbox that *starts* ticked — the user has to untick and retick to reveal it. The initial
  visibility must therefore be **derived from the checkbox value** at construction:
  `layout=widgets.Layout(width='320px', display='' if cb_resample.value else 'none')` (the checkbox is
  defined just above, so its `.value` is available). This is invisible while every such checkbox defaults
  to `False`, but the **Curry twin of tool 5 ticks `hp_check` ON by default** (DC-coupled data, injected by
  `_make_tool5_curry.py`) — which is where the bug surfaced. Applied to all four boxes: `txt_target_freq`
  + `hp_freq` in `5_quality_overview_voila`, `txt_target_freq` + `txt_notch_freq` in `6_preprocessing_voila`
  (tool 6's bandpass `txt_l_freq`/`txt_h_freq` are always visible — no toggle — and are unaffected). Edit
  the **EDF** notebooks and re-run the Curry generators; neither generator string-matches these blocks, so
  the pattern passes through. Apply the same derivation to any new checkbox-revealed widget.

## Tool descriptions

### 1. EDF Inspector (`1_inspect_edf_voila.ipynb`, `1_inspect_edf_perdataset.py`, `1_inspect_edf_perparticipant.py`)

Inspects EDF file parameters across an entire dataset **without loading signal data**. Reads EDF headers using a custom binary parser to handle encoding edge cases robustly (see design decisions in CLAUDE.md).

**Sampling frequency derivation**: the per-channel 8-byte header field is the *number of samples per data record*, **not** the sampling frequency. The parser stores it as `samples_per_record` and computes `sampling_frequency = samples_per_record / duration_data_record`. In classic EDF the data-record duration is 1 s, so the two values coincide; but EDF+ files frequently use a different record duration (e.g. `0.1 s` → `40 / 0.1 = 400 Hz`, or `2 s` → `512 / 2 = 256 Hz`), so dividing by `duration_data_record` is required to match the rate reported by MNE. The computed value is kept as a string (e.g. `'256'`, `'400'`) to preserve the existing `sorted(set(...))` grouping and avoid TSV round-trip type changes. This applies to every tool sharing the custom header parser: `1_inspect_edf_voila.ipynb`, `1_inspect_edf.ipynb`, `1_inspect_edf_perdataset.py`, `1_inspect_edf_perparticipant.py`, and `2_select&remap_channels_edf(_voila).ipynb`.

**Checks performed for EEG, EOG, ECG, and EMG channels:**
- Channel configuration and montage consistency across participants
- Sampling frequency consistency
- Filter settings consistency
- Signal units
- Inverted polarity (physical_min > physical_max)
- Signal clipping (dynamic range ≤ 500 µV)
- Poor resolution (dynamic range ≥ 0.1 µV per digital unit)

**Channel-type selection + robust detection**: which types are inspected is user-selectable, **default
EEG + EOG** (ECG and EMG are opt-in). Detection is harmonized across all four files to the same
*transducer-type OR curated channel-name list* convention as the shared `detect_channel_types` (see
CLAUDE.md / *Cross-cutting procedures*): EEG = transducer `EEG`/`AGAGCL ELECTRODE` OR `KNOWN_EEG_CHANNEL_RE`,
then **subtract** any `emg|ecg|eog|ekg|chin|menton` channel name (drops non-EEG sensors captured via the
generic AGAGCL-ELECTRODE transducer); EOG adds `LOC|ROC|E1|E2` names; ECG adds `EKG`; EMG = transducer
`EMG` OR `emg|chin|menton` (the `chin|menton` aliases cover Compumedics chin-EMG labels). The selection
surface differs per delivery form:
- **Voila** (`1_inspect_edf_voila.ipynb`): four checkboxes (`EEG`/`EOG`/`ECG`/`EMG`, EEG+EOG ticked) shown
  right after the folder chooser; `run_inspection` gates each per-type section on them (section 5 = EMG),
  so a run only produces the selected types' tables/report blocks. General dataset info and the
  anonymization check are type-independent and always run.
- **Jupyter** (`1_inspect_edf.ipynb`): **no checkboxes** — each type is its own `## N. Inspect X` section
  (section 5 = EMG) that the user chooses to run cell-by-cell, matching the notebook's manual-run model.
- **Batch scripts**: a top-of-file `INCLUDE_TYPES` dict selects types. `perparticipant.py` supports all
  four (EMG block mirrors ECG). `perdataset.py` intentionally aggregates **EEG + EOG only** — ECG/EMG were
  deliberately left out of the dataset-level report to keep it light; use `perparticipant.py` or the
  notebooks for per-participant ECG/EMG.

Selecting only the default types (or leaving EMG off) keeps outputs byte-compatible with the pre-selection
tool. The former standalone `inspect_edf_voila_EMG.ipynb` (a pre-skip/merge fork that hard-coded EMG) is
superseded by this and can be archived.

**Anonymization check** (section 1.3 in Voila / section 1.4 in Jupyter): inspects the EDF+ *Local Patient ID* field (80-byte header) to detect non-anonymized patient names. The EDF+ format encodes this field as `code sex birthdate name` (space-separated); Compumedics writes the name as `LASTNAME_FIRSTNAME` and replaces it with `X_X` on anonymized export. The check isolates the name sub-field (4th token onward), strips placeholder characters (`X`, `x`, `_`, `,`, `;`, whitespace), and flags the file if anything remains. The `,` and `;` separators are stripped so that headers anonymized by *other* systems using a `Lastname,Firstname` placeholder (e.g. `Xxxxxxx,Xxxx`) are still recognized as anonymized, not just the Compumedics `X_X` form. Additionally, each non-placeholder name token (≥ 3 characters, to avoid false positives from short codes or initials) is searched case-insensitively in the file stem to detect PII leaking into the file name. Two warning levels:
- `PII in header AND file name` — real name found in both header and file name.
- `header NOT anonymized (file name looks clean)` — real name in header but file name appears clean; the important case where the file was renamed but the header was forgotten.
Files where the check cannot be performed (read failures) are already captured in `failed_edf_read.tsv`. The `patient_name` field (raw name sub-field, before cleaning) is also stored as a column in `FULL_summary_table_edf.tsv` for quick cross-reference.

**Outputs** (all written to `<study_folder>/summary_inspection/`, created on first run with a `README.md`):
- `FULL_summary_table_edf.tsv` — all parameters for all channels/files, including `patient_name` column (raw name sub-field from the EDF+ `patient_id` field)
- `anonymization_check_edf.tsv` — per-file anonymization status; columns: `subject`, `path`, `patient_id`, `name_subfield`, `patient_id_format`, `header_anonymized`, `name_in_filename`, `anon_warning`
- `EEG_summary_table.tsv`, `EOG_summary_table.tsv`, `ECG_summary_table.tsv`, `EMG_summary_table.tsv` (only the selected types are written; ECG/EMG only when enabled)
- `EEG_inverted_polarity_edf.tsv`, `EEG_bad_dynamic_range_edf.tsv`, `EEG_bad_resolution_edf.tsv` (and the matching `EOG_*`/`ECG_*`/`EMG_*` set per selected type, plus `<TYPE>_missing_edf.tsv`/`<TYPE>_suspect_edf.tsv`)
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
2. Review table — one row per file to process: an include checkbox (pre-ticked for non-anonymized files), an editable new-name field, and a colour-coded badge (`name in header AND file name` / `name in header only` / `already anonymized`). Options: anonymize `recording_id` trailing fields (on), skip files already in the output folder / recompute everything (on), also list already-anonymized files (off). **Interruption-safe skip** (see *Skip + cumulative-merge*): a file is skipped only when **both** its anonymized copy **and** its `anonymization_log.tsv` row exist; a copy present without a log row (run interrupted before the end-of-loop log write) is **re-anonymized** with a `⚠` warning to restore the audit trail (copy + header patch is idempotent), rather than skipped and lost forever.
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

**Section 2bis — Identify context channels (EOG L/R, EMG, ECG)**: an **optional** step that declares, **per
channel configuration**, the non-EEG *context* channels — `eog_left`, `eog_right`, `emg`, `ecg` — so any
downstream tool can pull them in on demand without re-detecting them ad hoc. A "Run context channels" button
builds one accordion panel per configuration, each with four editable `Dropdown`s auto-filled by
`detect_context_channels` (see *Cross-cutting procedures*; pick `(none)` if a channel is absent). "Save
context selection" stores the choices in `context_by_config`, keyed by config label, which Section 5 fans out
to each participant like `remap`/`ref_channels`. The result is an **additive, backward-compatible**
`context_channels` block appended to each participant entry, its **keys the original EDF channel names** (so a
downstream tool can pass them straight to `include=` at read time):
```json
"73": { "config": "config. 1", "remap": {…}, "ref_channels": ["M2"],
        "context_channels": {"eog_left": "E1", "eog_right": "E2", "emg": "Menton", "ecg": "EKG"} }
```
The block is **omitted entirely when no context channel is selected**, so `remap_reref_persubject.json` files
produced without using Section 2bis stay byte-identical. Absent/`null` entries mean "not declared" and every
consumer must treat that as non-fatal. Analysis tools 5/6 ignore this block (they read only the `remap`
keys); only tools that explicitly ask for it (e.g. tool 7's per-epoch inspector) load these channels.

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

**Event sourcing (TXT-first / CSV / XML-fallback)**: a shared `load_events(edf_path, txt_suffix='_ScoredEvents_Export.txt', csv_suffix='_event_xml.csv')` helper reads three Compumedics companions in priority order — (1) the `*_ScoredEvents_Export.txt` **text export** (the "classic" Profusion/Curry French export: comma-separated, **no header**, `HH:MM:SS` clock time + `M:SS[.s]` duration + French label; **encoding sniffed** — UTF-16 with BOM or UTF-8/ANSI — and parsed inline by `load_events_from_txt`, mirroring `curry_io._parse_events_txt`), (2) the event CSV (`Name, Start, Duration` in seconds), (3) the `<ScoredEvents>` of the `*.edf.XML` (Profusion `CMPStudyConfig`; `<Input>` ignored). For the harmonization scan **only the names are needed**, so the `.txt` clock times are *not* converted (Start/Duration `NaN`, no EDF-header read); `read_edf_start_datetime()` (EDF header offsets 168/176, 2-digit-year clipping) supplies the recording-start datetime that converts them to seconds **only** in the 1bis check. CSV and XML were verified equivalent on ICEBERG. The `<ScoredEventSettings>` catalogue is **not** used for grouping.

**Configurable event suffixes**: two editable text fields — **`TXT suffix:`** (default `_ScoredEvents_Export.txt`) and **`CSV suffix:`** (default `_event_xml.csv`) — drive `event_companion_paths(edf_path, txt_suffix, csv_suffix)`, so datasets exported with different suffixes are supported. On folder selection **both** are **auto-detected** (mirroring the hypnogram-suffix detection of `5_quality_overview`, see *Cross-cutting procedures*): for the TXT the `.txt` files **whose name contains `event`** are scanned (so hypnogram `.txt` are not counted), for the CSV the `.csv` files; candidate suffixes are counted and each field auto-filled with the **most frequent** suffix (shortest on ties), each with its own colour-coded info line (green = all EDFs matched, orange = partial/none). The XML fallback (`.edf.XML`) is auto-derived, no field.

**Configuration grouping**: files are grouped by their `frozenset` of **unique event labels actually present** — two files with the same unique labels share one configuration even if their event counts/timing differ.

**Workflow (sections):**
1. **Scan** — select the data folder (recursive `rglob('*.edf')`); selecting the folder only refreshes an info line **and auto-detects both the TXT and CSV suffixes** (editable `TXT suffix:` / `CSV suffix:` fields, colour-coded detection lines — see *Configurable event suffixes* above), the scan runs on an explicit **Run scan** button. The scan summary reports the per-source counts (`N from TXT, N from CSV, N from XML fallback`). A **"Skip labels already mapped"** checkbox (on by default) hides labels already present in an existing `event_remap.json` (incremental harmonization when a new cohort is added).
1bis. **(Optional) text/CSV vs XML consistency check** — opt-in button; for every file having both a **primary text source** (the `.txt` if present, else the CSV) **and** the `*.edf.XML`, compares the two `<ScoredEvents>` descriptions **language-robustly** and writes `event_source_mismatch.tsv`. Both sources are first passed through `_canon_events()` — names normalized to the **canonical vocabulary** (`suggest_canonical`, raw-lowercase fallback) so English (XML) and French (`.txt`) labels compare equal, and start times floored to **whole seconds** (the `.txt` truncates clock times to the second; converted via `read_edf_start_datetime`). `_match_events()` then greedily pairs events **within ±`tol` seconds** (editable **`Match tol (s)`** `BoundedIntText`, default **1**, `0` = strict same-second) in two passes: (1) same canonical label within `tol` (nearest wins), (2) any leftover within `tol` → a **co-occurring pair with different labels**. Each event ends up **matched** (same label, within `tol`), **cooccur_difflabel** (paired in time but different canonical label — a candidate same-event whose names are not yet harmonized, e.g. a cross-language pair: the key discovery aid; if the paired labels are unrelated they may instead be two distinct events within `tol`), **only_in_primary** (no counterpart — e.g. an export that omits snoring) or **only_in_xml**. `status = match` only when all three "diff" buckets are empty. The default ±1 s absorbs the common case where the two Compumedics exports round a start one second apart (verified on ICEBERG: 2/53 hypopnea were exactly 1 s off). Opt-in because it reads both files per EDF.
2. **Visualize configurations** — because two files sharing the same label names can still form distinct configs (a config = the exact set of labels *present*, so a missing label splits it off), the configs are not stacked: a **dropdown** ("Show config:") selects one configuration to detail (its sorted unique labels + file/label counts), and a **"Show file ids" toggle button** (replacing the old `<details>` arrow) shows/hides that config's file-id list in a scrollable box. The global table of every raw label (file count + total occurrences + suggested canonical) is kept below.
3. **Harmonize labels** — one editable row per unique raw label (combobox pre-filled from `DEFAULT_EVENT_MAPPING`, free text allowed), with an **ignore** toggle (stored as `null`); filtered by the skip checkbox. Each row is a bordered, column-aligned line for readability (the row layout forces `overflow='hidden'` so the jupyter-widgets default `.widget-box { overflow:auto }` does not raise a stray per-row scrollbar — see CLAUDE.md). A **"Validate mapping & ignores"** button summarizes the choices (N mapped / N ignored / N left empty, warning on empties) and **unlocks** the Section 4 save button (which starts disabled). Editing any row after validating (or re-running the scan / toggling the skip checkbox) re-locks Section 4 and clears the previous save preview, so the saved JSON always reflects the latest Section 3 selection.
4. **Preview & save** — enabled only after Section 3 validation; builds a flat `{raw_label: canonical_label}` mapping and **merges** it into `config_param/event_remap.json` via the lenient JSON loader (labels mapped this session replace their old value, all others kept; keys sorted). Unmapped non-ignored labels are reported and not saved.
5. **Verify** — applies the saved mapping to every configuration, reports the resulting harmonized labels, and passes when no raw label is left unmapped (ignored labels count as handled). A scope dropdown can restrict the view to configs with unmapped labels.

**`DEFAULT_EVENT_MAPPING`** (editable suggestions, snake_case canonical vocabulary): apnea subtypes kept (`apnea_obstructive` / `apnea_central` / `apnea_mixed`), `hypopnea`, `spo2_desaturation`, arousal subtypes kept (`arousal_respiratory` / `arousal_spontaneous` / `arousal_limb` / `arousal`), limb laterality collapsed (`limb_movement`, `plm`), plus `snore` and `spo2_artifact` (from the French export's `Ronflement` / `Artéfact SpO2`). `suggest_canonical()` matches the English exact dict first (tolerating a trailing `(Left)`/`(Right)` marker), then falls back to **accent-insensitive French substring rules** (`FRENCH_EVENT_RULES`, via `_strip_accents`): informative tokens (`aro spont/res/plm`, `apnee obstructive/centrale/mixte`, `hypopnee`, `desaturation`, `artefact spo2`, `ronflement`, `plm`) so the variable-numbered `Micro-éveil N ARO …` labels still pre-fill. Shared verbatim with the Curry twin (the generator does not touch this block → French mapping benefits Curry directly).

**Outputs** (under `<data_folder>/config_param/`):
- `event_remap.json` — global flat `{raw_label: canonical_label}` (`null` = ignore), merged across runs
- `event_source_mismatch.tsv` — only if the 1bis text/CSV-vs-XML check is run; columns `file_id, primary_source` (`txt`/`csv`)`, match_tol_s` (the ±tolerance used)`, n_primary, n_xml, n_matched, n_cooccur_difflabel, only_in_primary` (per-type counts, e.g. `snore×426`; its total replaces the former `n_only_in_primary` column)`, only_in_xml, cooccur_label_pairs` (e.g. `bidule fr↔foobar en×1`)`, status` (`match` / `MISMATCH` / `missing text/CSV` / `missing XML` / `cannot read EDF start datetime` / `error: …`)
- `failed_event_read.tsv` — files with no readable event companion (only if any failed)

### 5. Quality overview (`5_quality_overview_voila.ipynb`)

*(Stable tool — formerly “Phase 1” of the preprocessing pipeline.)*
Implemented as `tools/5_quality_overview_voila.ipynb`. Produces one `mne.Report` HTML per participant. For each EEG channel: signal amplitude histogram with Savitzky-Golay smooth + peak detection, time series, metrics table, and a YASA hypnospectrogram (0.1–40 Hz bandpass applied per-channel just before plotting). Flags suspect channels for priority inspection. At the end of each run, generates `dataset_overview.html` — a single-page dataset-level summary with statistics and distribution plots per electrode, consumed by Phase 2 to identify channels to exclude.

**Signal preparation (analysis only — not saved)**: a UI section (named *Signal preparation*, **not** *Preprocessing*, to make clear nothing is written to disk — unlike tool 6; the transforms are applied transiently, only to compute the overview). When the **config JSON is selected** (`update_acq_info`, registered on `fc_config`; needs both the folder and the config), an **acquisition scan** reads the EDF headers directly (byte parser `read_edf_sf_highpass`, the same header approach as `1_inspect_edf`, **no MNE**) and reports, grouped by unique value with file counts, each file's **sampling frequency** and **acquisition high-pass** **restricted to the selected channels** (the config `remap` keys, matched against the header channel labels; falls back to all channels for a file absent from the config). The value shown is **all distinct values** among the kept channels (joined by ` / `), **not** the most common one — so a montage whose EEG channels ended up at **different sampling frequencies (or high-passes) because of a bad export** is made visible, not hidden: such a mixed entry is rendered in red with a ⚠ and a "check the export!" note. Sampling frequency comes from the montage channels (MNE's file-level `info['sfreq']` would instead be the *max*, biased by faster non-EEG channels such as a 512 Hz ECG), and the high-pass is the `HP:` value of the `prefiltering` header field (`none/DC` when absent). Two optional transforms follow, applied per file **in this order**: (1) **resampling** (`raw.resample(target_freq, npad='auto')`, default **OFF**, `cb_resample`/`txt_target_freq`, guarded to never upsample) then (2) **high-pass** (`raw.filter(l_freq=…, h_freq=None)`, default **OFF**, `hp_check`/`hp_freq`; used to harmonise a heterogeneous dataset to a common corner — choose a target ≥ the max acquisition high-pass shown — or to centre DC-coupled data). Resampling is applied **before** the high-pass: `raw.resample` already anti-aliases (FFT method), so resampling first introduces no aliasing, and any resampling edge transients stay at the recording ends (wake). The order does **not** change the result — both transforms are linear, so the signal interior is identical to <0.2 µV either way (verified empirically), and the resampler's anti-alias step rings identically on an abrupt amplitude change regardless of order (a low-corner high-pass does not touch the sharp edge that produces the ring). High-passing the already-downsampled signal is chosen purely because it runs on far fewer samples (markedly faster). Both are analysis-only — no filtered/resampled signal is saved, only the resulting metrics/plots reflect them. Off by default keeps results byte-identical. On the **Curry twin** the high-pass defaults **ON** (0.1 Hz, DC-coupled data has no hardware high-pass) and the acquisition scan reads the Curry header via MNE (`read_raw_curry`), reporting the high-pass as `none/DC`.

**Per-participant pipeline progress bar**: below the participants bar (`i/N`) a second bar spans the *current participant's whole pipeline* in arbitrary "time-cost" units (`COST_LOAD_DATA` + optional `COST_HIGHPASS`/`COST_RESAMPLE`, `COST_ANALYSE_CH`×channels, `COST_RENDER_CH`×channels). It advances continuously (no per-channel reset) through three phases — **Load** (data → resample → high-pass sub-ticks), **Per-channel** (one tick per analysed channel), **Report** (one tick per rendered channel) — with a 3-segment legend **glued directly under the bar** (`VBox([progress_ch, phase_legend])`, the legend drawn with no top border and square top corners so it reads as the bar's own labelled track) whose segment widths are proportional to those costs (≈ each phase's share of the run time). The `progress_ch_label` also names the current activity (`loading EDF signal…`, `high-pass 0.5 Hz…`, `analysing C3…`, `building report…`).

**EDF scan**: recursive (`rglob('*.edf')`), so datasets organized in subfolders (e.g. `group1/`, `group2/`) are fully covered without needing to run the tool per subfolder.

**Hypnogram suffix auto-detection** (see *Cross-cutting procedures*): tools 5 and 6 prefer the **longest** suffix among candidates appearing for ≥ 50% of the maximum count (more specific = remapped/processed version; count breaks ties on equal length). The 50% threshold prevents a rare accidental long suffix from winning when remapping is far from complete.

**Live output warnings**: if a hypnogram is not found, fails to load, has a length mismatch with the EDF, or contains unrecognised stage labels, a plain-text `⚠` warning is printed in the notebook output area immediately after the per-participant result line (in addition to the yellow banner already shown in the HTML report's Overview section).

**Hypnogram label validation**: after mapping labels to YASA integers (`W→0`, `N1→1`, `N2→2`, `N3→3`, `R→4`), the tool checks for unrecognised labels. `MT` (movement time) is silently tolerated — YASA treats it as an artifact epoch (NaN). Any other unrecognised label triggers a warning. Two severity levels: if > 10% of epochs are unrecognised, `hypno_vec` is set to `None` (spectrogram skipped, error-level warning pointing to `3_remap_hypno_voila`); if ≤ 10%, `hypno_vec` is kept but a warning lists the unrecognised labels and their count. This catches hypnograms that were not yet remapped to AASM convention (e.g. `S1/S2/S3/S4` or raw numeric labels).

**Custom (non-AASM) stages** (see *Cross-cutting procedures*): an editable `Custom stages` field (auto-filled from `config_param/custom_stages.json`) registers project-specific labels. Registered labels are treated as **recognised** — they no longer count toward the >10% skip and are not set to `None` — the YASA hypnospectrogram is replaced by the shared `plot_hypnospectrogram()` (custom labels stacked below N3, their own colours), and `quality_summary_by_stage.tsv` + the `dataset_overview.html` by-stage tables/boxplot insets iterate `['W','N1','N2','N3','R'] + custom_stages` (`stage_map` extends the YASA codes with `5+i`).

**Spectrogram fault tolerance**: the `yasa.hypno_upsample_to_data()` + `yasa.plot_spectrogram()` calls are wrapped in a `try/except`. If YASA raises an exception, the channel's spectrogram section in the HTML report shows the error message and processing continues to the next channel and participant without interruption.

**Report structure**: the MNE Report has one section per channel (e.g. "C3", "Fp1") plus an "Overview" section with the flag summary. The Overview section contains: flag summary → butterfly time series → PSD overlay → **across-channel median spectrogram** → amplitude-distribution overlay → averaged metrics table → **pooled metric distribution across channels** → **metric-by-electrode scatter**. Each channel section contains: histogram → time series → metrics table → spectrogram. The selector widget shows how many participants already have an existing report before the "Skip participants with an existing report" checkbox.

**Per-participant "All electrodes" figures** (`All electrodes — …`, Overview section): the per-participant analogue of `dataset_overview.html`'s pooled/grouped boxplots, restricted to the single recording's own channels so an outlier channel is spotted at a glance. Two figures, both over `OVERVIEW_KEY_METRICS` (`std_uV, flat_pct, bounds_pct, hist_extreme_pct, p99_abs_uV, p999_abs_uV`), each subplot one metric, each dot coloured **red = channel flagged**, grey = within thresholds:
- **Pooled metric distribution** (`plot_participant_pooled_boxplots`): one box per metric across the participant's channels (each dot = one channel), with a **per-stage inset** (box across channels per sleep stage, stage colours) when the participant has a valid hypnogram — built from the `file_stage_rows` collected during the run.
- **Metric-by-electrode scatter** (`plot_participant_electrode_scatter`): x-axis = electrode in **montage order**, one point per electrode (a boxplot would degenerate since each electrode carries a single value here). Purpose is fast outlier spotting, not distribution shape.

Both are built from the per-channel `file_rows` already collected for `quality_summary.tsv` (no recomputation), returned as matplotlib `Figure`s (added via `report.add_figure`, unlike the dataset-level functions which embed base64 `<img>`), and wrapped in non-fatal `try/except`. `KEY_METRICS`/`LABELS` were lifted to the module-level `OVERVIEW_KEY_METRICS`/`OVERVIEW_METRIC_LABELS` so `generate_dataset_overview` and these two functions share one source of truth. **Curry twin**: the generator (`_make_tool5_curry.py`) removes `bounds_pct` from `OVERVIEW_KEY_METRICS`/`OVERVIEW_METRIC_LABELS` (DC-coupled data has no EDF bounds), so the per-participant figures drop the bounds subplot too; the functions and their Overview wiring otherwise pass through unchanged. Because the Curry generator addresses cells by index, this code lives **inside the existing shared-functions cell** (no new notebook cell) to keep the indices stable.

**Across-channel median spectrogram** (`All channels — spectrogram`, Overview section): summarises the night's spectral structure for the whole montage — useful on dense 32/64-channel montages where reading one spectrogram per channel is impractical. **All channels are kept, flagged ones included**, like every other Overview plot. The other Overview plots are *overlays*, where a bad channel is an extra aberrant curve *beside* the good ones and the reader sees both; a spectrogram average is a single *aggregate* image with no such escape, so the robustness has to come from the estimator: the **median across channels of the dB spectrogram** tolerates up to half the channels being aberrant, which is what makes "keep everything" compatible with "still readable". Deliberately **not** normalised per channel (a dead channel has ~zero variance, so a z-score would explode instead of being absorbed) and taken in **dB, not linear power** (a linear mean is driven by the noisiest channels; the median of log-power is the scale-robust quantity, and absolute dB is the meaningful unit when inspecting *raw* data). The figure title reports the channel count and how many of them are flagged. Regional (frontal/central/posterior) panels were considered and rejected: the per-channel spectrograms already carry the topographic detail, and region assignment by channel name fails silently on non-10-20 exports.

- **`Sxx` computed once, used twice**: the median needs every channel's spectrogram simultaneously, so the stack is held in memory anyway (**~260 MB float32** on 58 ch / 7.9 h / 512 Hz; ~295 MB on 44 ch / 11.8 h / 1024 Hz; ~410 MB on 77 ch / 9.4 h / 2048 Hz). Given that, a dedicated pass placed **just before the Overview section is built** computes each channel's dB spectrogram into a `ch_sxx` cache, and the per-channel report figures **reuse** it — recomputing would cost ~8 min per file on a dense montage for no benefit. `plot_hypnospectrogram()` is therefore split into **`compute_hypnospectrogram_sxx()`** (signal → `(f, t, Sxx_dB)`) + **`render_hypnospectrogram()`** (`(f, t, Sxx_dB)` → Figure), with `plot_hypnospectrogram()` kept as a thin wrapper with an **unchanged signature**. Verified: the original function, the new wrapper, the explicit compute+render pair, and a float32 round-trip all produce **bit-identical PNGs** (float32 costs ≤ 1.5e-5 dB, four orders of magnitude below the ~0.2 dB colormap step). Failure to compute one channel's `Sxx` is non-fatal — that channel's report section shows a "computation failed" note and the median simply excludes it.
- **Frequency resolution is set by the 30 s window (1/30 Hz), not by `sfreq`**: the 0.5–40 Hz band always yields ~1186 bins, so **resampling shrinks the compute time (~6× from 1024 → 256 Hz) but not the stack size**. Frequency binning is the fallback if a much denser montage ever makes the stack too large.
- **Progress bar**: the spectrogram cost moves from the *Report* phase to the middle phase, which now covers per-channel metrics **then** the spectrogram pass (`COST_ANALYSE_CH = 1`, new `COST_SPECTRO_CH = 3`, `COST_RENDER_CH` lowered 3 → 2). The 3-segment legend still matches the three chronological phases.
- **Not propagated to tool 8**: `8_SpectralPower` keeps the monolithic `plot_hypnospectrogram()`. Behaviour is identical, but the two copies now differ *structurally* — when syncing a change to the plotting logic, apply it to tool 5's `render_hypnospectrogram()` body and to tool 8's single function.

**Metrics table — interpretation column**: for `Flat signal (%)`, `At EDF bounds (%)`, and `Extreme histogram (%)`, the interpretation cell shows the actual threshold value used for that run (read from the widget at run time, e.g. "flag if > 3.5%"), making each report self-documenting. The three threshold-based entries are generated dynamically inside `run_analysis` via a local `interpretations` dict that overrides the static `METRIC_INTERPRETATIONS` dict.

**Histogram axis scaling**: the shared Y-axis for histograms is computed from non-suspect channels only (those not flagged for `flat_pct` or `std_uV`). A flat or dead channel concentrates all samples in 1–2 bins and would otherwise crush healthy distributions. Fallback to global max if all channels are suspect. The shared X-axis is `±x_lim_hist` where `x_lim_hist = p99.9 of |amplitude| across all channels` (no cap on EDF; capped at `DISPLAY_YLIM_UV = 500 µV` on the Curry twin) — this allows cross-channel comparison: a clipped channel appears compressed, a dead channel appears as a narrow central spike.

**Distribution histogram — robust range for peak detection**: the histogram feeding the Savitzky-Golay smooth + `find_peaks` is bounded to robust percentiles (`HIST_CLIP_PCT = (0.5, 99.5)`), applied uniformly on EDF + Curry. Rare extreme samples (DC drift / non-clipped artifacts on DC-coupled Curry data; occasional spikes on EDF) otherwise stretch the amplitude range, inflating the Scott bin count and the `0.02·n_bins` smoothing window so the SG curve over-smooths into a dome or a flat line inside the display zoom. Clipping concentrates the bins on the bulk of the signal so the curve tracks the distribution and resolves genuine modes. `hist_extreme_pct` is deliberately kept on a **separate full-range histogram** (it is a clipping proxy on the outermost bins), so `flat_pct`/`bounds_pct`/`hist_extreme_pct` and all stats stay unchanged — **only `n_peaks` becomes more sensitive** (it now resolves bimodal/irregular channels the former over-smoothing hid). Degenerate/flat channels (`std < 1e-10` or `< 3` unique values) skip the clip. Because the bin count is set by Scott's rule on the *in-range* samples, its resolution scales with the recording's sample count and the shape of its bulk distribution — a higher-variance channel (e.g. a residually noisy DC recording after the high-pass) gets fewer, wider bins than a cleaner one. On the Curry twin the high-pass (0.1 Hz, ON by default) removes the DC drift before this metric, and the code passes through `_make_tool5_curry.py` unchanged.

**Flagging criteria and default thresholds:**

| Metric | Threshold | Detects |
|--------|-----------|---------|
| `flat_pct` | > 3.5% | Flat segments, dead channels |
| `bounds_pct` | > 1.0% | Saturation at EDF physical-range limits |
| `n_peaks` | ≥ 2 | Irregular / suspect amplitude distribution (bimodal or multimodal — inspect the channel) |
| `std_uV` | < 5.0 µV | Dead / near-dead channels |
| `hist_extreme_pct` | > 1.0% | In-range clipping (saturation within declared EDF range) |

`flat_pct` = fraction of consecutive sample pairs with \|diff\| < `max(2×ADC_step, 0.06 µV)`. `bounds_pct` = fraction of samples within 0.5 µV of the EDF physical_min/max header limits, reconstructed in µV by `get_phys_bounds_uV()`. **Unit handling**: MNE keeps `_raw_extras['physical_max']`/`['offsets']` in each channel's *native* EDF physical unit (µV, mV, V…), not in µV, whereas the signal is read as `raw.get_data() * 1e6` (always µV). `get_phys_bounds_uV()` therefore multiplies both bounds by `extras['units'][ch_idx] * 1e6` (MNE's native-unit→volts factor: `1e-6` µV, `1e-3` mV, `1.0` V) so the comparison is unit-consistent. This is essential for non-EEG channels: Compumedics EOG/EMG/ECG are declared in **mV** (e.g. `physical_max = 1.0 mV`), so without the conversion they are compared against a 1.0 µV bound — 1000× too small — and `bounds_pct` flags ~98–100 % of a perfectly healthy signal. EEG channels (declared in µV) are unaffected. `hist_extreme_pct` = fraction of samples in the outermost histogram bins (computed on a **full-range** histogram, independent of the percentile-clipped histogram that feeds the Savitzky-Golay curve + peak detection — see *Distribution histogram — robust range for peak detection* above). **Kurtosis is intentionally NOT used as a flagging criterion** — normal PSG EEG has physiologically high kurtosis (spindles, K-complexes produce values of 100–500), making it unreliable without per-subject normalization. `p99_abs_uV` and `p999_abs_uV` (99th and 99.9th percentile of |amplitude|) are recorded as informational metrics in `quality_summary.tsv` but are not currently used for flagging; they are useful for cross-channel and cross-dataset amplitude comparison.

**`dataset_overview.html` — dataset-level summary**: generated at the end of every run from the full cumulative `quality_summary.tsv` (reflects all participants processed to date, not just the current run). Contains two levels:
- **Global section (all electrodes pooled)**: stats table (mean / median / p5 / p25 / p75 / p95 per metric), **mean (median) by sleep stage table** (metrics as rows — `mean_uV`, `std_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `p99_abs_uV`, `p999_abs_uV` — stages W/N1/N2/N3/R as columns; shown only when `quality_summary_by_stage.tsv` is present), n_peaks frequency table by electrode (flags DC drift and quantization cases), pooled boxplots (one subplot per key metric; each subplot contains a **stage inset** in its upper-right corner showing per-stage boxplots with stage colours W/N1/N2/N3/R — inset shown only when stage data is available; stage colours: W `#969696`, N1 `#9e9ac8`, N2 `#807dba`, N3 `#6a51a3`, R `#c994c7`), grouped boxplots (one subplot per key metric, x-axis = electrode — compares electrodes side by side).
- **Per-electrode sections**: stats table for that electrode only, boxplots of each metric's distribution across participants with individual data points overlaid (reveals outlier participants for that channel).

Key metrics shown in plots: `std_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `p99_abs_uV`, `p999_abs_uV`. All numeric metrics (`mean_uV`, `kurtosis`, `skewness` included) appear in the stats tables. HTML is standalone (figures embedded as base64 PNG, no external dependencies).

**Outputs per run**:
- `<data_folder>/reports_quality_overview/<relative_subfolder>/<file_id>_quality_overview.html` — HTML reports mirror the EDF subfolder structure under `reports_quality_overview/`
- `<data_folder>/reports_quality_overview/<relative_subfolder>/<file_id>_quality_metrics.tsv` — the **per-file** numeric-metrics table (same columns as `quality_summary.tsv`), written next to the HTML **before** it (so "report exists ⇒ data exists"). This is the durable per-item data and the skip gate's *data present* marker (see *Skip + cumulative-merge*).
- `<data_folder>/reports_quality_overview/<relative_subfolder>/<file_id>_quality_by_stage.tsv` — the **per-file** by-stage table (same columns as `quality_summary_by_stage.tsv`); written only when the participant has a valid hypnogram.
- `<data_folder>/reports_quality_overview/quality_summary.tsv` — all numeric metrics for all channels, cumulative across runs (always at root); columns: `file_id`, `channel`, `mean_uV`, `std_uV`, `kurtosis`, `skewness`, `p99_abs_uV`, `p999_abs_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `n_peaks`, `suspect_reason`, `exclude` (last two columns). **Regenerated each run from all `*_quality_metrics.tsv` on disk** (`rglob`, sorted by `file_id, channel`) — not an in-memory merge — so a file processed earlier and skipped now is never dropped (interruption-safe).
- `<data_folder>/reports_quality_overview/dataset_overview.html` — dataset-level statistics and plots (always at root, regenerated each run)
- `<data_folder>/reports_quality_overview/quality_summary_by_stage.tsv` — key metrics split by sleep stage, one row per `file_id × channel × stage`; columns: `file_id`, `channel`, `stage`, `mean_uV`, `std_uV`, `flat_pct`, `bounds_pct`, `hist_extreme_pct`, `p99_abs_uV`, `p999_abs_uV`. Populated only for participants with a valid hypnogram. **Regenerated from all `*_quality_by_stage.tsv` on disk** (same rule as `quality_summary.tsv`, sorted by `file_id, channel, stage`).
- `<data_folder>/reports_quality_overview/failed_files.tsv` — files that could not be read (at root)

**End-of-run summary** (printed in the notebook output): participants processed, participants with ≥1 flagged channel, total flagged channels, files failed to load, path to `dataset_overview.html`.

### 6. Preprocessing + epoch rejection (`6_preprocessing_voila.ipynb`)

*(Stable tool — formerly “Phase 2” of the preprocessing pipeline.)*

Implemented as a Voila notebook with four sections: (1) path configuration, (2) preprocessing and rejection parameters, (3) participant selection, (4) processing loop. **Section 2 is organised into two headed sub-sections**: **Preprocessing** (resampling + notch + bandpass filter) and **Epoch rejection** (peak-to-peak amplitude per stage, flat signal & gradient, 1/f fit quality, and the optional event-based rejection — see below).

**Inputs**:
- `quality_summary.tsv` from Phase 1 — **optional** (strongly encouraged); its `exclude` column pre-fills the per-channel drop selection. Its `FileChooser` shows only `quality_summary.tsv` (`filter_pattern`), not the per-file `*_quality_metrics.tsv` that quality_overview also writes in that folder.
- `remap_reref_persubject.json` from `2_select&remap_channels_edf` — drives channel remapping and re-referencing per participant, **and is the fallback channel source** for participants absent from `quality_summary.tsv` (see *Section 3 discovery* below)
- Raw EDF files and remapped hypnograms (default suffix `_Hypnogram_remapped.txt`)

**Section 3 participant discovery (EDF folder = ground truth)**: the participant list is built from the **EDF files on disk** (`fc_edf` folder, recursive), *not* from `quality_summary.tsv` — so a recording present on disk with a config entry is never hidden just because it is missing from the QC file. Per participant, the channel list + `exclude` pre-selection comes from `quality_summary.tsv` if present, else from the config **`remap` values** (remapped names, nothing pre-excluded), else empty (shown, flagged *not found in JSON config*, unchecked). `quality_summary.tsv` is optional (removed from the Run guard); when absent or partial, `load_info` surfaces `⚠` warnings for every discrepancy: EDF-count vs quality-count, EDFs missing from quality (loaded from config), EDFs missing from config (unprocessable), and stale quality entries with no EDF on disk. All cross-source id matching is `os.path.normcase`-wrapped at the comparison only.
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
1. **Resampling** — `raw.resample(target_freq, npad='auto')`. Target frequency chosen by user; applied **before** filtering. `raw.resample` already anti-aliases (FFT method), so resampling first introduces no aliasing; the order is chosen because filtering the already-downsampled signal runs on far fewer samples (faster) — not to prevent aliasing (both steps are linear, so the order does not change the result). Step is skipped if checkbox is unchecked.
2. **Re-referencing** — applied as specified in JSON config per participant: `'average'` → common average reference; `[list]` → subtract listed channel(s) then drop them; empty → no re-referencing.
3. **Notch filter** *(optional, **ON by default**)* — removes power-line noise via `raw.notch_filter(freqs=notch_freq_val)` using **MNE's default method** (FIR; `method=` is left unset). Single editable frequency (`cb_notch` / `txt_notch_freq`, default **50 Hz**). Applied **after re-referencing, before the bandpass** ("notch then band-pass" convention). **Fatal** on failure (like the bandpass step). Defaults ON because a power-line notch is the norm for most databases; untick `cb_notch` to reproduce the former byte-identical no-notch output. The frequency box's initial visibility follows `cb_notch.value` (checkbox-revealed-widgets rule), so the pre-ticked box shows correctly.
4. **Bandpass filter** — FIR zero-double-pass Hamming window, defaults `l_freq=0.1 Hz, h_freq=50 Hz`. Applied via `raw.filter(..., method='fir', phase='zero-double', fir_window='hamming', fir_design='firwin')`.

**Epoching (configurable epoch length)**: epochs are created with `mne.make_fixed_length_epochs(raw, duration=epoch_sec)`. `epoch_sec` is user-selectable via a Section-2 **`Epoch length`** dropdown (`dd_epoch_len`) offering the **divisors of 30 ≥ 5 s** — `30 (classic, default) / 15 / 10 / 6 / 5` — so the classic 30 s scored epoch re-cuts into a whole number of sub-epochs and the 4 s Welch window used by the 1/f fit still fits. The hypnogram is scored at 30 s (one label per 30 s epoch): the 30 s length validation/trim against the recording is **unchanged**, then each 30 s label is **expanded** to fill its sub-epochs — `hypno_sub = np.repeat(expert_hypno, epoch_factor)` with `epoch_factor = 30 // epoch_sec` — so each sub-epoch simply inherits its parent 30 s stage (no re-scoring, no interpolation). All flagging/summary/plot code keys on `hypno_epochs`/`n_epochs` generically, so it adapts automatically; the Welch window is `n_per_seg = int(min(4, epoch_sec) * sf)` (byte-identical to `4·sf` for every allowed size). Amplitude p-p thresholds stay the **same absolute µV** across sizes (artefact p-p does not scale with window length). Default 30 s keeps every output byte-identical to the pre-feature tool. Sleep stage assigned to each epoch from the (expanded) hypnogram; epochs at the tail beyond the hypnogram length are discarded.

**Event flagging when epoch length ≠ 30 s** (see the event method below): the onset-only rule is **kept** (`compute_event_epoch_mask` takes `epoch_sec`; an event flags only the sub-epoch containing its `Start`), and a **prominent amber warning** (`epoch_warn`, shown whenever `dd_epoch_len.value != 30`, plus a per-file `⚠` printed in the run log when event flagging is on) reminds the user that — because annotated durations are unreliable — an event spanning several sub-epochs leaves the others un-flagged, a risk that grows as epochs shrink.

**Epoch rejection — six methods, each individually enable/disable-able** (the four signal-quality methods default ON, *Event overlap* is optional and off by default; see *Per-method enable/disable* below):

**Vocabulary — flagging, not rejection**: tool 6 only *flags* (epoch, channel) pairs; the actual **rejection**
(and channel interpolation) is delegated to the downstream decision tool (`7bis_reject_automatically`, see
§8.1 of the task plan). All **user-facing** strings (UI section headers, HTML report titles, progress phases,
prints) therefore say *flagged/flagging*; the **machine** identifiers are left unchanged for backward
compatibility (`.fif` `reject_flag`/`reject_method`, TSV columns `n_rejected`/`pct_rejected`/`reject_any`,
function names `compute_rejection_masks`/`build_rejection_summary`, file names `*_rejection*.tsv`, sidecar key
`rejection_thresholds`).

All methods operate on the epoch data in µV. The time-domain methods (amplitude / flat / gradient) read the signal **one channel at a time** from `epochs` (`epochs.get_data(picks=[ci])[:, 0, :] * 1e6`) so a high-density montage is never materialised as a second full `n_epochs × n_channels × n_times` float64 array (see *Per-channel rejection (memory)* under Curry 9 support — the same code runs in EDF and Curry); the 1/f method works on the small precomputed Welch PSD array. Flagging masks are boolean arrays of shape `(n_epochs, n_channels)` — a `True` entry means that (epoch, channel) pair was flagged. These per-(epoch, channel) masks are now **persisted** to `{file_id}_epoch_channel_rejection.tsv` (the fine-grained source of truth for the decision tool); the `.fif` per-epoch `reject_flag` still uses the naive rule *an epoch is flagged if any channel is flagged by any method* — a default flag, **not** the final rejection. The per-stage / global **flagging percentages are computed over (epoch × channel) pairs** (`pct = flagged_pairs / (n_epochs_in_stage × n_channels)`), replacing the former "≥1 bad channel ⇒ whole epoch" rate that saturated at ~100% as soon as every epoch had one bad channel. A single ordered registry — `METHOD_ORDER = ['amplitude', 'flat', 'gradient', '1f_error', '1f_r2', 'event']` with `METHOD_CODE`/`MULTIPLE_CODE` — is the one source of truth shared by the mask builder, heatmap, per-epoch log, per-(epoch, channel) log, per-stage summary and global table. Every consumer selects `[m for m in METHOD_ORDER if m in mask_dict]`, so **any method the user disables** (see *Per-method enable/disable* below) — like `event` when off — is simply absent from `mask_dict`, `methods_run`, and every downstream column/row; older or method-subset outputs keep their own column set.

**Per-method enable/disable** (checkboxes `cb_amplitude` / `cb_flat` / `cb_gradient` / `cb_1f` / `cb_event_reject`): each flagging method is its own UI subsection — a checkbox, a grey one-line description, and its threshold widget(s), the latter shown/hidden from the checkbox (checkbox-revealed-widgets rule, initial `display` derived from the value). The four signal-quality methods default **ON** (a default run stays byte-identical to before), `event` **OFF**. A disabled method is not computed and not added to `mask_dict` (so it vanishes from every table, exactly like `event`), and `methods_run` in the sidecar records the true subset. **Unticking `cb_1f` skips `compute_psd` entirely** (the caller passes `psds=None`) — the whole Welch-PSD + specparam fit is the dominant cost on short epochs, so this is the main speed-up (the fit scales ~linearly with epoch count *and* each fit is slower on the noisy short-epoch PSD). One checkbox `cb_1f` governs **both** `1f_error` and `1f_r2` (they share the single fit). A **guard** blocks the run (red message) if no method is ticked. `1f_error` and `1f_r2` remain gated by PSD presence (`do_1f = psds is not None`); the three time-domain methods take `do_amplitude`/`do_flat`/`do_gradient` args in `compute_rejection_masks` (default `True`, so the Curry twin and any old caller pass through). **Mixed-method safeguard**: at run start tool 6 reads `methods_run` from every existing `{file_id}_preprocessing_params.json` under `derivatives/` and, comparing the **signal-quality** subset only (`event` excluded — its presence legitimately varies per file with the annotations), shows a **non-blocking amber warning** if any already-processed file used a different method set, since mixing sets makes the per-method rates in `global_rejection_by_stage.tsv` non-comparable across files (see that file's `n_part_{m}` canary below).

| Method | Signal feature | Default threshold | Notes |
|--------|---------------|-------------------|-------|
| **Amplitude** | Peak-to-peak = `max(epoch) − min(epoch)` | W: 300, N1: 250, N2/N3: 200, REM: 250 µV | Per-stage threshold; W/REM more lenient because muscle and eye-movement artefacts are physiologically common in those stages. Equivalent to MNE's `drop_bad(reject=...)` criterion. |
| **Flat signal** | Peak-to-peak < threshold | 1 µV | Detects disconnected electrodes or amplifier saturation within a single epoch. Logically identical to MNE's `drop_bad(flat=...)` criterion: both compare `ptp` against a low-amplitude threshold. |
| **Gradient** | `max(|diff(epoch)|)` across time | 100 µV/sample | Maximum sample-to-sample absolute difference; sensitive to sudden jumps, electrode pops, and movement artefacts not captured by peak-to-peak. `diff` and `max` both operate on `axis=-1` (time axis) of the per-channel `(n_epochs, n_times)` slice. |
| **1/f fit quality** | Specparam aperiodic fit on Welch PSD (4 s windows, **configurable fit range, default 2–45 Hz** — see below, `aperiodic_mode='fixed'`, `peak_width_limits=[0.5, 20]`, `min_peak_height=0.3`, **`max_n_peaks=8`**) | MAE > 0.15 OR R² < 0.95 | Fit lower bound ≥ 2 Hz limits slow-wave influence. A **capped peak model is used** (`max_n_peaks=8`): periodic components (spindles, alpha…) are still modelled and removed *before* assessing the aperiodic fit quality, but the peak search is **bounded** so a noisy short-epoch PSD (e.g. 6 s epochs → only ~2 Welch segments) cannot send specparam into a runaway peak-fitting loop that made each fit ~2–3× slower. 8 ≥ the realistic neural peak count on 2–45 Hz, so this is **not** the `max_n_peaks=0` case, which would push all peak power into the aperiodic component, degrading R² and over-rejecting nearly every N2/REM epoch. **⚠ Not byte-neutral**: capping at 8 changes MAE/R² (hence the 1/f flag) wherever specparam previously fit > 8 peaks — on clean 30 s PSDs it usually fits few, so existing outputs are *nearly* unchanged, but **reprocess a dataset once to homogenise**. **Sync constraint**: the `SpectralModel(...)` lives in **three copies that must change together** — tool 6 EDF (`compute_rejection_masks`), its Curry twin (regenerated by `_make_tool6_curry.py`), and `qc_rejected_epochs_lib.py` (tools 7/7bis recompute attribution with it) — else tool 7's recomputed attribution diverges from tool 6's flags. Metrics read via `get_metrics('error','mae')` / `get_metrics('gof','squared')` (specparam 2.x). A failed fit is treated as a double flag (both error and R²). |
| **Event containment** *(optional, off by default)* | 30 s epoch **containing the onset** of any **selected** canonical scored-event type (arousal, apnea, hypopnea, limb movement, SpO2 desaturation…) | **onset-only** — the epoch holding the event `Start`; the annotated `Duration` is **intentionally ignored** (clinicians often score only the onset without a reliable duration), so each event flags exactly one epoch | **Epoch-level** flag, replicated across all channels → single `flag_event` column. Events read with the shared CSV-first / XML-fallback `load_events(edf, csv_suffix)`; raw labels mapped to canonical via `event_remap.json` (tool 4). UI: a checkbox to activate, **a wrapping row of checkboxes for the canonical types** (all shown at once, populated from the chosen `event_remap.json`), and an inline note explaining the onset-only rule so the choice is informed. A **"Count affected epochs"** button reports, over the participants currently checked in Section 3, how many epochs each selected type would flag (overall + per stage) using only the hypnogram length/stages and events — no signal is read. Missing/unreadable event companions are non-fatal (the file keeps the other 5 methods, never added to `failed`). |

**Configurable 1/f fit range**: the aperiodic fit window is user-editable via two widgets
(`txt_1f_fmin` / `txt_1f_fmax`, **default 2–45 Hz**); the fit uses
`freq_mask = (psd_freqs >= fit_fmin) & (psd_freqs <= fit_fmax)`. The Welch PSD ceiling **follows the fit
max** (`fmax_psd = min(fit_fmax_val, sf/2 - 0.5)`, was hardcoded 30 Hz), and the default bandpass is
**0.1–50 Hz** so the whole fit band is preserved. The range is persisted in the sidecar as
`rejection_thresholds.1f_fit_range_hz` and **read back by tool 7** (`qc_rejected_epochs_lib.load_params`
→ `info['fit_range']`, threaded into `compute_psds(fmax=)` / `fit_1f(fmin=)`; fallback `(2.0, 45.0)`), so
tool 7's recomputed per-channel attribution matches tool 6's. Wired across tool 6 (EDF + Curry),
`qc_rejected_epochs_lib.py`, and the tool-7 batch + Voila.

**Custom (non-AASM) stages** (see *Cross-cutting procedures*): an editable `Custom stages` field (Section 1, auto-filled from `config_param/custom_stages.json`) extends the per-stage logic. Each custom stage gets its **own amplitude-threshold widget** (default 250 µV, generated dynamically when the field changes) feeding `ptp_thresholds`; the per-participant summary, `global_rejection_by_stage.tsv`, the heatmap hypnogram strip and the **"Count affected epochs"** estimate all iterate `['W','N1','N2','N3','R'] + custom_stages`. Custom-stage epochs are still rejected by the other (stage-independent) methods regardless.

**Heatmap** (rendered into `{file_id}_preprocessing_report.html`, **not** saved as a standalone PNG): channels (Y-axis) × epochs (X-axis); each cell coloured by the flagging method with priority encoding when multiple methods fire. A hypnogram strip is drawn above the main heatmap. Colour scheme: dark purple = none, red = amplitude, blue = flat, orange = gradient, yellow = 1/f error, green = 1/f R², **magenta = event**, dark red = multiple. Title shows the overall **pair-based** flagging percentage (`(combined_matrix > 0).mean()`, i.e. fraction of (epoch × channel) cells flagged — matches the per-stage table denominator).

**Two-step QC approach** — Phase 2 does **not** drop epochs. It saves ALL epochs (including flagged ones) with an MNE `metadata` DataFrame attached, so downstream Phase 2b can inspect rejected epochs before finalising the rejection.

**Per-participant pipeline progress bar**: below the participants bar (`i/N` + the current phase from `set_phase`, e.g. `loading EDF…`, `resampling…`, `epoching…`, `flagging epochs…`, `building report…`) a second `IntProgress` (`progress_step`) spans the *current participant's whole pipeline* in arbitrary "time-cost" units (`COST_LOAD` + optional `COST_RESAMPLE`/`COST_NOTCH`/`COST_FILTER` for the *Load* group, then `COST_EPOCH`, `COST_REJECT`, `COST_SAVE`, `COST_REPORT`). It is sized per participant (the optional steps only count when enabled) and advanced **cumulatively** at each phase (`_reject_base` marks where the *Flag* segment begins), with a **4-segment legend** (`build_step_legend` → **Load / Epoch / Flag / Report**, Report = save+report) glued directly under the bar (`VBox([progress_step, step_legend])`, no top border, square top corners) whose widths are proportional to those costs. Mirrors tool 5's pipeline bar; UI-only, so all outputs stay byte-identical. Its *Reject* segment is **animated** by the per-channel/per-epoch `compute_rejection_masks` callback (`_rej_progress` → `progress_step`, per channel for the time-domain pass then per epoch for the 1/f fit, filling `[_reject_base, _reject_base + COST_REJECT]`), so a slow high-density participant is not mistaken for a crash — **identical in EDF and the Curry twin**, since both now run the rejection per channel (see *Curry 9 support → Per-channel rejection*). **Sync constraint**: this bar and the callback live in the EDF original and pass through / are string-matched by `tools_curry/_make_tool6_curry.py` — re-run that generator after editing the widgets / display / reject call site.

**Outputs per participant** — the `.fif` + its params sidecar + the optional context companion under `<output_folder>/derivatives/raw_epo/<edf_subtree>/` (a **sibling of `clean_epo_manual/`** written by tool 7 and **`clean_epo_auto/`** written by tool 7bis, so the three epoch stages sit side by side under `derivatives/`); the per-file TSV/HTML reports under `<output_folder>/reports_preprocessing/<edf_subtree>/` — **both mirror the EDF subfolder tree** (as `5_quality_overview` does for `reports_quality_overview/<edf_subtree>/`; an EDF sitting in the data-folder root has `<edf_subtree> == '.'` so its reports collapse back to the flat root, keeping subfolder-less datasets byte-identical to the pre-mirror tool). The skip-checks and the global-table rebuild therefore search **recursively** (`reports_dir.rglob(...)`), so a report written in a subfolder is still found. When a file previously processed by the old **flat** layout is reprocessed into its subfolder, tool 6 first deletes any stale **root-level** copies of that file_id's three report files (only when `<edf_subtree> != '.'`) so the globbed globals never double-count it:
- `{file_id}_all-epo.fif` — all epochs with `epochs.metadata` DataFrame (columns: `epoch_idx`, `stage`, `reject_flag`, `reject_method`, `flag_amplitude`, `flag_flat`, `flag_gradient`, `flag_1f_error`, `flag_1f_r2`, plus `flag_event` **when event flagging ran**). The `flag_<method>` columns are **per-epoch "any channel" booleans** (the naive epoch-level decision, read by tool 7); the full per-(epoch, channel) mask **is** now persisted separately, in `{file_id}_epoch_channel_rejection.tsv` (below). The per-epoch/per-stage TSVs and `global_rejection_by_stage.tsv` gain the matching `flag_event` / `event` entries the same way — additively, so event-free runs stay byte-compatible with earlier outputs.
- `{file_id}_preprocessing_params.json` — the resampling / notch / bandpass filter settings (the notch as an additive `notch: {applied, freq_hz}` key, provenance only) + the additive `epoch_length_s` (the selected epoch length; **30 assumed when absent**, so older sidecars stay compatible) + the per-stage rejection thresholds actually used (amplitude p-p per stage, flat, gradient, 1/f MAE/R²) + `methods_run`. Read back by **tool 7** (QC of rejected epochs) to draw threshold reference lines, recompute per-channel margins, and label the epoch length (`qc_rejected_epochs_lib.load_params → info['epoch_length_s']`). Written non-fatally.

  *Migration note*: earlier tool-6 versions wrote the `.fif`/sidecar/context directly under `derivatives/<edf_subtree>/`. When reprocessing such a file, tool 6 deletes the stale same-id copies at that old location so tools 7/7bis (which discover participants by recursively globbing `*_all-epo.fif`) never see the id twice. Tools 7/7bis strip a leading `raw_epo/` component when mirroring the subtree into `clean_epo_*/`, so a pre-`raw_epo` layout still maps to the same `clean_epo_*/<subtree>/`.
- `{file_id}_context-epo.fif` — **optional** EOG/EMG/ECG "context" companion (block `[H]`), written **only** when the participant's `sub_config` carries a tool-2 `context_channels` block. Tool 6 reads just those declared channels from the raw EDF **with `include=` at read time** (exactly like the EEG read — otherwise MNE upsamples every channel in the file to its file-wide max rate, which raised a `bad allocation` on mixed-rate montages), renames them to role labels (`EOG-L`/`EOG-R`/`EMG`/`ECG`), sets MNE channel types, then **unifies the sampling rate to the EEG working rate** (`epochs.info['sfreq']` — native, or the explicit resample target: a single `.fif` requires one rate for all channels, so EOG/EMG/ECG at different native rates are brought to the EEG's rate; display companion only), **display-filters per role** (AASM-like: EOG band-pass 0.3–35 Hz, EMG high-pass 10 Hz, ECG band-pass 0.5–40 Hz — so the tool-7 epoch montage is readable), and epochs them **identically to the EEG** (`make_fixed_length_epochs`, duration 30 s from t=0 → same epoch count regardless of sfreq, so 1:1 index alignment with `{file_id}_all-epo.fif`). Non-fatal: when no context is declared no companion is written and all other outputs stay byte-identical. Read on demand by **tool 7**'s per-epoch view (`load_context_epochs`).
- `{file_id}_epoch_rejection.tsv` — per-epoch flagging table (columns: `file_id`, `epoch_idx`, `stage`, `reject_flag`, one `flag_<method>` bool per method, incl. `flag_event` when event flagging ran). Unchanged (epoch-level "any channel" mirror of the `.fif` decision), kept for human/downstream convenience.
- `{file_id}_epoch_channel_rejection.tsv` — **per-(epoch, channel)** flagging table = the fine-grained **source of truth** for the decision tool. Dense (one row per (epoch, channel) pair, **epoch-major**: all channels of epoch 0, then epoch 1, …), columns: `file_id`, `epoch_idx`, `stage`, `channel`, one `flag_<method>` bool per method (incl. `flag_event` when event flagging ran — **replicated across all channels** of a flagged epoch, since `mask_dict['event']` is tiled), and `reject_any` (OR over the methods present). Self-sufficient: all channels appear for every epoch, so channel list and denominators are derivable without loading the `.fif`. Built by `build_epoch_channel_rejection_log()` (vectorised from the `(n_epochs, n_channels)` masks), written non-fatally before the report. **Per-file only** — no global concatenation (volume). It **also gates the skip check** (see below).
- `{file_id}_rejection_summary.tsv` — per stage per method flagging counts, now **over (epoch × channel) pairs** (columns unchanged: `file_id`, `stage`, `method`, `n_total`, `n_rejected`, `pct_rejected`; **values redefined**: `n_total = n_epochs_in_stage × n_channels`, `n_rejected` = flagged pairs, `pct_rejected` their ratio). Not byte-identical to the pre-change file — reprocess the dataset once to homogenise.
- `{file_id}_preprocessing_report.html` — MNE HTML report with the heatmap and the flagging tables. Its **"Flagged summary"** section opens with three global lines — `Total epochs: N`, `Flagged (epoch×channel): X%` (pair-based), `Epochs with ≥1 flagged channel: Y%` (epoch-level) — then a **per-stage table** titled *Summary per sleep stage and method* (glued grey subtitle *"(% (n) of flagged epoch×channel within each sleep stage)"*). That table has one row per stage (W/N1/N2/N3/R + custom); its **second column is `epochs×channels`** = the any-method pair-based rate for the stage, followed by one column per method (Amplitude, Flat, Gradient, 1/f error, 1/f R², and Event when event flagging ran); every cell shows **`% (n)` over (epoch × channel) pairs** (same denominator as `global_rejection_by_stage.tsv`). When event flagging is enabled, a second **"Event-based flagging by type"** table lists one row per selected canonical event type (+ a bold `(any selected)` union row) with the flagged epoch count, % of all epochs, then `%(n)` per stage — so event types that flag too many epochs can be spotted and de-selected. All tables are report-only (HTML); no TSV schema changes. Built by `build_stage_method_html()` / `build_event_type_html()`, with per-type epoch masks computed in the run loop via `compute_event_epoch_mask(..., [t], ...)`.

**Global output** (always at the `<output_folder>/reports_preprocessing/` **root**, not in a subfolder):
- `global_epoch_rejection.tsv` — concatenation of all `{file_id}_epoch_rejection.tsv` across participants, **regenerated by recursively globbing them from disk each run** (`rglob`, so the per-file tables in the `<edf_subtree>/` subfolders are all found; excluding the global file itself, which shares the `_epoch_rejection.tsv` suffix) — not an in-memory merge, so an interrupted-then-skipped participant is never dropped (see *Skip + cumulative-merge*)
- `global_rejection_by_stage.tsv` — concatenation of all `{file_id}_rejection_summary.tsv` across participants (also `rglob`-ed from disk each run). Counts are **pair-based** (inherited from the per-file summaries — `build_global_summary_table` sums `n_total`/`n_rejected` and recomputes the pooled `pct`). Each method also gets an additive **`n_part_{m}`** column = the number of participants that actually ran method `m`. It is a **canary, not a correction**: `pct_rej_{m}` keeps the pooled denominator `n_total` (summed over **all** participants), so when a method was run on only part of the database `pct_rej_{m}` is **diluted** and `n_part_{m} < n_participants` makes that visible. On a uniformly-processed database `n_part_{m} == n_participants` (constant column, no interpretation change); the mixed-method start-of-run warning above is meant to keep it that way.
- `preprocessing_failed.tsv` — participants that could not be processed (EDF not found, config missing, hypno mismatch, etc.)

### 7. Manually reject flagged epochs — Phase 2b (`7_reject_manually_voila.ipynb`, `7_reject_manually_batch.py`, `qc_rejected_epochs_lib.py`)

Manual quality control of the epochs that `6_preprocessing_voila` flagged: inspect each rejected epoch, override the keep/reject decision, and export a validated `{file_id}_clean-epo.fif`. Reads one participant's tool-6 outputs (`{file_id}_all-epo.fif` + optional `{file_id}_preprocessing_params.json`, and the optional `{file_id}_context-epo.fif` context companion) — **the raw EDF is never reloaded** (the EOG/EMG/ECG context traces come from tool 6's companion `.fif`, not the raw); the analysis channel set is exactly what the `.fif` holds (EEG-only in practice). Tool-6 outputs are never modified. The Voila app is the primary delivery (a code-visible Jupyter twin is planned); the `.py` batch twin produces the Section-2 report over a whole database.

**Shared library (`qc_rejected_epochs_lib.py`)**: the analysis + plotting used by both the notebook and the batch live in one module (imported by both) to avoid drift — `METHOD_ORDER`, the heatmap/method colour palette, the custom-stage helpers, the Welch-PSD config and the specparam 1/f fit are copied from tool 6 and kept in sync. The **per-epoch reject decision is authoritative from `epochs.metadata`**, while the **per-channel attribution** shown here (which channel drove a flag, margins to threshold) is **recomputed** from the signal with the same formulas + the persisted thresholds (fallback: tool-6 defaults when no params JSON is present). Tool 6 now *does* persist the per-(epoch, channel) boolean masks (`{file_id}_epoch_channel_rejection.tsv`), but tool 7 still recomputes because it needs the **continuous metric values** (p-p, gradient, 1/f MAE/R²) for its plots, not just the booleans — a future refactor could let it read the booleans from that TSV for the attribution while still recomputing the values (see task plan §8.3).

**Section 1 — Load**: **two folder pickers** (mirroring 7bis, minus the reports picker — tool 7 never reads `reports_preprocessing/`): an **optional Data folder** (`fc_data`, holds `derivatives/` + `reports_preprocessing/`) that on selection `reset(path=…)`s the raw picker to `<data>/derivatives`, and a **Raw-epochs folder** (`fc_raw`, holds the tool-6 `*_all-epo.fif`, e.g. `derivatives/raw_epo`). Selecting the raw folder explicitly (instead of one derivatives root scanned recursively) keeps renamed/versioned tool-6 runs apart on the participant dropdown. `find_participants(raw_root)` fills the dropdown; the params JSON + custom stages are read from the fif's folder / the data folder; a summary reports epoch/rejection counts, channels, sfreq, and threshold source. **Output derivation** (data vs reports split, precomputed at load): clean `.fif` → `raw_root.parent/clean_epo_manual/<subtree>/` (beside the raw folder), reports → `<data_root>/reports_rejection_manual/<subtree>/` (beside `reports_preprocessing/`), where `<subtree>` is the participant path relative to the raw folder and `data_root` is the selected Data folder or, if unset, derived by finding the `derivatives/` ancestor of the raw folder (its parent).

**Section 2 — Global per-stage report** (Run): per sleep stage (`W/N1/N2/N3/R` + custom) — a **PSD overlay** (per-epoch mean-across-channel PSD: clean epochs grey + median/IQR band, rejected epochs coloured by their reject method), **metric distributions** (p-p, gradient, 1/f MAE, 1/f R² — clean vs rejected, with threshold lines), a **p-p vs gradient scatter** (method-coloured), and a **stage × method rejection table** (from the metadata flags). 1/f is fitted per epoch (~1 min for a full night); the editable thresholds drive only the reference lines. Assembled into an `mne.Report`; **Save** writes `{file_id}_qc2b_report.html` into `reports_rejection_manual/` (see Section 4).

**Section 3 — Per-epoch navigator** (Run): walks the rejected epochs (filterable by method or stage). For the current epoch: a **stacked montage** (± context) with the current-epoch trace of each channel coloured by its recomputed flagging method, the gradient-max sample marked, and method-coloured channel labels. An optional **"Show EOG/EMG context"** checkbox (default off) stacks the EOG-L/EOG-R/EMG traces (own scale, dotted divider) below the EEG montage when the `{file_id}_context-epo.fif` companion exists — loaded on demand via `load_context_epochs` and aligned by epoch index (no raw-EDF reload); absent companion → the toggle is a no-op. Plus a **detail panel** — PSD + aperiodic fit (worst-R² channel), mean band power (δ/θ/α/σ/β) + 50 Hz ratio, an epoch spectrogram, and a per-channel metric table (value vs threshold). A **keep / reject** toggle overrides the decision in both directions — confirm a rejection or *rescue* a clean-looking flagged epoch.

**Section 4 — Manual override & save**: a review strip (hypnogram + final keep/reject per epoch, overridden epochs marked) and counts (kept / rejected / rescued / newly-rejected). **Save** follows the toolkit's **`derivatives/` (data) vs `reports_*` (reports) split** — both output folders are precomputed once at participant load (`S['out_folder']` / `S['reports_folder']`, reused by the Section-2 and Section-4 save handlers):
- **DATA** → **`derivatives/clean_epo_manual/<edf_subtree>/`** (kept separate from tool 7bis's `clean_epo_auto/`; mirrors the tool-6 subtree, stripping a leading `raw_epo/`): only `{file_id}_clean-epo.fif` (kept epochs only; metadata carries `manual_override` + `final_reject`).
- **REPORTS** → **`reports_rejection_manual/<edf_subtree>/`** (beside `reports_preprocessing/`, sibling of 7bis's `reports_rejection_auto/`; = `deriv_root.parent/reports_rejection_manual/`): `{file_id}_epoch_rejection_reviewed.tsv` (per-epoch metadata + the two override columns), `{file_id}_qc2b_review_log.tsv` (one row per overridden epoch: `epoch_idx`, `stage`, `orig_reject`, `final_reject`, `action` ∈ rescued/added), and the Section-2 `{file_id}_qc2b_report.html`.

Override granularity is **whole-epoch** (only 3–4 EEG channels, and `clean-epo.fif` drops whole epochs anyway).

**Batch twin (`7_reject_manually_batch.py`)**: `python 7_reject_manually_batch.py <derivatives_root> [--no-1f] [--limit N] [--out DIR]`. Runs Section 2 for every participant (writing each `{file_id}_qc2b_report.html`) plus a **database-level aggregate** — rejection rate per participant, rejection rate by stage across participants, pooled metric distributions + scatter — in `qc2b_database_report.html`, alongside `qc2b_database_rejection_summary.tsv` (one row per `file_id × stage`: `n_total`, `n_rejected`, per-method counts, `pct_rejected`). `--no-1f` skips the slow 1/f fitting. The batch is report-only (no `.fif`/override written). Outputs default to `<derivatives_root>/qc2b_reports/`.

**Naming scheme — the three post-flagging tools**: `6_preprocessing` (**flagging** only — no rejection) →
`7_reject_manually` (**manual** rejection, this section) → `7bis_reject_automatically` (**automatic**
rejection by thresholds, next section). Tools 7 and 7bis are **alternative** routes to a `_clean-epo.fif`,
written to sibling folders (`clean_epo_manual/` vs `clean_epo_auto/`) so both can coexist.

### 7bis. Automatic epoch rejection — channel-first then epoch (`7bis_reject_automatically_voila.ipynb`)

Automatically cleans the epochs `6_preprocessing` flagged, with a **channel-first then epoch** decision (à la
PREP / FASTER): a globally-bad channel is **dropped before** the epoch vote, so it no longer condemns every
epoch it appears in. Voila-only (no batch twin — the run loop already processes the whole database). It is a
pure **decision + write** tool: it reads tool-6 outputs **read-only** and never modifies them.

**Inputs — Section 1 selects two folders explicitly, then a `Scan` button** (deliberately *not* one root
scanned recursively — tool 6 is often re-run with different parameters and its outputs renamed/versioned side
by side, e.g. `raw_epo`, `raw_epo_v2`, `reports_preprocessing_v2`; a recursive `rglob` on `derivatives/`
would mix every version, duplicating `file_id`s):
- **Data folder** (`fc_data`, **optional** convenience, first chooser — holds `derivatives/` +
  `reports_preprocessing/` as in tool 6): on selection it `reset(path=…)`s `fc_raw` to `<data>/derivatives`
  (where the `raw_epo*` folders sit) and `fc_reports` to `<data>` (parent of the `reports_preprocessing*`
  folders), so the two choosers below open at the right place. The user still picks the specific
  (possibly versioned) subfolder; it does not select anything on its own.
- **Raw-epochs folder** (`fc_raw`) — holds the tool-6 `{file_id}_all-epo.fif` (e.g. `derivatives/raw_epo`);
  participants are `qc_rejected_epochs_lib.find_participants(raw_root)` and each participant's subtree is
  `fif.parent.relative_to(raw_root)`.
- **Reports folder** (`fc_reports`) — holds `{file_id}_epoch_channel_rejection.tsv` (e.g.
  `reports_preprocessing`); each TSV is located by `reports_root.rglob(...)` **within this folder only**
  (cached at scan time as `S['tsv_by_fid']`, reused by the run). `reject_any` = all methods incl. `event`.
- **Scan** lists the participants, warns (orange) when some have no matching TSV in the reports folder
  (version mismatch), and **derives the stage set from the TSVs' `stage` column** (self-sufficient; AASM +
  any custom stage — no `custom_stages.json` needed), populating one checkbox per stage.
- **Output root** = `raw_root.parent / clean_epo_auto/` — i.e. `clean_epo_auto/` is written **beside the
  chosen raw-epochs folder** (so it lands next to `raw_epo/` under `derivatives/`, and versioned raw folders
  keep their outputs separate). No `resolve_tool6_roots` / `raw_epo`-stripping is needed anymore: the subtree
  is relative to the explicitly-chosen raw folder.

**Decision** (`auto_reject_decision`, both thresholds editable, defaults **20 %** via
`DEFAULT_CHANNEL_REJECT_PCT` / `DEFAULT_EPOCH_REJECT_PCT`; computed over the **stages of interest** only):
1. **Channel-first** — `badness[c] = fraction of in-scope epochs where (e, c) is flagged`; a channel with
   `badness > channel_pct` is **dropped**.
2. **Epoch** — among the **remaining good** channels, an in-scope epoch with
   `fraction of good channels flagged > epoch_pct` is **rejected**.

No interpolation — bad channels are dropped (interpolation needs an electrode montage tool 6 does not set;
deferred). Uses `reject_any` (event, being epoch-level, raises every channel's badness equally — acceptable).

**Edge cases** (non-fatal, add to `failed`, no clean-epo written): none of the selected stages present; **all
channels** would be dropped; **all in-scope epochs** rejected. (On sparse 3-channel montages with high
flagging the 20 % channel threshold can drop every channel — expected; the thresholds are editable.)

**Outputs — the toolkit's `derivatives/` (data) vs `reports_*` (reports) split** (`<subtree>` = the
participant path relative to the chosen raw folder):
- **DATA** — `<raw_root.parent>/clean_epo_auto/<subtree>/{file_id}_clean-epo.fif`: `clean_epo_auto/` beside
  the chosen raw-epochs folder, so when `raw_root == derivatives/raw_epo` it lands at
  `derivatives/clean_epo_auto/` (a sibling of tool 6's `raw_epo/` and tool 7's `clean_epo_manual/`). The
  **only** file here is the `.fif` (selected-stage kept epochs, dropped channels removed; `metadata` gains
  `auto_reject_channel_pct` / `auto_reject_epoch_pct` / `auto_reject_stages` provenance).
- **REPORTS** — `<reports_root.parent>/reports_rejection_auto/<subtree>/`: `reports_rejection_auto/` beside
  the chosen reports folder (i.e. next to `reports_preprocessing/`). Per participant:
  - `{file_id}_autoreject_decision.tsv` — one-row durable record **and** the global-summary row source
    (thresholds, stages, `rejected_channels`, per-channel `channel_badness_pct`, `n_epochs` /
    `n_epochs_rejected` / `n_epochs_kept` / `pct_epochs_rejected`). Written **before** the report.
  - `{file_id}_autoreject_report.html` — an `mne.Report`: a channels × epochs flagged-pair heatmap (dropped
    channels' flagged cells greyed + a red bold label; a rejected-epoch strip; a hypnogram strip — the matrix
    height scales with the channel count while the two strips keep a fixed height) + the per-stage table.

**Global** (at the `reports_rejection_auto/` root): `global_autoreject_summary.tsv` (full schema, **rebuilt
each run by globbing the per-file `_autoreject_decision.tsv` from disk** — interruption-safe), an end-of-run
**compact summary table** shown in the notebook (`file_id`, `n_channels`, `n_channels_rejected` and
`n_epochs_rejected` each with the % in parentheses, `n_epochs`) + a `% epochs rejected per participant` bar
and the **full** table in `autoreject_database_report.html`, and `autoreject_failed.tsv`. **Skip +
cumulative-merge**: a participant is skipped when **both** its clean-epo (in `clean_epo_auto/`) and its
decision TSV (in `reports_rejection_auto/`) exist (uncheck *Skip* to reprocess).

**Curry twin (verbatim copy)**: 7bis reads only format-agnostic MNE `.fif` + the TSV, both produced
identically by the Curry tool 6, so there is **no EDF-specific code to swap** (unlike tools 5/6). The Curry
twin `tools_curry/7bis_reject_automatically_curry_voila.ipynb` is therefore a **verbatim copy** of the EDF
notebook with only the title retitled, produced by `tools_curry/_make_tool7bis_curry.py` (re-run it after
editing the EDF notebook — it re-copies + syntax-checks; the generator asserts every code cell is byte-equal
to the source). The shared-library import search in the setup cell probes `cwd`, `cwd/tools`,
`dirname(cwd)/tools` and `cwd/../tools`, so `tools/qc_rejected_epochs_lib.py` is found whether Voila is
launched from the repo root **or** from `tools_curry/` (where the twin lives, the lib in `../tools`).

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

## Curry 9 (`.cdt`) support — experimental

**Status**: an exploratory port of a subset of the EDF tools to Neuroscan **Curry 9** recordings,
living in a **separate `tools_curry/` folder**. It shares the same `config_param/*.json` conventions
and output layout as the EDF tools (a project is assumed to be *either* EDF *or* Curry, never mixed).
The main line of the toolkit remains EDF-first; this section exists so the adaptation can be **reused
and extended later** without re-deriving the Curry-specific deltas. The EDF tools are untouched.

### What a Curry 9 recording looks like

Each recording is a small directory of sibling files sharing one stem (e.g. `y_S005`):

| File | Role |
|---|---|
| `{stem}.cdt` | binary signal (float32, multi-GB — a full night of 44 ch @ 1024 Hz ≈ 5 GB) |
| `{stem}.cdt.dpo` | **plain-text** parameter sidecar (~7 kB): channel labels, sampling rate, start datetime, per-channel impedances, 3-D sensor positions |
| `{stem}.cdt.ceo` | native events (binary; **not** used — see events below) |
| `{stem}_Hypnogram_Export.txt` | one stage label per 30 s epoch (same format as the EDF hypnograms) |
| `{stem}_ScoredEvents_Export.txt` | scored events, comma-separated text export (often French labels). **Encoding varies between exports** — UTF-16 (with BOM) on some, plain UTF-8/ANSI on others: `load_events_curry()` sniffs the BOM (`\xff\xfe`/`\xfe\xff` → UTF-16, else UTF-8 with `errors="replace"`). Never assume UTF-16 |

Unlike Compumedics EDF, Curry uses a **single global sampling rate** for all channels, stores signal as
**float** (no digital→physical scaling, so no EDF physical bounds), never appends MNE `-0`/`-1` duplicate
suffixes, and carries **real electrode positions**.

### Shared Curry modules (`tools_curry/`, the analogue of the EDF header parser)

- **`curry_header.py` — header-only `.cdt.dpo` parser** (the Curry counterpart of the custom EDF header
  parser; **never** reads the `.cdt` signal). `read_curry_header(cdt_path)` returns a dict with `sfreq`,
  `n_samples`, `n_epochs_30s`, `start_datetime`, `data_unit`, EEG vs "other" channel groups
  (`eeg_group_size` / `ch_labels` from the `LABELS` block vs `LABELS_OTHERS`), `sensor_xyz` (3-D positions,
  mm), and raw impedances. `get_impedance_summary(hdr)` returns per-channel kΩ values — **the file stores
  impedances in Ohms, divide by 1000**; sentinels `-1` = not measured, `-2` = disabled. Motivation:
  `mne.io.read_raw_curry(preload=False)` is **not** header-only (it takes ~15 s because `curryreader` loads
  the whole signal) and exposes **no impedances**, so a dedicated parser is required for tools that only
  need metadata (tool 1).
- **`curry_io.py` — signal / hypnogram / events loading**:
  - `read_curry_signal(cdt_path, include=None, preload=False)` — thin wrapper over `mne.io.read_raw_curry`;
    channel selection is a plain `raw.pick(include)` (no include-at-read trick needed — single global rate).
  - `load_hypnogram_curry(txt_path)` — identical format to the EDF hypnograms.
  - `load_events_curry(txt_path, rec_start_dt)` — parses the UTF-16 export into a `Name/Start/Duration`
    (seconds) DataFrame. The export gives **clock times** (`HH:MM:SS`), converted to seconds-from-start via
    the header start datetime, with **midnight rollover** (an event > 1 h *before* the recording start is on
    the next calendar day). Two parsing gotchas handled: single-digit hours after midnight
    (`0:05:11` — `datetime.time.fromisoformat` rejects these, so split manually) and `M:SS[.s]` durations.
  - `rec_start_from_header(hdr)` — builds the recording-start `datetime` used by `load_events_curry`.

### Recipe: adapting an EDF tool to Curry

Each Curry tool is generated from its EDF twin by a small, re-runnable `tools_curry/_make_toolN_curry.py`
script. `_make_tool{2,4,5,6}_curry.py` do **parsed-JSON edits** (join the cell source, string-replace,
re-split; validated by `json.load` + `ast.parse` of every code cell), so their patterns are ordinary Python
source with real newlines and quotes.

> **`_make_tool3_curry.py` is the exception: it string-replaces the RAW notebook JSON.** In raw JSON a
> newline is the two characters `\n` and a quote is `\"`, so **every pattern in it must be a plain,
> escape-free substring** — a pattern containing a real newline or a bare `"` silently never matches and is
> reported only as a `⚠ NOT FOUND` line. Always check the run prints **5/5 replacements applied**: a lower
> count means a replacement was skipped and EDF wording leaked into the twin. (Such a miss was found in
> Jul 2026 — the "No EDF files" → "No .cdt files" pattern had never matched, and the twin on disk carried a
> hand-edit that regeneration reverted.)

The recurring deltas:

1. **Discovery**: `rglob('*.edf')` → `rglob('*.cdt')` (match bare `f.suffix == '.cdt'`, exclude `._*`).
2. **Header reads** → `read_curry_header()`; classify channels by the `.dpo` **EEG group** (authoritative)
   instead of EDF transducer type / name regex.
3. **Signal loading**: replace the EDF `read_raw_edf(preload=False, include=list(remap.keys()))` +
   `drop_suffix_duplicates` + `adapt_remap_dict_to_suffixes` pattern with
   `read_raw_curry(preload=False)` → `raw.pick(montage)` → `raw.rename_channels(remap)` →
   drop deselected → `load_data()`. The two MNE suffix helpers are **removed** (no `-0`/`-1` in Curry).
4. **Drop EDF-only concepts**: no physical bounds → remove `get_phys_bounds_uV`, the `bounds_pct` metric,
   its threshold widget, and its report/summary columns (tool 5). No `patient_id`/`recording_id` header
   fields → **no anonymizer** (tool 1bis has no Curry twin for now).
5. **Events**: the EDF tool 4 now reads the same `*_ScoredEvents_Export.txt` format (among TXT/CSV/XML);
   the Curry twin swaps that three-source `load_events()` for a wrapper over `load_events_curry`
   (single source, recording start from the `.cdt.dpo` header instead of the EDF header), keeping the
   **same call signature** so downstream scan/rejection code is untouched. Because there is only one event
   source, the 1bis consistency check is **removed**, and so are the EDF-only helpers it replaced
   (`read_edf_start_datetime`, `load_events_from_txt/csv/xml`, `_canon_events`) and the extra `TXT suffix`
   widget (Curry keeps a single suffix field). Suffix auto-detection scans `.txt` files whose name contains
   `event`. The **French label mapping** (`FRENCH_EVENT_RULES` / `suggest_canonical`) is format-agnostic and
   passes through **unchanged**, so the twin gets it for free.
6. **User-facing strings** (English per project rule) and re-runnable generators only; all format-agnostic
   logic — skip/cumulative-merge, custom stages, `plot_hypnospectrogram`, rejection methods + heatmap +
   `METHOD_ORDER`, the tool-6 `_preprocessing_params.json` sidecar read by tool 7 — is copied **unchanged**.

### Tools ported (per-tool deltas beyond the recipe)

- **1 — Inspect** (`1_inspect_curry_voila.ipynb`, purpose-built, not generated): per-channel summary from
  `read_curry_header` (no signal read); adds an **impedance check** with an editable **kΩ threshold**
  (default 20 kΩ). Outputs to `summary_inspection_curry/`: `FULL_summary_table_curry.tsv`,
  `impedance_flags_curry.tsv`, `failed_cdt_read.tsv`.
- **2 — Channel selection & remap**: channel mask is `channel_group == 'EEG'` from the `.dpo`; Section 6
  loads via `read_raw_curry`. Same `remap_reref_persubject.json` output.
- **3 — Hypnogram remap**: discovery only (`.edf`→`.cdt`); the `_Hypnogram_Export.txt` companion and remap
  logic are unchanged.
- **4 — Event harmonization**: events from `*_ScoredEvents_Export.txt`; §1bis removed; same
  `config_param/event_remap.json` output (consumed by tool 6). The `FRENCH_EVENT_RULES` suggestion layer
  (shared with the EDF tool) now pre-fills the common French raw labels (e.g. `Micro-éveil 1 ARO SPONT`
  → `arousal_spontaneous`, `Hypopnée obstructive` → `hypopnea`, `Ronflement` → `snore`); unusual labels
  still need a manual choice.
- **5 — Quality overview**: drops `bounds_pct` and `get_phys_bounds_uV`; keeps flat/std/n_peaks/percentiles/
  histograms/PSD/time-series/hypnospectrogram. Same `reports_quality_overview/` outputs (minus the bounds
  column). Memory note: EEG channels are picked before `load_data()` to avoid holding the full multi-GB file.
- **6 — Preprocessing + epoch rejection**: Curry load block (recipe step 3); event-based rejection reads the
  scored events via the **shared TXT-first / CSV / XML** `load_events` chain (the twin swaps only
  `read_edf_start_datetime` for a `.cdt`-header read; the `_events_df_from_txt/csv/xml` parsers + dispatcher
  and the two suffix widgets pass through, so a Curry `_event_xml.csv` is picked up too). All rejection
  methods, per-stage summaries, `_all-epo.fif`,
  `_epoch_rejection.tsv`, `_rejection_summary.tsv`, and the `_preprocessing_params.json` sidecar are
  identical to tool 6, so **tool 7 (QC of rejected epochs) works on Curry outputs unchanged**.
  - **Per-channel rejection (memory) — shared EDF + Curry**: `compute_rejection_masks` receives the
    `epochs` object and reads the signal **one channel at a time** for the time-domain methods
    (`epochs.get_data(picks=[ci])[:, 0, :] * 1e6`) instead of materialising a full `epochs.get_data() * 1e6`
    array (which would be a 3rd full copy after `raw._data` and `epochs._data` — ~3× ~9.5 GiB for 32 ch @
    1024 Hz, enough to overflow RAM); with `del raw` right after epoching this drops the sustained
    rejection-time peak from ~3× to ~1× (~9.8 GiB). **Formulas are unchanged**, so the per-(epoch, channel)
    masks — and therefore every output (`.fif` metadata, `_epoch_rejection.tsv`, `_rejection_summary.tsv`,
    heatmap) and tool-7 compatibility — are **byte-identical** to the former full-array version (verified
    end-to-end: the sha256 of each output file matches on the test EDFs, old vs new). The transient epoching
    peak stays ~2× (`raw._data` + `epochs._data` coexist during the copy); for extreme density without enough
    RAM the built-in **resample** is the lever (÷4). This now lives in the **EDF source** (`tools/6_preprocessing_voila.ipynb`)
    and is format-agnostic, so it **passes through into the Curry twin unchanged** — the generator no longer
    injects it.
  - **Progress feedback — shared EDF + Curry**: the per-participant **pipeline bar** (`progress_step` + the
    4-segment `build_step_legend`, see tool 6 above) and the `_rej_progress` callback that animates its
    *Reject* segment (per channel for the time-domain pass, per epoch for the 1/f fit, filling
    `[_reject_base, _reject_base + COST_REJECT]`) both live in the **EDF source** and are inherited by the
    twin unchanged, so EDF and Curry share the exact same two-bar layout and animation — a slow high-density
    participant is never mistaken for a crash in either. The generator adds **no** separate rejection sub-bar.

### Environment

Curry support requires **MNE ≥ 1.12** (`environment.yml` pins `mne=1.12.*`) plus **`curryreader`** (pip-only,
under the `pip:` block). These were bumped for Curry; the EDF tools continue to run on the same environment.

### Not ported / known limitations

- **Tool 1bis (anonymizer)** — no equivalent: the `.cdt.dpo` has no patient-name field to scrub (revisit if
  a PII field is found in other exports).
- **Tools 7, 8, 9** — not ported. Tool 7 already consumes Curry tool-6 outputs as-is (no `.cdt` reload).
- **Electrode positions unused** — `curry_header` parses `sensor_xyz` and `read_raw_curry` loads a montage,
  but no tool uses them yet; a genuine Curry advantage over position-less Compumedics EDF (enables topomaps,
  geometry-based bad-channel interpolation, spatial re-referencing) — future enhancement.
- **French event vocabulary** — see tool 4 caveat above.

## Planned modules (in development)

Modules still in development. The quality-overview (5), preprocessing (6) and QC-of-rejected-epochs (7) tools above are stable and were promoted out of this section.

### Event-based epoch rejection (Phase A) — **implemented in tool 6**

Promoted out of this section: the **"Event overlap"** rejection method is now part of `6_preprocessing_voila.ipynb` (Section 2 → *Epoch rejection*). See the tool-6 description above for the UI, the `event_remap.json` chooser + auto-detected `Event CSV suffix`, the "Count affected epochs" button, and the additive `flag_event` outputs. (The mirror in `8_live_explore_1file`'s quick rejection remains a possible follow-up.)

### Event-epoch visualizer (Phase B)

Helps decide whether a given event *type* is worth feeding into Phase A. Extends `8_live_explore_1file`'s Section 3 (which already overlays scored-event spans via the shared CSV-first / XML-fallback `load_events()`):
- an **event-type filter** on the navigator: "jump to next/previous epoch overlapping event type X", with a per-type count;
- a per-type **accept/reject decision** widget that writes a small `event_type_decisions.tsv` feeding Phase A's default selection.

### Interactive QC of rejected epochs (Phase 2b) — **implemented as tool 7**

Promoted out of this section: implemented as `7_reject_manually_voila.ipynb` (+ `7_reject_manually_batch.py`, `qc_rejected_epochs_lib.py`) — see the **tool-7 section above**. It loads `{file_id}_all-epo.fif` (+ optional `{file_id}_preprocessing_params.json`), lets the user navigate flagged epochs and override the whole-epoch keep/reject decision, and saves a validated `{file_id}_clean-epo.fif`. **Done**: optional EOG/EMG context in the per-epoch view — declared once in tool 2 (`context_channels`), persisted by tool 6 as `{file_id}_context-epo.fif`, and shown behind the "Show EOG/EMG context" toggle (see the tool-6/tool-7 sections above). **Possible follow-ups**: the code-visible Jupyter twin; per-(epoch, channel) override if a future tool-6 run persists the channel-level mask.

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
