# Dream-Toolkit — Inspect_EDF

For the full project description, tool inventory, design rationale, byte offsets, verification notes, and
planned modules, see [SPEC.md](SPEC.md). **Read SPEC.md once at the start of each task** — the decisions
below are short reminders; the rationale, exact offsets/thresholds, widget names, defaults, and per-tool
specifics live in SPEC.md (each reminder points to the relevant SPEC section).

## Key design decisions

The "what to do / what not to break" reminders, grouped by theme. Each points to SPEC.md for the detail.

**Loading & parsing** (→ SPEC *Cross-cutting procedures* + §1)
- **EDF headers = header-only custom binary parser**: read headers with the hand-written binary parser
  (robust to encoding edge cases), **never** for signal data. Derive `sampling_frequency =
  samples_per_record / duration_data_record` (the 8-byte field is samples *per record*, not the rate).
  Keep in sync across `1_inspect_edf*`, `2_select&remap_channels_edf*`, `8_live_explore_1file*`.
- **Signal = MNE, export = edfio**: load signal with `mne.io.read_raw_edf()` (never pyedflib); write EDF
  with `edfio.EdfSignal()`/`edfio.Edf()` directly (not `mne.export.export_raw()`) for per-channel
  physical-range control.
- **MNE duplicate-channel suffixes**: MNE ≥ 1.8 appends `-0`/`-1` to duplicate names. After every
  `read_raw_edf` call `drop_suffix_duplicates()`, plus `adapt_remap_dict_to_suffixes()` when a remap dict
  is involved.

**Channel typing** (→ SPEC *Cross-cutting procedures* + §1)
- **`detect_channel_types`**: classify non-EEG channels by transducer type OR name (incl. `chin|menton`
  for EMG); reuse the helper from `8_live_explore_1file`, pre-filled as an editable selection so the user
  can correct misses. Tool 1 inspects EEG/EOG/ECG/EMG with the same *transducer-type OR curated
  name-list* convention, **default selection = EEG + EOG** (ECG/EMG opt-in) across all four files.
- **`get_phys_bounds_uV`**: scale MNE's physical bounds by `units[ch]*1e6` — otherwise mV channels
  (Compumedics EOG/EMG/ECG) falsely flag ~100 % `bounds_pct`. Duplicated in tools 5 & 8 — **keep in sync**.

**Shared conventions** (→ SPEC *Cross-cutting procedures*)
- **Scored events (`load_events`)**: read the Compumedics `*_event_xml.csv` first, then the
  `<ScoredEvents>` of `*.edf.XML`; labels harmonized to canonical via `config_param/event_remap.json`.
- **Custom (non-AASM) sleep stages**: declared once in `config_param/custom_stages.json` (written only by
  `3_remap_hypno`, read by tools 5/6/7). Three duplicated helpers (`load_custom_stages`,
  `parse_custom_field`, `custom_stage_style`); tools 5/7 use a custom `plot_hypnospectrogram()` because
  YASA's plotting hard-rejects non-AASM labels. Reading the JSON is non-fatal.
- **Flat/dead-epoch colour scaling (`plot_hypnospectrogram()`)**: exclude near-zero (dead-epoch) columns
  from the vmin/vmax percentiles and render them grey; keep in sync across tools 5, 8, 8-voila.
- **Time-series display cap for DC-coupled data**: DC sources (Curry `.cdt`, any DC export with no
  clipping) → cap the shared amplitude limit at `y_lim = min(max_p999, 500.0)` (`DISPLAY_YLIM_UV = 500.0`),
  never a fixed window. Currently **Curry-only**; EDF tools keep the uncapped autoscale on purpose (full
  range helps spot export clipping). Extend to any future DC-source tool.
- **Dual delivery**: every user-facing tool ships a code-visible Jupyter notebook **and** a code-hidden
  Voila app (kept in sync); some add a batch `.py`. Outputs are TSV (machine) + HTML (human).

**Per-tool** (→ SPEC §1bis / §5 / §6 / §7)
- **Anonymization (`1bis_anonymize_edf*`)**: copy the file, then overwrite **only** `patient_id` and
  `recording_id` in the 256-byte header; everything from byte 256 on stays byte-identical (verified by
  `sha256(file[256:])`); originals are **never** modified.
- **Preprocessing (`6_preprocessing`)**: `METHOD_ORDER` is the single source of truth so optional
  additions (event rejection, notch, resample, configurable 1/f fit range) keep event/feature-free runs
  **byte-compatible**. Writes the `{file_id}_preprocessing_params.json` sidecar (thresholds actually used
  + `1f_fit_range_hz` + `methods_run`), read back by tool 7. Notch (`cb_notch`, MNE default FIR method,
  50 Hz), resampling (also in tool 5, `cb_resample`), and the 1/f fit range (default 2–45 Hz) are all
  optional and **off/neutral by default**.
- **QC of rejected epochs (`7_inspect_rejected_epochs*`)**: reads tool-6 `{file_id}_all-epo.fif`
  (+ optional params JSON), **never reloads the raw EDF**, **never modifies tool-6 outputs**. Per-epoch
  reject decision is authoritative from `epochs.metadata`; per-channel attribution is **recomputed** with
  tool-6 formulas + persisted thresholds. Analysis + plotting live in a **shared module
  `qc_rejected_epochs_lib.py`** — a deliberate exception to the "duplicate helpers" rule; the copied bits
  (`METHOD_ORDER`, palette, custom-stage helpers, Welch-PSD, 1/f fit) must stay in sync with tool 6.

**Curry twins** (→ SPEC *Curry 9 support*)
- **Generated, not hand-edited**: `tools_curry/_make_tool{5,6}_curry.py` regenerate the Curry notebooks
  from the EDF originals by string replacement. Edit the EDF notebook, then **re-run the generator**. New
  code passes through automatically **unless** it sits inside a block the generator string-matches or
  wholesale-replaces (e.g. tool 6's memory-efficient per-channel `compute_rejection_masks`), in which case
  the change must be mirrored in the generator.

## Working agreements

- **Developer background**: sleep research engineer, PSG/EEG expert, limited software-engineering experience.
  Prioritize readability and correctness over abstraction or cleverness.
- **No unnecessary abstraction**: explicit, readable functions — three clear lines beat a clever one-liner.
  No new helpers/classes unless the complexity clearly justifies it.
- **Comments**: only where the EEG/PSG domain logic isn't obvious to a non-sleep-researcher (e.g. why a
  500 µV threshold flags clipping, or why stage 4 is remapped to N3).
- **Proactive error handling**: add `try/except` in every new feature/notebook. Fatal step (file I/O,
  epoching) → add to a `failed` list and `continue`; non-fatal step (re-referencing, report) → `⚠` warning
  and continue. In Voila, wrap every button callback and per-item loop so one failure never crashes the run
  or freezes the UI; always surface errors via a widget or `print()`.
- **Skip + cumulative-merge**: every per-participant processing tool has a "Skip already processed" checkbox
  (on by default), an "N / M already done" info line, and merge/replace output semantics (cumulative per-row
  files merged on the item id; aggregated summaries regenerated from all per-item files). See *Cross-cutting
  procedures* in SPEC.md; reference impls: tools 5 & 6.
- **Normalize path comparisons**: wrap **both** sides of any path/stem/filename string comparison
  (`==`, `in`, `.isin()`, set/dict membership) in `os.path.normcase(...)` **at the comparison only** (keep
  the stored/displayed value original). Prevents skip checks silently failing on `C:`/`c:` and `/`/`\`.
- **Language**: all user-facing strings in notebooks and all of SPEC.md must be in English.
- **No formal test suite**: validate against real EDF files from the dataset; don't add pytest/unittest
  unless explicitly asked.

## Operational notes (this machine)

- **OS**: Windows 10 + Conda, PowerShell default. Use forward slashes in Python path strings (MNE/pathlib
  handle them on Windows).
- **Running Python**: `conda` is **not** on PATH; miniforge3 is at `$env:LOCALAPPDATA\miniforge3`. Always use
  the explicit interpreter — not `conda run`, `conda activate`, or bare `python`:
  ```powershell
  & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" tools/my_script.py
  ```
- **Editing large notebooks**: `Read`/`Edit`/`NotebookEdit` are blocked or size-capped on big `.ipynb`
  files — use the **rename-to-`.txt`** method (and the JSON-surgery fallback) in SPEC.md →
  *Repository structure → Editing the larger notebooks*.

## Process for new tasks

- **Read SPEC.md first**, then **ask clarifying questions** (inputs/outputs, edge cases, UI layout, backward
  compatibility, interaction with existing tools, output file names) before planning — resolve ambiguities
  first.
- **Plan-first**: draft the plan with the most capable available model, write it to a temporary
  `tools/plan_<task>.md` for the user to read/annotate, and start implementing **only** after explicit
  confirmation. Delete the plan file when the task is done and confirmed.
- **Keep docs in sync — details go to SPEC.md, not CLAUDE.md**: add any directly-imported package to
  `environment.yml` (conda-forge if available, else pip). Put the *details* of any design decision (offsets,
  thresholds, widget names, defaults, output files) in **SPEC.md** — the *Cross-cutting procedures* section
  or the relevant tool section — and **propose** those SPEC edits for confirmation. Only update **CLAUDE.md**
  when a *cross-cutting rule* or a *sync constraint* changes, and then only as a short reminder + pointer to
  SPEC (keep CLAUDE.md short).
