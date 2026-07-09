# Dream-Toolkit — Inspect_EDF

For the full project description, tool inventory, design rationale, byte offsets, verification notes, and
planned modules, see [SPEC.md](SPEC.md). **Read SPEC.md once at the start of each task** — most of the
decisions below are short reminders that point to its *Cross-cutting procedures* section for the detail.

## Key design decisions

The "what to do / what not to break" reminders. Rationale, exact offsets/thresholds, and per-tool
specifics live in SPEC.md.

- **EDF headers: custom binary parser, header-only**: read EDF headers with the hand-written binary parser
  (robust to encoding/header edge cases), **never** for signal data. Derive
  `sampling_frequency = samples_per_record / duration_data_record` (the 8-byte per-channel field is samples
  *per record*, not the rate; they coincide only at 1 s records). Keep in sync across `1_inspect_edf*`,
  `2_select&remap_channels_edf*`, `8_live_explore_1file*`.
- **Signal loading = MNE; export = edfio**: load signal with `mne.io.read_raw_edf()` (never pyedflib for
  signal). Write EDF with `edfio.EdfSignal()`/`edfio.Edf()` directly (not `mne.export.export_raw()`) for
  per-channel physical-range control.
- **MNE duplicate channel suffixes**: MNE ≥ 1.8 appends `-0`/`-1` to duplicate channel names. After every
  `read_raw_edf` use `drop_suffix_duplicates()`, plus `adapt_remap_dict_to_suffixes()` when a remap dict is
  involved.
- **Physical bounds (`get_phys_bounds_uV`)**: reconstruct EDF physical bounds from MNE internals **scaled by
  the channel's native unit** (`units[ch]*1e6`) — otherwise mV channels (Compumedics EOG/EMG/ECG) falsely
  flag ~100 % `bounds_pct`. Duplicated in `5_quality_overview_voila` + `8_live_explore_1file*`; keep in sync.
- **EOG/EMG/ECG detection (`detect_channel_types`)**: classify non-EEG channels by transducer type OR name
  (incl. `chin|menton` for EMG). Reuse the helper from `8_live_explore_1file`, pre-filled as an editable
  selection so the user can correct misses.
- **In-place header anonymization (`1bis_anonymize_edf*`)**: copy the file, then overwrite **only**
  `patient_id` (→ `X X 30-DEC-1899 X_X`) and `recording_id` (keep the real `Startdate` token, blank
  admin/tech/equipment) in the 256-byte header. Everything from byte 256 on stays byte-identical (verified
  by `sha256(file[256:])`); originals are **never** modified; companions are copied + renamed, content
  unchanged.
- **Scored events: CSV-first / XML-fallback (`load_events`)**: read the Compumedics `*_event_xml.csv`
  (`Name, Start, Duration`) first, then the `<ScoredEvents>` of `*.edf.XML`. The CSV suffix is
  configurable/auto-detected (tools 4, 6); labels are harmonized to canonical via
  `config_param/event_remap.json`. Tool 4 returns tuples; tool 7's overlay returns a DataFrame.
- **Event-based epoch rejection (`6_preprocessing_voila`)**: optional 6th method, **onset-only** containment
  (annotated `Duration` ignored). `METHOD_ORDER` is the single source of truth so `event` outputs are purely
  additive (event-free runs stay byte-compatible); event loading is non-fatal per file.
- **Custom (non-AASM) sleep stages**: declared once in `config_param/custom_stages.json` (written only by
  `3_remap_hypno`, read by tools 5/6/7). Three helpers duplicated like `get_phys_bounds_uV`:
  `load_custom_stages`, `parse_custom_field`, `custom_stage_style`. Tools 5/7 use a custom
  `plot_hypnospectrogram()` because YASA's plotting hard-rejects non-AASM labels. Reading the JSON is
  non-fatal.
- **Time-series display cap for DC-coupled data (±500 µV physiological ceiling)**: DC-coupled recordings
  (Curry `.cdt`, and any acquisition exported in DC with no clipping) carry no export clipping, so slow drift
  or artifacts can push the p99.9-based autoscale far past physiological range and crush the real EEG in a
  time-series/butterfly plot. For such tools, **cap** the shared amplitude limit at a wide physiological
  ceiling — `y_lim = min(max_p999, 500.0)` (constant `DISPLAY_YLIM_UV = 500.0`) — never a hard fixed window,
  so clean low-amplitude channels still auto-zoom below the cap and the shared cross-channel scale is kept.
  Applied to the per-channel + butterfly time series **and** the histogram X-axis (`x_lim_hist` follows
  `y_lim_ts`, so capping the one variable covers all three). Currently **Curry-only** (injected by
  `tools_curry/_make_tool5_curry.py`, block "cap time-series y-limit…" — re-run the generator after editing);
  the EDF tools keep the uncapped autoscale on purpose (full range helps spot export clipping). Extend the
  same cap to any future DC-source tool.
- **Flat/dead-epoch colour scaling (`plot_hypnospectrogram()`)**: exclude near-zero (dead-epoch) columns from
  the `vmin/vmax` percentiles and render them grey, else the spectrogram washes out once >2.5 % of epochs are
  fully flat. Kept in sync across tools 5, 8, 8-voila; clean channels stay byte-identical.
- **Dual delivery + outputs**: every user-facing tool ships a code-visible Jupyter notebook **and** a
  code-hidden Voila app (kept in sync); some add a batch `.py`. Outputs are TSV (machine-readable) + HTML
  (human-readable).
- **Tool numbering**: `7` = QC of rejected epochs (Phase 2b); `8` = live single-file explorer
  (`8_live_explore_1file*`, was tool 7); `9` = spectral analysis (`9_SpectralPower*`, was tool 8).
- **QC of rejected epochs (tool 7, `7_inspect_rejected_epochs*`)**: reads tool-6 `{file_id}_all-epo.fif`
  (+ optional `{file_id}_preprocessing_params.json`) and **never reloads the raw EDF** (EEG-only, the `.fif`
  channels as-is) and **never modifies tool-6 outputs**. The per-epoch reject decision is authoritative from
  `epochs.metadata`; the per-**channel** attribution shown (which channel drove a flag, margins) is
  **recomputed** from the signal with tool-6's formulas + the persisted thresholds (fallback: tool-6
  defaults). Analysis + plotting live in a **shared module `qc_rejected_epochs_lib.py`** imported by both the
  Voila notebook and the batch `.py` — a deliberate exception to the "duplicate helpers" rule, justified by
  the heavy report code shared across the notebook/batch pair; the copied bits (`METHOD_ORDER`, heatmap
  palette, custom-stage helpers, Welch-PSD config, specparam 1/f fit) must stay in sync with `6_preprocessing`.
  Override is **whole-epoch** keep/reject (only 3–4 EEG channels). Voila-only for now (Jupyter twin planned).
- **Tool-6 threshold sidecar**: `6_preprocessing_voila` writes `{file_id}_preprocessing_params.json` next to
  the `.fif` (resample/filter + per-stage rejection thresholds actually used + `1f_fit_range_hz` +
  `methods_run`); additive and non-fatal, read back by tool 7. The real 1/f fit uses a **full peak model**
  (`peak_width_limits=[0.5,20]`, `min_peak_height=0.3`), **not** `max_n_peaks=0` (which over-rejects N2/REM)
  — SPEC's table reflects this.
- **Configurable 1/f fit range (`6_preprocessing`)**: the aperiodic fit window is user-editable via two
  widgets, **default 2–45 Hz** (`txt_1f_fmin`/`txt_1f_fmax`); `freq_mask = (psd_freqs >= fit_fmin) &
  (psd_freqs <= fit_fmax)`. The Welch PSD ceiling **follows the fit max** (`fmax_psd = min(fit_fmax_val,
  sf/2-0.5)`, was hardcoded 30 Hz), and the default **bandpass is 0.1–50 Hz** (was 0.1–40) so the fit band is
  preserved. The range is persisted as `rejection_thresholds.1f_fit_range_hz` and **read back by tool 7**
  (`qc_rejected_epochs_lib.load_params` → `info['fit_range']`, threaded into `compute_psds(fmax=)` /
  `fit_1f(fmin=)`; fallback `(2.0, 45.0)`) so its recomputed per-channel attribution matches. Keep the fit
  range wired across tool 6 (EDF + Curry), `qc_rejected_epochs_lib.py`, tool-7 batch + voila.
- **Optional notch filter (`6_preprocessing`, `cb_notch`/`txt_notch_freq`, default OFF, 50 Hz)**: removes
  power-line noise via `raw.notch_filter(freqs=notch_freq_val)` with **MNE's default method** (FIR — do not
  pass `method=`). Applied as pipeline block `[D2]` **after re-referencing, before the bandpass** (matching
  the "notch then band-pass" convention); **fatal** on failure like the bandpass block. Single editable
  frequency (no harmonics/list). Off by default keeps `.fif`/TSV outputs byte-identical. Persisted additively
  to the sidecar as `notch: {applied, freq_hz}` for **provenance only** — **tool 7 needs no change** (the
  `.fif` is already notched and the 1/f fit band 2–45 Hz excludes 50 Hz). Wired to the Curry twin: it passes
  through the generator untouched (no `_make_tool6_curry.py` edit), so **re-run the generator** after editing.
- **Optional resampling in `5_quality_overview` (`cb_resample`/`txt_target_freq`, default OFF, 256 Hz)**:
  mirrors tool 6 — `raw.resample(target_freq, npad='auto')` applied right after load/rename (guarded to never
  upsample, non-fatal). All quality metrics/spectrograms are then computed on the resampled signal. Off by
  default keeps results byte-identical. Placed **outside** the block the Curry generator replaces, so it
  carries through to the Curry twin unchanged.
- **Curry twins are generated, not hand-edited**: `tools_curry/_make_tool{5,6}_curry.py` regenerate the
  Curry notebooks from the EDF originals by string replacement. Edit the EDF notebook, then **re-run the
  generator**. New code passes through automatically **unless** it sits inside a block the generator
  string-matches or wholesale-replaces — e.g. tool 6's generator replaces `compute_rejection_masks` (with the
  memory-efficient **per-channel** read) and the rejection call site, so any signature change there must be
  mirrored in `_make_tool6_curry.py` (`NEW_REJECT` + the `cell8 rejection call` OLD/NEW strings).
  Tool 6 does **not** persist a per-(epoch, channel) mask; per-participant files are `_all-epo.fif`,
  `_epoch_rejection.tsv`, `_rejection_summary.tsv`, `_preprocessing_report.html` (heatmap embedded, no PNG).

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
- **Editing large notebooks (rename-to-`.txt` method)**: `Read` ignores `offset`/`limit` on `.ipynb` and
  fails once total size passes ~25 k tokens; `Edit` is blocked on the `.ipynb` extension and `NotebookEdit`
  needs a prior whole-file `Read`. So rename around the extension blocks:
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

## Process for new tasks

- **Read SPEC.md first**, then **ask clarifying questions** (inputs/outputs, edge cases, UI layout, backward
  compatibility, interaction with existing tools, output file names) before planning — resolve ambiguities
  first.
- **Plan-first**: draft the plan with the most capable available model, write it to a temporary
  `tools/plan_<task>.md` for the user to read/annotate, and start implementing **only** after explicit
  confirmation. Delete the plan file when the task is done and confirmed.
- **Keep docs in sync as part of the task**: add any directly-imported package to `environment.yml`
  (conda-forge if available, else pip); update **CLAUDE.md** directly when a design decision or dev rule
  changes; **propose** SPEC.md updates for any new feature/workflow/constant/output/module — but don't write
  them without explicit confirmation.
