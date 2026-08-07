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
- **Context channels declared once in tool 2 (`2_select&remap_channels_edf*`, Section 2bis)**: identifies,
  **per channel configuration**, the non-EEG *context* channels — `eog_left`, `eog_right`, `emg`, `ecg` —
  via a dataframe-based `detect_context_channels(df_sub)` (same *transducer-type OR name* convention, but on
  the header scan, **not** the raw-based `detect_channel_types`; **must tolerate a missing `transducer_type`
  column** so the generated Curry twin still works — name-only there). Auto-detected into editable dropdowns
  (`(none)` if absent), fanned out to participants like `remap`/`ref_channels`, and saved in an **additive,
  optional** `context_channels` block of `remap_reref_persubject.json` (**omitted when empty** → JSONs
  without context stay byte-identical). Analysis tools 5/6 ignore it; only tools that explicitly ask (e.g.
  tool 7's per-epoch inspector) load these channels on demand. EDF Voila + Jupyter kept in sync; the Curry
  twin is **regenerated** via `tools_curry/_make_tool2_curry.py` (new code passes through — no generator
  edit needed; re-run it after editing the EDF Voila).
- **`get_phys_bounds_uV`**: scale MNE's physical bounds by `units[ch]*1e6` — otherwise mV channels
  (Compumedics EOG/EMG/ECG) falsely flag ~100 % `bounds_pct`. Duplicated in tools 5 & 8 — **keep in sync**.

**Shared conventions** (→ SPEC *Cross-cutting procedures*)
- **Scored events (`load_events`)**: tool 4 reads three Compumedics companions **TXT-first** — the
  `*_ScoredEvents_Export.txt` text export (French, UTF-16-or-UTF-8, no header, clock-time; parsed inline like
  `curry_io._parse_events_txt`), then `*_event_xml.csv`, then the `<ScoredEvents>` of `*.edf.XML`. The scan
  needs **only names** (Start/Duration left NaN, no EDF-header read); `read_edf_start_datetime` (header
  offsets 168/176) converts `.txt` clock times to seconds **only** in the 1bis check, which compares the
  primary text source vs XML **language-robustly**: `_canon_events` normalizes names to the canonical vocab
  (so EN/FR compare equal) + floors starts to the second, then `_match_events` greedily pairs events
  **within ±`tol` s** (editable `Match tol (s)`, default 1, 0=strict; same-label pass then any-label pass)
  classifying them as matched / cooccur-difflabel (paired in time, different name — cross-language
  discovery aid, or two distinct events within tol) / only-in-one-source. French
  labels get `FRENCH_EVENT_RULES`/`suggest_canonical` (shared verbatim with the Curry twin). **Tool 6 joins
  the same TXT-first / CSV / XML `load_events`** but returns a DataFrame and — because it flags epochs by
  event **onset**, not just names — **reads the EDF start datetime to convert `.txt` clock times to real
  `Start` seconds** (two suffix widgets: Event TXT + Event CSV). Its Curry twin mirrors the whole chain,
  swapping only `read_edf_start_datetime` for a `.cdt`-header read — so `_make_tool6_curry.py`'s match
  strings (`OLD_STARTDT`, the cell-4 detection loop vars) must track any edit to tool 6's event section;
  re-run the generator. **Tool 8 is unchanged** (CSV-first/XML-fallback). Re-run
  `tools_curry/_make_tool4_curry.py` after editing the EDF Voila.
  Labels harmonized to canonical via `config_param/event_remap.json`. → SPEC *Cross-cutting → Event sourcing* + §4.
- **Per-event-type persistence (tool 6 → 7bis/7)**: when event flagging runs on a file that has events, tool 6
  writes two **optional/additive** sidecars beside `_epoch_channel_rejection.tsv` — `{id}_event_epoch_flags.tsv`
  (per-epoch `evt_<type>` flags, for 7bis event sub-selection) and `{id}_event_counts.tsv` (raw `n_events` per
  type, for 7bis's per-event participant-exclusion threshold) — **and** adds the same `evt_<type>` columns to
  the `.fif` epoch metadata (so tool 7's inspector names the culprit event). All three derive from
  `event_type_masks` (`build_event_epoch_flags`/`build_event_counts`); **no events → nothing written, outputs
  byte-identical**. Format-agnostic → pass through the Curry generator (re-run `_make_tool6_curry.py`).
- **Rejection-method palette = single source in `qc_rejected_epochs_lib.HEATMAP_COLORS`** (CVD-validated;
  index 0 `none` dark, `multiple` cyan; `CTX_COLOR` deliberately outside the six method hues). Duplicated
  **verbatim** in tool 6's `plot_rejection_heatmap` and in tools 8/8-voila (whose array is shorter, has no
  `event` and uses a different order → remap per method, never copy-paste). 7bis imports the lib. Six
  categories can't be told apart by colour alone → always ship a legend / method-coloured labels. Changing
  it changes rendered PNGs; re-run the Curry generators. → SPEC *Cross-cutting → Rejection-method colour palette*.
- **Event onsets (tool 6 → tool 7)**: tool 6 persists `{file_id}_event_onsets.tsv` (`type`, `onset_s`,
  `duration_s`) so tool 7 can mark event onsets *inside* an epoch without reloading the raw. It is a
  **montage companion** → written **beside the `.fif` in `derivatives/raw_epo/`** (like `_context-epo.fif`),
  **not** in `reports_preprocessing/` (tool 7 never reads the reports tree). Optional/additive; passes
  through the Curry generator. → SPEC *Cross-cutting → Event-onset sidecar* + §6.
- **Per-epoch montage scales are FIXED, never autoscaled** (`DISPLAY_SCALE_UV` = EEG 150 / EOG 300 /
  EMG 100 / ECG 1000 µV per row, clinical conventions): window-based autoscaling made the same waveform
  change height between epochs and destroyed visual amplitude criteria. Overflow into the neighbouring row
  is intended (clinical-viewer behaviour) — never clip or hide amplitude. → SPEC *Cross-cutting → Fixed
  clinical display scales*.
- **Custom (non-AASM) sleep stages**: declared once in `config_param/custom_stages.json` (written only by
  `3_remap_hypno`, read by tools 5/6/7). Three duplicated helpers (`load_custom_stages`,
  `parse_custom_field`, `custom_stage_style`); tools 5/7 use a custom `plot_hypnospectrogram()` because
  YASA's plotting hard-rejects non-AASM labels. Reading the JSON is non-fatal.
- **Flat/dead-epoch colour scaling (`plot_hypnospectrogram()`)**: exclude near-zero (dead-epoch) columns
  from the vmin/vmax percentiles and render them grey, **then cap the colour span
  (`vmin = max(vmin, vmax − 45)`)** so a continuum of partly-flat/clipped low-power epochs can't drag
  `vmin` and wash the plot to red (clean channels span < 45 dB → untouched, still byte-identical). Keep
  both in sync across tools 5, 8, 8-voila (+ tool 5's Curry twin via `_make_tool5_curry.py`). → SPEC
  *Cross-cutting → Flat/dead-epoch colour scaling*.
- **Time-series display cap for DC-coupled data**: DC sources (Curry `.cdt`, any DC export with no
  clipping) → cap the shared amplitude limit at `y_lim = min(max_p999, 500.0)` (`DISPLAY_YLIM_UV = 500.0`),
  never a fixed window. Currently **Curry-only**; EDF tools keep the uncapped autoscale on purpose (full
  range helps spot export clipping). Extend to any future DC-source tool.
- **Distribution `n_peaks` histogram = robust-clipped**: in tool 5 the histogram feeding the Savitzky-Golay
  curve + `find_peaks` is bounded to robust percentiles (`HIST_CLIP_PCT`, **uniform EDF + Curry**) so rare
  extremes don't over-smooth the curve into a dome/flat line on DC data; `hist_extreme_pct` stays on a
  **separate full-range** histogram (only `n_peaks` changes, now more sensitive). Passes through the Curry
  generator unchanged. → SPEC §5 *Distribution histogram — robust range for peak detection*.
- **Checkbox-revealed widgets**: a parameter box shown/hidden by a checkbox must derive its **initial**
  `display` from that checkbox (`display='' if cb.value else 'none'`), never hard-code `'none'` — the
  `observe` handler only fires on a *change*, so a pre-ticked checkbox would leave its box hidden (as it did
  in the Curry twin of tool 5, where `hp_check` defaults ON). Applies to tools 5 & 6; use it for any new one.
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
  50 Hz) now **defaults ON** (a power-line notch is the norm; untick to reproduce the byte-identical no-notch
  output — its `methods_run` entry keeps tool 7 consistent); resampling (also in tool 5, `cb_resample`) and
  the 1/f fit range (default 2–45 Hz) stay optional and **off/neutral by default**. **PSD smoothing before
  the 1/f fit** (`cb_ssmatch` + optional `cb_ss_median`, two-stage median→LOWESS on the linear Welch PSD,
  oscip parity) is built in and **defaults ON** (the one *not* byte-neutral default — untick to reproduce the
  pre-smoothing output; reprocess to homogenise), recorded in the additive `psd_smoothing` sidecar block and
  **honoured by tool 7** when it recomputes PSDs (`smooth_psd_median`/`smooth_psd_lowess` duplicated in tool 6
  + `qc_rejected_epochs_lib.py` — keep in sync; 7bis untouched, it never recomputes a PSD). In the 1/f box the
  **fit parameters** (fit range + smoothing) are grouped one indent shallower than the **thresholds**, and all
  methods' thresholds share that deeper indent so they stay aligned (Variant A; layout-only). `compute_rejection_masks`
  reads the signal **one channel at a
  time** from `epochs` (+ `del raw` after epoching) to bound memory on dense montages; formulas unchanged so
  every output is **byte-identical** to the former full-array version, and the code is format-agnostic (EDF
  source, passes through to the Curry twin). Keep the formulas in sync with tool 7's `qc_rejected_epochs_lib.py`.
  **Configurable epoch length** (`dd_epoch_len`, divisors of 30 ≥ 5 s, **default 30 = classic**): the 30 s
  hypnogram is validated/trimmed at 30 s then **expanded** (`np.repeat(expert_hypno, 30 // epoch_sec)`) so each
  sub-epoch inherits its parent stage; `make_fixed_length_epochs(duration=epoch_sec)`, event mask + context
  companion + Welch window (`min(4, epoch_sec)·sf`) all take `epoch_sec`, sidecar gains `epoch_length_s`
  (30 when absent). Event flagging stays onset-only with an **amber warning** when ≠ 30 s. Thread `epoch_sec`
  consistently across tool 6 → context companion → sidecar → tools 7/7bis/`qc_rejected_epochs_lib` (re-run the
  Curry generator). → SPEC §6 *Epoching (configurable epoch length)*.
  **Output location**: tool 6 writes its `.fif`/sidecar/context under **`derivatives/raw_epo/<subtree>/`**
  (sibling of tool 7's `clean_epo_manual/` and 7bis's `clean_epo_auto/`); tools 7/7bis **strip a leading
  `raw_epo/`** when mirroring the subtree into their `clean_epo_*/` (pre-`raw_epo` layouts pass through), and
  tool 6 deletes stale same-id copies at the old `derivatives/<subtree>/` location on reprocess.
- **Context-channels companion (`6_preprocessing` block `[H]`, → tool 7)**: when the participant's
  `sub_config` carries a tool-2 `context_channels` block, tool 6 reads **only** those declared EOG/EMG/ECG
  channels from the raw EDF, renames them to role labels (`EOG-L`/`EOG-R`/`EMG`/`ECG`), sets MNE channel
  types, applies the same optional resample, **display-filters per role** (EOG band-pass 0.3–35 Hz, EMG
  high-pass 10 Hz, ECG band-pass 0.5–40 Hz — it's a display companion, so filtering in place is fine and
  makes the tool-7 epoch montage readable), epochs them **identically to the EEG** (`make_fixed_length_epochs`
  duration=30 from t=0 — same epoch count regardless of sfreq, so 1:1 index alignment), and saves
  `{file_id}_context-epo.fif`. **Optional + non-fatal**: no `context_channels` → no companion written, all
  other outputs byte-identical. The Curry generator swaps the reader (`read_raw_edf`→`read_raw_curry`) and
  drops the suffix-dedup line via a dedicated `_make_tool6_curry.py` replacement — re-run the generator after
  editing the `[H]` block.
- **Manually reject flagged epochs (`7_reject_manually*`)**: reads tool-6 `{file_id}_all-epo.fif`
  (+ optional params JSON), **never reloads the raw EDF**, **never modifies tool-6 outputs**. **Section 1 =
  optional data folder + Raw-epochs folder** (like 7bis **minus** the reports picker — tool 7 never reads
  `reports_preprocessing/`); `find_participants` runs on the chosen raw folder so versioned tool-6 runs stay
  apart on the participant dropdown. `data_root` = the selected data folder, else the `derivatives/`
  ancestor's parent. An **"already processed" badge** (clean-epo **and** a review/decision TSV on disk) is
  informative only — it blocks nothing.
  **The decision is RECOMPOSED, not read**: Section 1 checkboxes (stages / methods / event types, all on by
  default) feed `recompute_reject`, which rebuilds `base_reject` + `reject_method` the way 7bis's
  `build_pair_matrix` does — all ticked ⇒ **identical to tool 6's `reject_flag`** (keep it that way). It
  drives the navigator, Section 2 **and** the clean-epo (**out-of-scope stages are excluded from the
  `.fif`**), and gates the cost (`fit_mask=in_scope`, `do_1f=False` when no 1/f method is ticked).
  Per-channel attribution is **recomputed** with tool-6 formulas + persisted thresholds. Section 4 also
  writes a **7bis-style decision record** (`_manualreject_decision.tsv` → `_manualreject_report.html` →
  `global_manualreject_summary.tsv` **rebuilt by globbing the per-file TSVs**; data before report). Analysis + plotting live in a **shared module
  `qc_rejected_epochs_lib.py`** — a deliberate exception to the "duplicate helpers" rule; the copied bits
  (`METHOD_ORDER`, palette, custom-stage helpers, Welch-PSD, **PSD smoothing `smooth_psd_median`/`smooth_psd_lowess`**,
  1/f fit) must stay in sync with tool 6. `compute_psds(…, smoothing=info['psd_smoothing'])` re-applies the
  tool-6 smoothing so the recomputed plots/attribution match the flags (`smoothing=None` → byte-identical).
  An **optional "Show EOG/EMG context" toggle** (default off) stacks the EOG-L/EOG-R/EMG traces under the
  per-epoch montage, loaded on demand from the `{file_id}_context-epo.fif` companion (`load_context_epochs`,
  aligned by epoch index) — still no raw-EDF reload; absent companion → toggle is a no-op. The navigator
  header names the **event type(s)** that flagged the epoch, read from the `.fif` `evt_<type>` metadata
  columns (no reports-folder access; present only when tool-6 event flagging ran). **Data/reports
  split** (both folders precomputed at load: `S['out_folder']` / `S['reports_folder']`): `_clean-epo.fif` →
  **`derivatives/clean_epo_manual/<subtree>/`**; reviewed TSVs + `qc2b_report.html` → **`reports_rejection_manual/`**
  (`deriv_root.parent/…`, beside `reports_preprocessing/`).
- **Automatic epoch rejection (`7bis_reject_automatically_voila`)**: **flagging → rejection** decision tool
  (naming scheme: 6 flag → 7 manual → 7bis auto; 7/7bis are alternatives). Voila; it reads only
  format-agnostic `.fif` + `{file_id}_epoch_channel_rejection.tsv` (both identical from Curry tool 6), so its
  **Curry twin is a verbatim copy** — `tools_curry/_make_tool7bis_curry.py` re-copies + retitles (no
  EDF-specific code to swap; re-run after editing the EDF notebook). **Channel-first then epoch**
  (`auto_reject_decision`): drop a channel flagged in
  > `channel_pct` of in-scope epochs, then reject an in-scope epoch flagged in > `epoch_pct` of the *remaining
  good* channels (both default 20 %, editable; computed over the **stages of interest** only). Reads tool-6
  outputs **read-only**. **Section 2 selectable flagging methods + event sub-selection**: `build_pair_matrix`
  recomposes the (epoch×channel) matrix from the ticked `flag_<method>` columns (all on = the stored
  `reject_any`, byte-identical default); ticking `event` reveals per-canonical-type rows (checkbox +
  `exclude if > N events` threshold + a mean/median-per-file hint) fed by the tool-6 `_event_epoch_flags.tsv`
  (per-type epoch flags; absent → fall back to `flag_event`) and `_event_counts.tsv` (raw counts). A ticked
  type over its threshold **excludes the whole participant** (no clean-epo; row/report marked `excluded`).
  Provenance columns `methods_used`/`event_types_used`/`excluded`/`exclude_reason` (+ `.fif`
  `auto_reject_methods`/`_event_types`); run tally gains `excluded`. Still format-agnostic → **Curry twin stays a
  verbatim copy** (re-run `_make_tool7bis_curry.py`). **Section 1 = two explicit folder pickers (raw-epochs + reports) + a Scan button**
  (+ an optional data-folder chooser that pre-points them) — *not* one root scanned recursively, so
  renamed/versioned tool-6 runs (`raw_epo_v2`, `reports_preprocessing_v2`) stay apart. **Data/reports split**:
  `{file_id}_clean-epo.fif` (selected stages, dropped channels removed) → **`clean_epo_auto/`** (`raw_root.parent/…`,
  beside the raw folder); `_autoreject_decision.tsv` (durable record + global-summary source, written **before**
  the report) + `_autoreject_report.html` + the global summary/report → **`reports_rejection_auto/`**
  (`reports_root.parent/…`, beside the reports folder). No interpolation (channels dropped; deferred). Skip gate
  = clean-epo **and** decision TSV present; global summary rebuilt by globbing the per-file decision TSVs.
  All-channels-/all-epochs-rejected → still get a decision row + report (100 % in the summary), no clean-epo. See SPEC §7bis.

**Curry twins** (→ SPEC *Curry 9 support*)
- **Generated, not hand-edited**: `tools_curry/_make_tool{5,6}_curry.py` regenerate the Curry notebooks
  from the EDF originals by string replacement. Edit the EDF notebook, then **re-run the generator**. New
  code passes through automatically **unless** it sits inside a block the generator string-matches or
  wholesale-replaces (e.g. tool 6's Curry load block, event loader, or `[H]` context reader), in which case
  the change must be mirrored in the generator. (`7bis` is the exception — format-agnostic, so
  `_make_tool7bis_curry.py` is a **verbatim copy + retitle**.)
- **Canonical launch cwd = repo root; shared-module imports must be cwd-independent**: every tool is run
  from `Inspect_EDF/` (see SPEC *How to run*), but the shared-module import must still resolve from either
  cwd. The Curry twins (1,2,4,5,6) probe `cwd`, `cwd/tools_curry`, `dirname(cwd)/tools_curry` for
  `curry_header.py`; 7bis (EDF + Curry) probes the `tools/` variants for `qc_rejected_epochs_lib.py`. This
  block lives in the **generators** (except tool 1 = hand-edited, no generator) — fix it there and
  **re-run the generator**, never the old `os.path.dirname(os.path.abspath('__file__'))` (= cwd only,
  breaks from the repo root).

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
  files merged on the item id; aggregated summaries regenerated from all per-item files). **Interruption-safe**:
  skip an item only when **both** its report **and** its durable per-item data are on disk (mismatch → ⚠ warn +
  reprocess); write per-item data **before** the report; rebuild every cumulative/aggregated table by **globbing
  the per-item files on disk**, never from an in-memory `attempted_ids` merge. See *Cross-cutting procedures* in
  SPEC.md; reference impls: tools 5 & 6.
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
