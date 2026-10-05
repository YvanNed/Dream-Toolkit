# Dream-Toolkit: Inspect_EDF

For the full project description, tool inventory, design rationale, byte offsets, verification notes, and
planned modules, see [SPEC.md](SPEC.md). **Read SPEC.md once at the start of each task**: the decisions
below are short reminders. The rationale, exact offsets/thresholds, widget names, defaults, and per-tool
specifics live in SPEC.md (each reminder points to the relevant SPEC section).

## Key design decisions

The "what to do / what not to break" reminders, grouped by theme. Each points to SPEC.md for the detail.

**Loading & parsing** (→ SPEC *Cross-cutting procedures* + §1)
- **EDF headers = header-only custom binary parser**: read headers with the hand-written binary parser
  (robust to encoding edge cases), **never** for signal data. Derive `sampling_frequency =
  samples_per_record / duration_data_record` (the 8-byte field is samples *per record*, not the rate).
  Keep in sync across `1_inspect_edf*`, `2_select&remap_channels_edf*`, `0_live_explore_1file*`.
- **Signal = MNE, export = edfio**: load signal with `mne.io.read_raw_edf()` (never pyedflib). Write EDF
  with `edfio.EdfSignal()`/`edfio.Edf()` directly (not `mne.export.export_raw()`) for per-channel
  physical-range control.
- **MNE duplicate-channel suffixes**: MNE ≥ 1.8 appends `-0`/`-1` to duplicate names. After every
  `read_raw_edf` call `drop_suffix_duplicates()`, plus `adapt_remap_dict_to_suffixes()` when a remap dict
  is involved.

**Channel typing** (→ SPEC *Cross-cutting procedures* + §1)
- **`detect_channel_types`**: classify non-EEG channels by transducer type OR name (incl. `chin|menton`
  for EMG). Reuse the helper from `0_live_explore_1file`, pre-filled as an editable selection so the user
  can correct misses. Tool 1 inspects EEG/EOG/ECG/EMG with the same *transducer-type OR curated
  name-list* convention, **default selection = EEG + EOG** (ECG/EMG opt-in) across all four files.
- **Context channels declared once in tool 2 (`2_select&remap_channels_edf*`, Section 2bis)**: identifies,
  **per channel configuration**, the non-EEG *context* channels (`eog_left`, `eog_right`, `emg`, `ecg`)
  via a dataframe-based `detect_context_channels(df_sub)` (same *transducer-type OR name* convention, but on
  the header scan, **not** the raw-based `detect_channel_types`). It **must tolerate a missing `transducer_type`
  column** so the generated Curry twin still works (name-only there). Auto-detected into editable dropdowns
  (`(none)` if absent), fanned out to participants like `remap`/`ref_channels`, and saved in an **additive,
  optional** `context_channels` block of `remap_reref_persubject.json` (**omitted when empty** → JSONs
  without context stay byte-identical). Analysis tools 6/7 ignore it. Only tools that explicitly ask (e.g.
  tool 8's per-epoch inspector) load these channels on demand. EDF Voila + Jupyter kept in sync. The Curry
  twin is **regenerated** via `tools_curry/_make_tool2_curry.py` (new code passes through, no generator
  edit needed). Re-run it after editing the EDF Voila.
- **`get_phys_bounds_uV`**: scale MNE's physical bounds by `units[ch]*1e6`, otherwise mV channels
  (Compumedics EOG/EMG/ECG) falsely flag ~100 % `bounds_pct`. Duplicated in tools 0 & 6: **keep in sync**.

**Shared conventions** (→ SPEC *Cross-cutting procedures*)
- **Rejected data is computed, flagged, left out of the outputs** (→ SPEC *Cross-cutting → Rejected-data
  policy* + *Participant exclusion registry* + *Participant selector*): no tool skips a rejected participant,
  epoch or channel. Tables of tools 6-9 carry `participant_excluded` / `participant_exclude_reason` (**tool 5
  knows nothing of exclusions**: raw extraction, the decisions belong to 5bis), and tool 9's per-epoch
  tables `in_scope` / `rejected` / `channel_dropped`; only group statistics and figures (and per-participant
  figures) use the kept data. Participant exclusions live in **one registry**
  (`config_param/participant_exclusions.tsv`, `tools/participant_selection_lib.py`, shared module): each tool
  rewrites **only its own `source` rows**, a forced inclusion (5bis manual) overrides every exclusion, an
  unreadable registry is never overwritten, identical rows are never rewritten (its date = "a decision
  changed"). A registry change reaches a tool's tables at its next run (Skip ticked = rebuild only); tools
  6/7/8bis/9 warn at scan when the registry is newer than their tables (`registry_newer_than`). Epoch / channel decisions are **tables** (`derivatives/rejection_auto|manual/`
  `_epoch_decision.tsv` + `_channel_decision.tsv`): **no `clean-epo.fif` any more**, tool 9 joins them to
  tool 7's `_all-epo.fif`. Tools 5/6/7/8bis/9 use the shared **`ParticipantSelector`** (all / subset;
  Skip applies to a subset too, with an orange warning box), never one checkbox per participant.
- **Scored events (`load_events`)**: tool 4 reads three Compumedics companions **TXT-first**: the
  `*_ScoredEvents_Export.txt` text export (French, UTF-16-or-UTF-8, no header, clock-time, parsed inline like
  `curry_io._parse_events_txt`), then `*_event_xml.csv`, then the `<ScoredEvents>` of `*.edf.XML`. The scan
  needs **only names** (Start/Duration left NaN, no EDF-header read). `read_edf_start_datetime` (header
  offsets 168/176) converts `.txt` clock times to seconds **only** in the 1bis check, which compares the
  primary text source vs XML **language-robustly**: `_canon_events` normalizes names to the canonical vocab
  (so EN/FR compare equal) + floors starts to the second, then `_match_events` greedily pairs events
  **within ±`tol` s** (editable `Match tol (s)`, default 1, 0=strict, same-label pass then any-label pass)
  classifying them as matched / cooccur-difflabel (paired in time, different name: cross-language
  discovery aid, or two distinct events within tol) / only-in-one-source. French
  labels get `FRENCH_EVENT_RULES`/`suggest_canonical` (shared verbatim with the Curry twin). **Tool 7 joins
  the same TXT-first / CSV / XML `load_events`** but returns a DataFrame and (because it flags epochs by
  event **onset**, not just names) **reads the EDF start datetime to convert `.txt` clock times to real
  `Start` seconds** (two suffix widgets: Event TXT + Event CSV). Its Curry twin mirrors the whole chain,
  swapping only `read_edf_start_datetime` for a `.cdt`-header read, so `_make_tool7_curry.py`'s match
  strings (`OLD_STARTDT`, the cell-4 detection loop vars) must track any edit to tool 7's event section.
  Re-run the generator. **Tool 0 is unchanged** (CSV-first/XML-fallback). Re-run
  `tools_curry/_make_tool4_curry.py` after editing the EDF Voila.
  Labels harmonized to canonical via `config_param/event_remap.json`. → SPEC *Cross-cutting → Event sourcing* + §4.
- **Manual annotations (`{id}_manual_events.tsv`, tool 8)**: events added by hand in tool 8's navigator go
  **beside the recording**, next to the scored companions (an annotation belongs to the *recording*, not to
  one tool-7 run): a **sibling** file, never inside tool 7's `{id}_event_onsets.tsv`. First three columns =
  `_event_onsets.tsv`'s, so the two concatenate. Tool 7 **appends** them only when `cb_manual_events` is
  ticked (**off by default** → byte-identical). Their labels are already canonical, so
  `augment_remap_for_manual` identity-maps **only** the rows tagged `Source='manual'` (an unmapped *scored*
  label must stay unmapped). Tools 0/4 untouched. Re-run both Curry generators. → SPEC *Cross-cutting →
  Event sourcing → Manual annotations* + §8.
- **Per-event-type persistence (tool 7 → 8bis/8)**: when event flagging runs on a file that has events, tool 7
  writes two **optional/additive** sidecars beside `_epoch_channel_rejection.tsv`: `{id}_event_epoch_flags.tsv`
  (per-epoch `evt_<type>` flags, for 8bis event sub-selection) and `{id}_event_counts.tsv` (raw `n_events` per
  type, for 8bis's per-event participant-exclusion threshold). It **also** adds the same `evt_<type>` columns to
  the `.fif` epoch metadata (so tool 8's inspector names the culprit event). All three derive from
  `event_type_masks` (`build_event_epoch_flags`/`build_event_counts`). **No events → nothing written, outputs
  byte-identical**. Format-agnostic → pass through the Curry generator (re-run `_make_tool7_curry.py`).
- **Rejection-method palette = single source in `qc_rejected_epochs_lib.HEATMAP_COLORS`** (CVD-validated,
  index 0 `none` dark, `multiple` cyan, `CTX_COLOR` deliberately outside the six method hues). Duplicated
  **verbatim** in tool 7's `plot_rejection_heatmap` and in tools 0/0-voila (whose array is shorter, has no
  `event` and uses a different order → remap per method, never copy-paste). 8bis imports the lib. Six
  categories can't be told apart by colour alone → always ship a legend / method-coloured labels. Changing
  it changes rendered PNGs. Re-run the Curry generators. → SPEC *Cross-cutting → Rejection-method colour palette*.
- **Event onsets (tool 7 → tool 8)**: tool 7 persists `{file_id}_event_onsets.tsv` (`type`, `onset_s`,
  `duration_s`) so tool 8 can mark event onsets *inside* an epoch without reloading the raw. It is a
  **montage companion** → written **beside the `.fif` in `derivatives/raw_epo/`** (like `_context-epo.fif`),
  **not** in `reports_preprocessing/` (tool 8 never reads the reports tree). Optional/additive. Passes
  through the Curry generator. → SPEC *Cross-cutting → Event-onset sidecar* + §7.
- **High-density montages (32–64 ch) = channel-count-keyed, never format-keyed**: every adaptation lives in
  `qc_rejected_epochs_lib` and is **inert at or below `HD_CHANNEL_THRESHOLD` (12)**, so the 3–6 channel PSG
  path stays byte-identical: keep it that way when touching the montage geometry, the detail panel, the
  channel triage, the flagging heatmap or the 1/f subsample. Two traps: `recompute_reject` must reproduce
  tool 7's `reject_flag` exactly when every channel is kept and the rule is `any` (methods with no
  per-channel column are **broadcast** across kept channels, else the 1/f contribution vanishes on the
  recomputed fallback), and the channel layer must be engaged **only** when the user actually uses it (the
  fallback flags depend on thresholds that differ from the tool-7 run when no params JSON exists).
  → SPEC *Cross-cutting → High-density montage support* + §8.
- **Per-epoch montage scales are FIXED, never autoscaled** (`DISPLAY_SCALE_UV` = EEG 150 / EOG 300 /
  EMG 100 / ECG 1000 µV per row, clinical conventions): window-based autoscaling made the same waveform
  change height between epochs and destroyed visual amplitude criteria. Overflow into the neighbouring row
  is intended (clinical-viewer behaviour): never clip or hide amplitude. → SPEC *Cross-cutting → Fixed
  clinical display scales*.
- **Custom (non-AASM) sleep stages**: declared once in `config_param/custom_stages.json` (written only by
  `3_remap_hypno`, read by tools 6/7/8). Three duplicated helpers (`load_custom_stages`,
  `parse_custom_field`, `custom_stage_style`). Tools 6/8 use a custom `plot_hypnospectrogram()` because
  YASA's plotting hard-rejects non-AASM labels. Reading the JSON is non-fatal.
- **Flat/dead-epoch colour scaling (`plot_hypnospectrogram()`)**: exclude near-zero (dead-epoch) columns
  from the vmin/vmax percentiles and render them grey, **then cap the colour span
  (`vmin = max(vmin, vmax − 45)`)** so a continuum of partly-flat/clipped low-power epochs can't drag
  `vmin` and wash the plot to red (clean channels span < 45 dB → untouched, still byte-identical). Keep
  both in sync across tools 0, 0-voila, 6 (+ tool 6's Curry twin via `_make_tool6_curry.py`). → SPEC
  *Cross-cutting → Flat/dead-epoch colour scaling*.
- **Time-series display cap for DC-coupled data**: DC sources (Curry `.cdt`, any DC export with no
  clipping) → cap the shared amplitude limit at `y_lim = min(max_p999, 500.0)` (`DISPLAY_YLIM_UV = 500.0`),
  never a fixed window. Currently **Curry-only**. EDF tools keep the uncapped autoscale on purpose (full
  range helps spot export clipping). Extend to any future DC-source tool.
- **Distribution `n_peaks` histogram = robust-clipped**: in tool 6 the histogram feeding the Savitzky-Golay
  curve + `find_peaks` is bounded to robust percentiles (`HIST_CLIP_PCT`, **uniform EDF + Curry**) so rare
  extremes don't over-smooth the curve into a dome/flat line on DC data. `hist_extreme_pct` stays on a
  **separate full-range** histogram (only `n_peaks` changes, now more sensitive). Passes through the Curry
  generator unchanged. → SPEC §6 *Distribution histogram (robust range for peak detection)*.
- **Checkbox-revealed widgets**: a parameter box shown/hidden by a checkbox must derive its **initial**
  `display` from that checkbox (`display='' if cb.value else 'none'`), never hard-code `'none'`: the
  `observe` handler only fires on a *change*, so a pre-ticked checkbox would leave its box hidden (as it did
  in the Curry twin of tool 6, where `hp_check` defaults ON). Applies to tools 6 & 7. Use it for any new one.
- **Dual delivery**: every user-facing tool ships a code-visible Jupyter notebook **and** a code-hidden
  Voila app (kept in sync). Some add a batch `.py`. Outputs are TSV (machine) + HTML (human).

**Per-tool** (→ SPEC §1bis / §5 / §7 / §8 / §8bis / §9)
- **Anonymization (`1bis_anonymize_edf*`)**: copy the file, then overwrite **only** `patient_id` and
  `recording_id` in the 256-byte header. Everything from byte 256 on stays byte-identical (verified by
  `sha256(file[256:])`). Originals are **never** modified.
- **Sleep macrostructure (`5_sleep_macrostructure_voila`)**: hypnogram + scored events + **EDF header only**
  (never the epochs, never the EEG signal, single-channel `include=` reads for the optional `Light`/`SpO2`
  only): the documented exception to *feature tools start from the clean epochs*. Binary header read, **not
  edfio** (crashes on 73.edf). Policies live in SPEC §5 (do not re-derive them): 0-based epochs with an
  **exclusive** lights-on, stage durations within SPT, no sleep → NaN (never 0). **MT/custom = YASA
  convention** (in TIB, out of TST/WASO, `SPT = TST + WASO + other`, `OTHER_STAGE_POLICY` commented alt).
  Lights source order txt → `Light` channel → participant table → bounds with `lights_source` provenance and
  the *last-ON-before-first-sleep / first-ON-after-last-sleep* rule. `hypopnea*`/`arousal*` matched by prefix.
  Unmapped labels **not counted** + warned. Copies tool 7's event chain (its XML parser gains an additive
  `Desaturation` column). Skip marker = `{file_id}_sleep_metrics.tsv` (checks TSV written first). Globals
  globbed. Every metric is declared once in the `METRICS` registry: add there, never ad hoc.
- **Participant selection (`5bis_select_participants_voila`)**: reads only tool 5's `global_sleep_metrics.tsv`,
  writes a decision record (`participant_selection.tsv`) beside it **and the registry**: Validate replaces all
  `5bis_criteria` rows, Section 3 saves `5bis_manual` exclusions / forced inclusions at once. Tool 5 stays
  decision-free. `METRICS` / `PROVENANCE_KEYS` / `REFERENCE_RANGES` are **verbatim copies of tool 5's**: edit
  both notebooks together. The unicorn GIF on Validate is deliberate (`SHOW_UNICORN`). → SPEC §5bis.
- **Preprocessing (`7_preprocessing`)**: `METHOD_ORDER` is the single source of truth so optional
  additions (event rejection, notch, resample, configurable 1/f fit range) keep event/feature-free runs
  **byte-compatible**. Writes the `{file_id}_preprocessing_params.json` sidecar (thresholds actually used
  + `1f_fit_range_hz` + `methods_run`), read back by tool 8. Notch (`cb_notch`, MNE default FIR method,
  50 Hz) now **defaults ON** (a power-line notch is the norm). Untick it to reproduce the byte-identical no-notch
  output, whose `methods_run` entry keeps tool 8 consistent. Resampling (also in tool 6, `cb_resample`) and
  the 1/f fit range (default 2–45 Hz) stay optional and **off/neutral by default**. **PSD smoothing before
  the 1/f fit** (`cb_ssmatch` + optional `cb_ss_median`, two-stage median→LOWESS on the linear Welch PSD,
  oscip parity) is built in and **defaults ON** (the one *not* byte-neutral default: untick it to reproduce the
  pre-smoothing output, and reprocess to homogenise), recorded in the additive `psd_smoothing` sidecar block and
  **honoured by tool 8** when it recomputes PSDs (`smooth_psd_median`/`smooth_psd_lowess` are duplicated in tool 7
  + `qc_rejected_epochs_lib.py`, keep them in sync. 8bis is untouched since it never recomputes a PSD). In the 1/f box the
  **fit parameters** (fit range + smoothing) are grouped one indent shallower than the **thresholds**, and all
  methods' thresholds share that deeper indent so they stay aligned (Variant A, layout-only). `compute_rejection_masks`
  reads the signal **one channel at a
  time** from `epochs` (+ `del raw` after epoching) to bound memory on dense montages. The formulas are unchanged, so
  every output is **byte-identical** to the former full-array version, and the code is format-agnostic (EDF
  source, passes through to the Curry twin). Keep the formulas in sync with tool 8's `qc_rejected_epochs_lib.py`.
  **Configurable epoch length** (`dd_epoch_len`, divisors of 30 ≥ 5 s, **default 30 = classic**): the 30 s
  hypnogram is validated/trimmed at 30 s then **expanded** (`np.repeat(expert_hypno, 30 // epoch_sec)`) so each
  sub-epoch inherits its parent stage. `make_fixed_length_epochs(duration=epoch_sec)`, event mask + context
  companion + Welch window (`min(4, epoch_sec)·sf`) all take `epoch_sec`, sidecar gains `epoch_length_s`
  (30 when absent). Event flagging stays onset-only with an **amber warning** when ≠ 30 s. Thread `epoch_sec`
  consistently across tool 7 → context companion → sidecar → tools 8/8bis/`qc_rejected_epochs_lib` (re-run the
  Curry generator). → SPEC §7 *Epoching (configurable epoch length)*.
  **Output location**: tool 7 writes its `.fif`/sidecar/context under **`derivatives/raw_epo/<subtree>/`**
  (sibling of tool 8's `rejection_manual/` and 8bis's `rejection_auto/`), and deletes stale same-id copies at
  the old `derivatives/<subtree>/` location on reprocess.
  **Deselected channels are kept as `bads`** (→ SPEC §7): in the `.fif` and the sidecar (`bad_channels`),
  out of the average reference and of the flagging (dropped right after epoching, put back at save), but
  **re-referenced by hand** (MNE leaves bads in their original reference). Every tool-7 table stays
  byte-identical. The Curry generator's `OLD_LOAD`/`NEW_LOAD` carries these lines: mirror any edit there.
- **Context-channels companion (`7_preprocessing` block `[H]`, → tool 8)**: when the participant's
  `sub_config` carries a tool-2 `context_channels` block, tool 7 reads **only** those declared EOG/EMG/ECG
  channels from the raw EDF, renames them to role labels (`EOG-L`/`EOG-R`/`EMG`/`ECG`), sets MNE channel
  types, applies the same optional resample, **display-filters per role** (EOG band-pass 0.3–35 Hz, EMG
  high-pass 10 Hz, ECG band-pass 0.5–40 Hz: it is a display companion, so filtering in place is fine and
  makes the tool-8 epoch montage readable), epochs them **identically to the EEG** (`make_fixed_length_epochs`
  duration=30 from t=0: same epoch count regardless of sfreq, hence 1:1 index alignment), and saves
  `{file_id}_context-epo.fif`. **Optional + non-fatal**: no `context_channels` → no companion written, all
  other outputs byte-identical. The Curry generator swaps the reader (`read_raw_edf`→`read_raw_curry`) and
  drops the suffix-dedup line via a dedicated `_make_tool7_curry.py` replacement: re-run the generator after
  editing the `[H]` block.
- **Manually reject flagged epochs (`8_reject_manually*`)**: reads tool-7 `{file_id}_all-epo.fif`
  (+ optional params JSON), **never reloads the raw EDF**, **never modifies tool-7 outputs**. **Section 1 =
  optional data folder + Raw-epochs folder** (like 8bis **minus** the reports picker: tool 8 never reads
  `reports_preprocessing/`). `find_participants` runs on the chosen raw folder so versioned tool-7 runs stay
  apart on the participant dropdown. `data_root` = the selected data folder, else the `derivatives/`
  ancestor's parent. An **"already processed" badge** (`_epoch_decision.tsv` **and** a review/decision TSV on
  disk) is informative only: it blocks nothing. The **Exclude this participant** button writes `8_manual`.
  **The decision is RECOMPOSED, not read**: Section 1 checkboxes (stages / methods / event types, all on by
  default) feed `recompute_reject`, which rebuilds `base_reject` + `reject_method` the way 8bis's
  `build_pair_matrix` does: all ticked ⇒ **identical to tool 7's `reject_flag`** (keep it that way). It
  drives the navigator, Section 2 **and** the saved decision (**out-of-scope stages = `in_scope False`**), and
  gates the cost (`fit_mask=in_scope`, `do_1f=False` when no 1/f method is ticked).
  Per-channel attribution is **recomputed** with tool-7 formulas + persisted thresholds. Section 4 also
  writes an **8bis-style decision record** (`_manualreject_decision.tsv` → `_manualreject_report.html` →
  `global_manualreject_summary.tsv` **rebuilt by globbing the per-file TSVs**, data written before the report). Analysis + plotting live in a **shared module
  `qc_rejected_epochs_lib.py`**: a deliberate exception to the "duplicate helpers" rule. The copied bits
  (`METHOD_ORDER`, palette, custom-stage helpers, Welch-PSD, **PSD smoothing `smooth_psd_median`/`smooth_psd_lowess`**,
  1/f fit) must stay in sync with tool 7. `compute_psds(…, smoothing=info['psd_smoothing'])` re-applies the
  tool-7 smoothing so the recomputed plots/attribution match the flags (`smoothing=None` → byte-identical).
  **Epoch identity**: the UI shows `Epoch 123 (#124) — 01:01:30 — 23:23:54`. Compumedics numbers from 1 at
  the recording start and tool 7 epochs from t=0 with no crop, so `#N = index + 1` (editable offset). The
  clock comes from the `.fif`'s `meas_date` (**no raw reload**). `epoch_idx` **stays the 0-based MNE index**
  everywhere: it is the join key with tools 7/8bis/9. **Review tracking**: `visited` feeds a status banner
  (loud when the decision was *changed*), a `reviewed n/N` HTML strip and the Section-4 strip
  (`plot_review_strip(..., visited=)`, `None` ⇒ byte-identical), persisted as the additive `seen` column +
  `n_seen`/`pct_seen` and re-read at load.
  An **optional "Show EOG/EMG context" toggle** (**default ON**, a no-op when tool 7 wrote no companion) stacks the EOG-L/EOG-R/EMG traces under the
  per-epoch montage, loaded on demand from the `{file_id}_context-epo.fif` companion (`load_context_epochs`,
  aligned by epoch index): still no raw-EDF reload. Absent companion → toggle is a no-op. The navigator
  header names the **event type(s)** that flagged the epoch, read from the `.fif` `evt_<type>` metadata
  columns (no reports-folder access, present only when tool-7 event flagging ran). **Data/reports
  split** (both folders precomputed at load: `S['out_folder']` / `S['reports_folder']`): decision tables →
  **`derivatives/rejection_manual/<subtree>/`**. Reviewed TSVs + `qc2b_report.html` → **`reports_rejection_manual/`**
  (`deriv_root.parent/…`, beside `reports_preprocessing/`).
- **Automatic epoch rejection (`8bis_reject_automatically_voila`)**: **flagging → rejection** decision tool
  (naming scheme: 7 flag → 8 manual → 8bis auto, where 8 and 8bis are alternatives). Voila app. It reads only
  format-agnostic `.fif` + `{file_id}_epoch_channel_rejection.tsv` (both identical from Curry tool 7), so its
  **Curry twin is a verbatim copy**: `tools_curry/_make_tool8bis_curry.py` re-copies + retitles (no
  EDF-specific code to swap). Re-run it after editing the EDF notebook. **Channel-first then epoch**
  (`auto_reject_decision`): drop a channel flagged in
  > `channel_pct` of in-scope epochs, then reject an in-scope epoch flagged in > `epoch_pct` of the *remaining
  good* channels (both default 20 %, editable, computed over the **stages of interest** only). Reads tool-7
  outputs **read-only**. **Section 2 selectable flagging methods + event sub-selection**: `build_pair_matrix`
  recomposes the (epoch×channel) matrix from the ticked `flag_<method>` columns (all on = the stored
  `reject_any`, byte-identical default). Ticking `event` reveals per-canonical-type rows (checkbox +
  `exclude if > N events` threshold + a mean/median-per-file hint) fed by the tool-7 `_event_epoch_flags.tsv`
  (per-type epoch flags, absent → fall back to `flag_event`) and `_event_counts.tsv` (raw counts). A ticked
  type over its threshold **excludes the whole participant** (registry `8bis_auto`, like all-channels-dropped /
  all-epochs-rejected; a re-run replaces these rows). Provenance columns
  `methods_used`/`event_types_used`/`excluded`/`exclude_reason`. Run tally gains `excluded`. Still format-agnostic → **Curry twin stays a
  verbatim copy** (re-run `_make_tool8bis_curry.py`). **Section 1 = two explicit folder pickers (raw-epochs + reports) + a Scan button**
  (+ an optional data-folder chooser that pre-points them): *not* one root scanned recursively, so
  renamed/versioned tool-7 runs (`raw_epo_v2`, `reports_preprocessing_v2`) stay apart. **Data/reports split**:
  `_epoch_decision.tsv` + `_channel_decision.tsv` → **`rejection_auto/`** (`raw_root.parent/…`, beside the raw
  folder; 8bis no longer reads the `.fif`). `_autoreject_decision.tsv` (durable record + global-summary source, written **before**
  the report) + `_autoreject_report.html` + the global summary/report → **`reports_rejection_auto/`**
  (`reports_root.parent/…`, beside the reports folder). No interpolation (channels are dropped, interpolation is deferred). Skip gate
  = `_epoch_decision.tsv` **and** decision TSV present. Global summary rebuilt by globbing the per-file decision TSVs.
  All-channels-/all-epochs-rejected → still get tables + row + report (100 % in the summary). See SPEC §8bis.

- **Spectral features (`9_spectral_features_voila`)**: first **feature tool**: it reads tool 7's `_all-epo.fif`
  + a tool 8/8bis decision folder, never the raw, and computes **every** epoch and channel (aggregates over the
  kept epochs). It saves `_psd_epoch.npz` so **Re-aggregate only** rebuilds everything after a decision change
  without PSD / fit; `finish_participant` is shared by both paths. Night thirds span the whole night. Tool 7's params sidecar is **provenance only** (auto-located, absence
  non-fatal). Three invariants: the **PSD smoothing feeds the specparam fit only** (band powers always use
  the raw PSD). **Band integrals are rectangular sums** (`Σ PSD × df`, so tiling bands give relative powers
  summing to 1). And the **peak model follows specparam's defaults** (`[0.5, 12]`, `min_peak_height 0.1`),
  *not* tool 7's wider `[0.5, 20]` / `0.3`: this divergence is **deliberate, do not align them** (tool 7
  uses the fit as an artefact detector, tool 9 as a decomposition. Aligning them moves 9–12 % of tool 7's
  (epoch × channel) 1/f flags and forces a full reprocess). Tool 9 is **not** part of tool 7's three-copy
  `SpectralModel` sync constraint. Outputs split data → `derivatives/features_spectral/`, reports
  → `reports_features_spectral/`. Database tables are globbed from disk and **padded to the channel union**
  (database figures: `figure_rows` = not excluded, not dropped). → SPEC §9.
- **Any tool computing a PSD**: work in **µV²/Hz** and guard logs with `np.where(psd > 0, psd, np.nan)`,
  and **never** with `psd + 1e-10` on a V²/Hz array (that is a 100 µV²/Hz floor, above most of the sleep spectrum).
  With **multitaper**, pass `normalization='full'`: MNE's `'length'` default is not a density (off by
  `sfreq`, so recordings at different rates stop being comparable). → SPEC *Cross-cutting → PSD units and
  log guards*.

**Curry twins** (→ SPEC *Curry 9 support*)
- **Tool 8's twin = 3 replacements** (`_make_tool8_curry.py`): H1, one prose sentence, and the
  **shared-library import probe** (the EDF notebook probes only `cwd`/`cwd/tools`, which does not resolve
  from `tools_curry/`). Nothing Curry-specific is injected: the high-density support is inherited from the
  shared lib. Re-run the generator after editing the EDF notebook. It asserts `3/3` and that no `.edf` /
  `raw EDF` / `read_raw_edf` string survives.
- **Generated, not hand-edited**: `tools_curry/_make_tool{6,7}_curry.py` regenerate the Curry notebooks
  from the EDF originals by string replacement. Edit the EDF notebook, then **re-run the generator**. New
  code passes through automatically **unless** it sits inside a block the generator string-matches or
  wholesale-replaces (e.g. tool 7's Curry load block, event loader, or `[H]` context reader), in which case
  the change must be mirrored in the generator. (`8bis` is the exception: format-agnostic, so
  `_make_tool8bis_curry.py` is a **verbatim copy + retitle**.)
- **Canonical launch cwd = repo root, and shared-module imports must be cwd-independent**: every tool is run
  from `Inspect_EDF/` (see SPEC *How to run*), but the shared-module import must still resolve from either
  cwd. The Curry twins (1,2,4,6,7) probe `cwd`, `cwd/tools_curry`, `dirname(cwd)/tools_curry` for
  `curry_header.py`. 8bis (EDF + Curry) probes the `tools/` variants for `qc_rejected_epochs_lib.py`. This
  block lives in the **generators** (except tool 1 = hand-edited, no generator): fix it there and
  **re-run the generator**, never the old `os.path.dirname(os.path.abspath('__file__'))` (= cwd only,
  breaks from the repo root).

## Working agreements

- **Developer background**: sleep research engineer, PSG/EEG expert, limited software-engineering experience.
  Prioritize readability and correctness over abstraction or cleverness.
- **No unnecessary abstraction**: explicit, readable functions, three clear lines beat a clever one-liner.
  No new helpers/classes unless the complexity clearly justifies it.
- **Comments**: only where the EEG/PSG domain logic isn't obvious to a non-sleep-researcher (e.g. why a
  500 µV threshold flags clipping, or why stage 4 is remapped to N3).
- **Proactive error handling**: add `try/except` in every new feature/notebook. Fatal step (file I/O,
  epoching) → add to a `failed` list and `continue`. Non-fatal step (re-referencing, report) → `⚠` warning
  and continue. In Voila, wrap every button callback and per-item loop so one failure never crashes the run
  or freezes the UI. Always surface errors via a widget or `print()`.
- **Skip + cumulative-merge**: every per-participant processing tool has a "Skip already processed" checkbox
  (on by default), an "N / M already done" info line, and merge/replace output semantics (cumulative per-row
  files merged on the item id, aggregated summaries regenerated from all per-item files). **Interruption-safe**:
  skip an item only when **both** its report **and** its durable per-item data are on disk (mismatch → ⚠ warn +
  reprocess). Write per-item data **before** the report. Rebuild every cumulative/aggregated table by **globbing
  the per-item files on disk**, never from an in-memory `attempted_ids` merge. See *Cross-cutting procedures* in
  SPEC.md. Reference impls: tools 6 & 7.
- **Normalize path comparisons**: wrap **both** sides of any path/stem/filename string comparison
  (`==`, `in`, `.isin()`, set/dict membership) in `os.path.normcase(...)` **at the comparison only** (keep
  the stored/displayed value original). Prevents skip checks silently failing on `C:`/`c:` and `/`/`\`.
- **Language**: all user-facing strings in notebooks and all of SPEC.md must be in English.
- **Tests** (`tests/`, pytest, → SPEC *How to run → Tests*): the notebooks are driven headlessly through
  their widgets on two datasets. **synthetic** (generated, golden versioned) and **real** (tools/test_data, a
  real participant's PSG: **never versioned**, `tests/golden/real/` git-ignored). Three levels: **quick**
  (`-m quick`), **standard** (`-m "not full"`), **full** (everything). **At each implementation, propose a level
  with a recommendation before running**, and select the tools concerned with the `toolN` markers (modified tool
  + downstream + users of a modified shared module). The chain is reused and rerun from the first changed tool
  (`tests/chaincache.py`): a writing test works on a **copy**, never on the chain's folders. Regenerate a golden
  (`tests/make_golden.py`) **only on trusted code, before** a change, and list an intended output change in
  `INTENDED_CHANGES` with its own dedicated test. Beyond that, validate on real EDF files; no new test framework.

## Operational notes (this machine)

- **OS**: Windows 10 + Conda, PowerShell default. Use forward slashes in Python path strings (MNE/pathlib
  handle them on Windows).
- **Running Python**: `conda` is **not** on PATH. miniforge3 is at `$env:LOCALAPPDATA\miniforge3`. Always use
  the explicit interpreter, not `conda run`, `conda activate`, or bare `python`:
  ```powershell
  & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" tools/my_script.py
  ```
- **Editing large notebooks**: `Read`/`Edit`/`NotebookEdit` are blocked or size-capped on big `.ipynb`
  files: use the **rename-to-`.txt`** method (and the JSON-surgery fallback) in SPEC.md →
  *Repository structure → Editing the larger notebooks*.

## Process for new tasks

- **Read SPEC.md first**, then **ask clarifying questions** (inputs/outputs, edge cases, UI layout, backward
  compatibility, interaction with existing tools, output file names) before planning: resolve ambiguities
  first.
- **Plan-first**: draft the plan with the most capable available model, write it to a temporary
  `tools/plan_<task>.md` for the user to read/annotate, and start implementing **only** after explicit
  confirmation. Delete the plan file when the task is done and confirmed.
- **Keep docs in sync (details go to SPEC.md, not CLAUDE.md)**: add any directly-imported package to
  `environment.yml` (conda-forge if available, else pip). Put the *details* of any design decision (offsets,
  thresholds, widget names, defaults, output files) in **SPEC.md** (the *Cross-cutting procedures* section
  or the relevant tool section) and **propose** those SPEC edits for confirmation. Only update **CLAUDE.md**
  when a *cross-cutting rule* or a *sync constraint* changes, and then only as a short reminder + pointer to
  SPEC (keep CLAUDE.md short).
