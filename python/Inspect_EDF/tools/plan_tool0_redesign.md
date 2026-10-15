# Plan: tool 0 redesign, the whole pipeline on one participant

**Status**: validated by the user on 2026-10-09 (decisions D1-D8 below), **in progress: phase 0**.

**Goal**: rebuild `tools/0_live_explore_1file_voila.ipynb` so that it walks **one participant** through every
step of the pipeline (tools 1 to 9), each step in its own section, with all outputs shown **inside the
app**. Purpose: **diagnose one participant** and **train a new user** on the pipeline.

**Core rule**: tool 0 **never re-implements a computation**. Each section calls the very functions the
batch tool calls, moved into shared libraries (`tools/*_lib.py`). Same code, so same numbers, by
construction. The current tool 0 shows why: its 1 925 lines of copies drifted from tools 6 and 7 within a
few weeks (no `HIST_CLIP_PCT`, no notch / PSD smoothing / event method, a 2-30 Hz 1/f fit instead of
2-45 Hz, custom stages looked up in the wrong folder since the study-folder layout).

---

## 0. Decisions

### Taken (user, 2026-10-09)

1. Purpose: diagnosis of one participant + training of new users.
2. Architecture B (shared libraries), **on condition of before/after tests** proving that the batch tools
   behave exactly as before.
3. All sections, built **step by step**, with a user validation after each phase.
4. Rescoring and annotations stay in tool 0, annotations in tool 8's format.
5. Interactive canvas (ipympl / Plotly) later. The tool must stay interactive (fast redraws, buttons,
   keyboard).

### Design decisions (validated 2026-10-09, user notes integrated)

| # | Question | Decision |
|---|---|---|
| D1 | Where does the epoch explorer sit? | A permanent panel under the sections, not a step (see §1.2) |
| D2 | Can tool 0 change keep/reject decisions? | No, read-only: changing a decision is tool 8's job. **But the decision in force feeds the next sections**: §8 (spectral features on the kept epochs) and the explorer's decision strip. §7 offers its source: **live** (8bis rule computed in §7, default) or **saved on disk** (the participant's tool 8 or 8bis decision tables, offered only when they exist and their epoch count matches §6) |
| D3 | What does tool 0 write? | `_rescored.txt` + `_rescore_log.tsv` (beside the hypnogram), `{id}_manual_events.tsv` (beside the recording), and a **progressive inspection report** (§1.4). The `_live_rejection_mask.tsv` export is dropped (no tool reads it) |
| D4 | What goes into a library? | Format-agnostic computation only (metrics, flags, decisions, features, plots). EDF/Curry reading stays in the notebooks where the Curry generators expect it (§3.3) |
| D5 | EDF reading helpers (header parser, events chain, hypnogram loading) | A new `edf_io_lib.py`, used first by tools 0 and 5 (no Curry twin). **Then the other tools switch to it in a dedicated phase (phase 8)**, adapting their Curry generators. Until then `test_tool0.py` catches any divergence through the numbers |
| D6 | Show the participant's exclusion status? | Yes, in the header: included / excluded, by which tool, why (registry, read-only) |
| D7 | Training aid | Yes: each section opens with a short box: what this step does, which batch tool runs it on the whole database, what to look at |
| D8 | Signal in the explorer before §6 has run | "As recorded" + optional display band-pass. After §6: Raw / Preprocessed toggle (tool-7 settings) |

---

## 1. Target layout

### 1.1 Page structure

```
Tool 0: one participant, the whole pipeline
  Study folder [picker]   Participant [dropdown]   status: included | EXCLUDED (8bis_auto: ...)
  ┌ Tabs ───────────────────────────────────────────────────────────────────────┐
  │ 1 Load │ 2 Header │ 3 Scoring & events │ 4 Macro │ 5 Quality │ 6 Preproc │ 7 Decision │ 8 Spectral │
  │  [training box]  [parameters]  [▶ Run]  badge: not run / up to date / OUT OF DATE  │
  │  outputs (tables, figures, "Show in explorer" buttons)                             │
  └─────────────────────────────────────────────────────────────────────────────┘
  ▼ Epoch explorer (collapsible, always available once §1 is loaded)
     navigator: spectrogram + hypnogram + event ticks + [flag strip] + [decision strip]
     ◀ Prev  [slider]  Next ▶  Next stage ⇥    Epoch 123 (#124), 01:01:30, 23:23:54, N2
     montage (fixed clinical scales) + per-epoch PSD / detail
     Rescore: W N1 N2 N3 R MT ?   Annotate: [type] [comment] [+ Add]   💾 Save rescored
```

Built with `ipywidgets.Tab` (one tab per step) so the page stays short. Each tab is computed **only when
its Run button is pressed** (lazy, keeps the app responsive). A shared `STATE` dict holds the participant
in memory (raw signal, hypnogram, events, then each section's results).

### 1.2 Why the explorer is not a step (answer to the ordering question)

The explorer serves two moments of the pipeline:
- **upstream**: scoring review (rescoring) and adding missed events. Their results feed §4 (macrostructure)
  and §6 (per-stage amplitude thresholds, event flagging);
- **downstream**: looking at the epochs that §6 flagged and §7 rejected.

So it cannot sit at one place in the sequence. It becomes a **permanent viewer** that gains layers as the
sections run: events (§1), GSSC hypnodensity strip (§3), flag heatmap strip (§6), decision strip (§7).
Every list in a section (events, flagged epochs, rejected epochs, dropped channels) has a **"Show in
explorer"** button.

### 1.3 Out-of-date sections

Rescoring changes the hypnogram, so the sections computed from it are no longer valid. A dependency map
marks them **OUT OF DATE** (orange badge, "Run again"), nothing is recomputed automatically (keeps the app
responsive):

| Change | Sections marked out of date |
|---|---|
| hypnogram rescored | 4, 6, 7, 8 |
| annotation added / deleted | 4 (event indices), 6 (when event flagging includes manual events), 7, 8 |
| §6 run again | 7, 8 |
| §7 run again, or its decision source changed | 8 |
| new participant | everything reset |

### 1.4 Progressive inspection report (D3)

A record of what was looked at and changed during the inspection, kept up to date as you go (checkbox
**Keep an inspection report**, on by default):
- `{file_id}_inspection_log.tsv` (data, **appended across sessions**): one timestamped row per action
  (section run + its parameters, rescoring saved, annotation added / deleted, decision source changed);
- `{file_id}_inspection_report.html` (human, **rewritten** after each action): participant header +
  exclusion status, the whole log, then each section's latest tables and figures with its status
  (up to date / OUT OF DATE / not run in this session).

Both go to `<study>/reports_live_explore/<subtree>/` (mirrors the recordings' tree, like every other
report). Data first, report second, as everywhere in the toolkit. A failure to write them is a `⚠`
warning, never a crash.

---

## 2. The sections

| § | Title (batch tool) | Shows | Computation from |
|---|---|---|---|
| 1 | Load & configuration (tools 2, 3, 4) | Channels with their roles: EEG from the tool-2 remap, EOG-L/R / EMG / ECG from `context_channels`, detection fallback when no JSON. Hypnograms found by suffix (expert / remapped / GSSC, user picks one). Events TXT → CSV → XML (+ manual), canonical labels via `event_remap.json`. Custom stages from `<study>/config_param/` | `edf_io_lib` |
| 2 | Recording header (tool 1, 1bis check) | Per-channel table: sampling rate, unit, physical / digital range, prefiltering. Duration, start date-time, anonymization check | `edf_io_lib` (tool 1's parser + `check_anonymization`) |
| 3 | Scoring & events (tools 3, 3bis, 4) | Hypnogram (raw labels → canonical, unmapped labels listed). Event timeline + count per canonical type. **Optional** GSSC: hypnodensity, agreement with the expert (kappa, confusion matrix), and "use as rescoring suggestion" in the explorer | `autoscoring_lib` (from 3bis) |
| 4 | Sleep macrostructure (tool 5, 5bis ranges) | The `METRICS` table with the reference ranges and out-of-range flags, lights-off/on source, AHI / ODI / T90 when the companions exist | `macrostructure_lib` (from 5; 5bis imports its `METRICS` / `REFERENCE_RANGES`) |
| 5 | Signal quality (tool 6) | Current §2: per-electrode metrics + flags, histogram, time series, hypnospectrogram, whole-night and per-stage PSD | `quality_lib` (from 6) |
| 6 | Preprocessing & flagging (tool 7) | Tool-7 parameters with tool-7 defaults (notch, band-pass, resample, epoch length, methods incl. `event`, 1/f range, PSD smoothing). Flag heatmap, per-stage and per-method counts, per-channel badness | `preprocessing_lib` (from 7) + `qc_rejected_epochs_lib` |
| 7 | Rejection decision (tool 8bis, tool 8 = manual alternative) | Decision source (D2): live 8bis rule (20 % / 20 % defaults) or the tool 8 / 8bis tables on disk. Dropped channels, rejected epochs per stage, decision heatmap, would-be exclusion reasons. The decision in force feeds §8 and the explorer | `autoreject_lib` (from 8bis) |
| 8 | Spectral features (tool 9) | On the kept epochs: PSD per stage, band powers, aperiodic exponent / offset per stage and over the night, peaks | `spectral_lib` (from 9) |
| - | Epoch explorer (tool 0's own + tool 8's display) | See §1.1. Montage, scales, epoch identity, annotations: `qc_rejected_epochs_lib` (`plot_epoch_montage`, `DISPLAY_SCALE_UV`, `epoch_labels`, `save_manual_events`) | `qc_rejected_epochs_lib` |

Annotations written by tool 0 use `MANUAL_EVENT_COLUMNS` with `source = 'manual_tool0'`, so tool 7
(`cb_manual_events`) and tool 8 read them like tool 8's own.

Not brought into tool 0: anonymization itself (1bis, a file operation), 5bis's database distributions,
the database-level reports of tools 6 and 9, manual keep/reject (tool 8, see D2).

---

## 3. Shared libraries (architecture B)

### 3.1 Rules

- A library holds **functions and constants only** (no widget, no file picker). Moved **verbatim**: the move
  commit changes no logic, only where the code lives.
- Imported through the cwd-independent probe already used for `qc_rejected_epochs_lib`
  (`cwd`, `cwd/tools`, `dirname(cwd)/tools`).
- `tests/chaincache.py` `TOOL_DEPS`: each new library is added to the tools importing it (so a library edit
  reruns the right part of the chain).
- Every "duplicated, keep in sync" rule that becomes an import is **removed** from CLAUDE.md / SPEC.md
  (e.g. 5bis's verbatim copy of `METRICS`, `get_phys_bounds_uV` in tools 0 & 6, the hypnospectrogram in
  0 & 6, the palette remapped in tool 0, `detect_channel_types` "defined in tool 0").
- `qc_rejected_epochs_lib.load_participant(fif_path)` is split into `participant_from_epochs(epochs)` +
  the file wrapper, so tool 0 can feed its in-memory epochs to the same plotting code (tool 8 unchanged).

### 3.2 The libraries

| Library | Moved from | Content (main functions) | Imported by |
|---|---|---|---|
| `edf_io_lib.py` | tool 1 (parser, anonymization check), tool 5 (events chain, the superset with the `Desaturation` column), tool 0 (`detect_channel_types`), tools 3/6/7 (hypnogram suffix detection, custom-stage helpers) | `read_edf_header_custom`, `read_edf_start_datetime`, `load_events`, `load_event_remap`, `augment_remap_for_manual`, `detect_channel_types`, suffix detection, `drop_suffix_duplicates`, `adapt_remap_dict_to_suffixes`, `get_phys_bounds_uV` | 0, 5 |
| `autoscoring_lib.py` | 3bis cell 2 | `gssc_combination_logits`, `run_gssc`, `compare_scoring`, `plot_hypnodensity`, `plot_confusion` | 0, 3bis |
| `macrostructure_lib.py` | 5 cells 1, 3, 4 | `METRICS`, `REFERENCE_RANGES`, `compute_macrostructure`, lights resolution, `compute_event_indices`, `compute_odi`, `compute_t90`, plots | 0, 5, 5bis |
| `quality_lib.py` | 6 cell 1 (per-channel part) | `compute_signal_metrics`, `flag_channel`, `HIST_CLIP_PCT`, per-channel plots, hypnospectrogram | 0, 6 |
| `preprocessing_lib.py` | 7 cell 2 (computation part) | `METHOD_ORDER`, `compute_rejection_masks`, `compute_event_epoch_mask`, `smooth_psd_median` / `smooth_psd_lowess`, `build_*` tables | 0, 7 (+ `qc_rejected_epochs_lib`, which then imports the smoothing instead of copying it) |
| `autoreject_lib.py` | 8bis cell 1 | `build_pair_matrix`, `auto_reject_decision`, `build_epoch_decision`, `build_channel_decision`, `plot_decision_heatmap` | 0, 8bis |
| `spectral_lib.py` | 9 cells 2, 3 | `compute_psd_array`, `fit_spectra`, `compute_band_measures`, `aggregate_*`, per-participant plots | 0, 9 |

The UI cells (pickers, Run loops, reports, Skip / merge logic) stay in each batch notebook.

### 3.3 Curry twins: what can move

| Tool | Curry generator | Impact |
|---|---|---|
| 6 | `_make_tool6_curry.py` rewrites `drop_suffix_duplicates`, `adapt_remap_dict_to_suffixes`, `get_phys_bounds_uV` and the signature of `compute_signal_metrics` | The generator is adapted to the lib (the twin calls `compute_signal_metrics` without physical bounds), then rerun |
| 7 | `_make_tool7_curry.py` rewrites the suffix helpers, `read_edf_start_datetime`, the load block, the `[H]` reader | Those blocks **stay in the notebook**. The computation functions of cell 2 already pass through unchanged, so they can move |
| 8bis | `_make_tool8bis_curry.py` copies the notebook verbatim | The twin imports `autoreject_lib` through the same probe as `qc_rejected_epochs_lib` |
| 1, 3bis, 5, 5bis, 9 | none (tool 1's twin is hand-edited and is not touched) | Free |

---

## 4. Tests: proving the batch tools did not change

| Layer | What | When |
|---|---|---|
| A. Baseline | Full suite on today's code (synthetic + real). Must be green before anything moves, otherwise we stop and look | Phase 0 |
| B. Before / after check (new `tests/refactor_check.py`) | Runs the affected tools on synthetic + real **before** the move and keeps the whole output tree, runs them again **after**, and compares **everything**: TSV / JSON / TXT values, `.fif` signal arrays + metadata + info, `.npz`, xlsx sheets, PNG bytes, HTML with dates and paths masked. Also covers the tools outside the chain (1, 3bis). Expected report: **empty**. Floats use the golden tolerance (tool 9's MKL noise, commit 844cee2). Any difference stops the phase until it is explained | Each library move |
| C. Existing checks | Goldens (`test_pipeline.py`), `test_curry_twins.py`, quick (every notebook executes) | Each phase |
| D. Tool 0 consistency (new `tests/test_tool0.py`, marker `tool0`) | Drives tool 0 on a synthetic participant section by section and compares with the chain outputs of **the same participant**: §4 = tool 5's `sleep_metrics`, §5 = tool 6's `quality_metrics`, §6 = tool 7's `epoch_channel_rejection`, §7 = 8bis's epoch / channel decision, §8 = tool 9's tables. This is the permanent proof "same code, same numbers" | Grows section by section |

Layer B goes further than the goldens on purpose: they summarize a `.fif` as channels + metadata and skip
HTML / PNG, but moving tool 7's code must not change a single signal value.

---

## 5. Phases (each ends with your validation)

Every phase follows the same steps:
1. snapshot "before" (`refactor_check.py`);
2. move the code verbatim into the library, the batch notebook imports it (+ Curry generator if needed);
3. snapshot "after" → empty diff, goldens + Curry twins green;
4. build the tool-0 section, extend `test_tool0.py`;
5. you test in Voila and validate;
6. I stage (`git add`), you commit. Two commits suggested per phase: "extract X lib (no behaviour
   change)" then "tool 0 §N".

| Phase | Content | Libraries | Tests |
|---|---|---|---|
| 0 | Baseline (full suite from scratch on HEAD), `refactor_check.py` + its calibration (a second identical run must give an empty report), `tool0` marker in `pytest.ini` | none | full, twice |
| 1 | Skeleton (header, tabs, explorer panel, `STATE`, out-of-date badges, inspection log + report), §1 Load, §2 Header, explorer moved to tool 8's display (fixed scales, epoch identity, context channels, annotations), rescoring kept. `test_tool0.py` created | `edf_io_lib` (tool 0 only at this stage), split of `load_participant` in `qc_rejected_epochs_lib` | standard tool 8 + tool 8bis + tool 9 (qc lib) + quick |
| 2 | §3 Scoring & events, GSSC optional (lazy torch import, the app works without torch) | `autoscoring_lib` | refactor check on 3bis + quick |
| 3 | §4 Macrostructure | `macrostructure_lib`, tool 5 switches to `edf_io_lib` | standard tool 5 + 5bis + downstream |
| 4 | §5 Quality | `quality_lib` + tool 6 Curry generator | standard tool 6 + downstream + Curry twins |
| 5 | §6 Preprocessing & flagging | `preprocessing_lib` (+ `qc_rejected_epochs_lib` imports the smoothing) | full (tool 7 = the heart of the chain) |
| 6 | §7 Rejection decision | `autoreject_lib` + 8bis Curry twin | standard 8bis + 8 + 9 + Curry twins |
| 7 | §8 Spectral features | `spectral_lib` | standard tool 9 |
| 8 | The other tools switch to `edf_io_lib` (D5): events chain / start date-time / suffix helpers / hypnogram suffix detection in tools 3bis, 4, 6, 7; header parser in tools 1, 1bis, 2; custom-stage helpers in tools 3, 6, 7, 9 and `qc_rejected_epochs_lib`;
tool 3's mid-night `?` rule; tool 6's load block (EDF helpers + bounds lines → `edf_io_lib.phys_bounds_for_remap`); the I/O helpers duplicated in `autoscoring_lib` (`drop_suffix_duplicates`,
`adapt_remap_dict_to_suffixes`, `read_edf_start_datetime`, `load_json_lenient`, `load_custom_stages`). Curry generators 2, 4, 6, 7 adapted, twins regenerated. Where two tools expose different interfaces (tool 4's events loader returns tuples, tool 7's a DataFrame) the tool keeps a thin adapter, decided then | `edf_io_lib` everywhere | full + refactor check on tools 1, 3bis + Curry twins |
| 9 | Docs: SPEC §0 rewritten, *Cross-cutting → Shared computation libraries*, CLAUDE.md (rule change + sync rules removed), TESTS.md, TODO file trimmed to the canvas item. Plan deleted | none | quick |

The test level is proposed again at the start of each phase, as usual.

### Phase 0 results (2026-10-09 / 10-12)

- Baseline: full suite green on HEAD 7b8c6b4 (157 passed, 19 skipped by design), outputs saved
  (`refactor_check.py save baseline`, 302 files, 614 MB).
- Found and fixed on the way: `test_tool9.same()` applied the 1e-5 dB tolerance of commit 844cee2 only to
  columns *ending* in `_ap_removed_db`, not to the per-stage aggregates (`..._ap_removed_db_from_log_median`),
  so MKL noise (2e-7 dB) failed the real dataset.
- Calibration (two identical fresh runs compared): every data file byte-identical (tables, JSON, `.fif`
  signal, `.npz`, xlsx), tool 9 included. Only tool 6's pooled-distribution figures differed: their dots were
  jittered by an unseeded `np.random.uniform`. **Fixed** (user decision): `np.random.default_rng(0)`, the
  report is now reproducible (Curry twin regenerated).
- The first run hung for 1 h: Windows connected standby (screen off, even on AC) pauses the programs.
  `conftest.py` now keeps the display on during a pytest session (`SetThreadExecutionState`, user decision).
- Measured durations with the PC awake: both chains from scratch + `test_pipeline.py` = 8 min, full suite
  with the chains reused = 17 min.

### Phase 1 results (2026-10-12)

- `qc_rejected_epochs_lib.load_participant` split into `participant_from_epochs` (+ the file wrapper):
  tools 8bis / 8 / 9 rerun on both datasets, `refactor_check compare seeded` = **no difference** (302 files).
- New `tools/edf_io_lib.py` (tool 0 only for now): tool 1's header parser (+ optional `fallback_encoding`,
  tool 0 passes latin-1) and its header table / 500 µV dynamic-range rule, tool 0's `detect_channel_types`,
  tool 2's `detect_context_channels`, tool 5's TXT/CSV/XML event chain, tool 7's reading helpers, suffix
  rules, hypnogram length rule and `[H]` context reader. Found on the way (left as is, documented in SPEC):
  the shared header parser overwrites the global `reserved` field (EDF+C/D) with the per-channel list.
- Tool 0 rewritten (12 cells): header (study folder, participant, exclusion status, report checkbox),
  tabs 1 Load and 2 Header built, tabs 3-8 = their *About this step* box + "phase N", explorer panel (tool 8's
  montage + detail panel, rescoring, annotations `manual_tool0`, display band-pass), out-of-date badges,
  inspection log + report, unsaved-rescoring guard (first reload / participant change refused).
- Explorer trace colours = **preview** at tool 7's default thresholds on the displayed signal (said in the
  panel) until step 6 exists.
- Tests: `test_tool0.py` (step 1 = tool 7's epochs / stages / channels / mapped events on every participant;
  step 2 = tool 1's known answer; explorer writes only its files; context traces = tool 7's companion),
  `test_edf_io_lib.py` (quick).

### Phase 2 results (2026-10-12)

- `tools/autoscoring_lib.py`: tool 3bis's cell 2 verbatim (checked by script) + its guarded GSSC import +
  `score_night` (the scoring lines of `score_one`, unchanged). Tool 3bis imports them by name; its notebook
  lost 667 lines. Before / after on tool 3bis (new `tests/outside_chain.py`, both datasets, 86 files):
  calibration (two identical runs) and the move both = **no difference**.
- `edf_io_lib.mid_recording_unscored` = tool 3's mid-night `?` rule (tool 3 switches to it in phase 8).
- Tool 0 step 3: overview figure (hypnogram + events by type), hypnogram checks + crosstab against the other
  hypnogram files of the recording, events per tool-4 label and per stage; GSSC read from tool 3bis or redone
  here (`autoscoring_lib` imported on first use), tool 3bis's figures and agreement; explorer GSSC layer
  (dashed line, status, Next disagreement, Use GSSC's stage). Rescoring / annotation now also mark step 3.
- Tests: step 3 events count; GSSC = tool 3bis (stages, probabilities, κ; marker `tool3bis`).
- User feedback after phase 2 (done): explorer laid out like tool 8 (navigation above the figures, decision
  slot, rescoring, annotations, display settings under them), navigator at the montage's width and margins.

### Phase 3 results (2026-10-12)

- `tools/macrostructure_lib.py`: tool 5's constants + helpers (setup cell), computations (cell 3), figures /
  HTML (cell 4), verbatim, + `compute_recording` (the computation lines of `process_participant`, lifted with
  its UI references turned into parameters). Tool 5's event chain was identical to `edf_io_lib`'s (checked
  function by function): the library imports it, so tool 5 now uses `edf_io_lib`. Tool 5 imports 56 names
  explicitly; tool 5bis imports `METRICS` / `METRIC_KEYS` / `METRIC_INFO` / `PROVENANCE_KEYS` / `REFERENCE_RANGES`
  (checked identical before). `chaincache.TOOL_DEPS`: tools 5 / 5bis depend on both libraries.
- Tool 0 step 4: tool 5's computation on the hypnogram in use, metrics with ranges and a "tool 5 on disk"
  column, checks, provenance, tool 5bis's last decision. Annotations no longer mark step 4 out of date (tool 5
  never reads them).
- Tests: step 4 = tool 5 on every participant (marker `tool5`).

### Phase 4 results (2026-10-13)

- `tools/quality_lib.py`: tool 6's shared-function cell verbatim (metrics, flags, figures, custom-stage helpers,
  hypnospectrogram, "All electrodes" figures; the EDF helpers and `generate_dataset_overview` stay in the
  notebook) + six functions holding the computation lines of `run_analysis` (`hypno_vector`, `prepare_signal`,
  `check_hypno_length`, `compute_channel_metrics`, `display_limits`, `compute_channel_spectrograms`): checked by
  script, every line is tool 6's except signatures, returns and the progress callbacks. Only UI change: the Load
  segment of the progress bar ticks once after resampling + high-pass (was once after each).
- Curry generator: the library is **copied** into the twin's cell 1 in place of the import, the `bounds_pct`
  removals and the ±500 µV cap target that copy (indentation adapted), the load block defines empty bounds.
  Old twin vs new twin: 68 definitions AST-identical, only `run_analysis` differs, as in the EDF notebook.
- Before / after: whole chain (302 files) and tool 6 with resampling 128 Hz + high-pass 0.5 Hz
  (`outside_chain.py 6variant`, new, 62 files) = **no difference**.
- `edf_io_lib.phys_bounds_for_remap`: tool 6's bounds lines (header-only read) for tool 0; tool 6 keeps its
  copy until phase 8.
- Tool 0 step 5: tool 6's run (thresholds and options exposed with tool 6's defaults), overview box, table with a
  "tool 6 on disk" column, per-stage tables, overview figures, one-channel panel + "Show in the explorer". A
  rescoring now marks step 5 out of date.
- Tests: step 5 = tool 6 on every participant (marker `tool6`).

### Between phases 4 and 5 (2026-10-13 / 14)

- Display harmonization (user request): REM in thin red on every hypnogram, one spectrogram palette (viridis) and
  colour rule (tool 6's, robust to data gaps) with a colour scale on each, explorer events strip + GSSC legend /
  toggle; tool 5's per-recording night report (hypnogram + scored events), linked from its database report.
- The user's commit 7eaec8d renamed `derivatives/raw_epo` → `preproc_epo` and `reports_preprocessing` →
  `reports_preproc` (tool 7 migrates the old folder once).

### Phase 5 results (2026-10-14)

- `tools/preprocessing_lib.py`: tool 7's shared-function cell minus the EDF reading (methods registry, custom
  stages, event mask, PSD smoothing, `compute_rejection_masks`, heatmap, every table / HTML builder), verbatim, +
  the computation lines of `run_preprocessing` in small functions (`preprocess_raw`, `prepare_hypnogram` →
  `hypnogram_epochs`, `make_epochs`, `set_bads_aside`, `epoch_psd_for_1f`, `event_flags`, `epoch_metadata`) whose
  fatal paths raise `StepFailed(reason, message)` with the former strings. Checked by script both ways. The Curry
  twin imports the library (generator unchanged; 137 definitions unchanged, 25 moved ones AST-identical).
  `qc_rejected_epochs_lib` imports the two smoothing functions (they were AST-identical copies).
- `select_montage_channels`: tool 7's selection lines (its load block stays in the notebook), for tool 0.
- Before / after: whole chain and `outside_chain.py 7variant` (resampling 128 Hz, 10 s epochs, no notch,
  LOWESS-only, events).
- Tool 0 step 6: every tool-7 setting with tool 7's defaults, channel pre-selection from tool 6's table (or step 5),
  heatmap, tables, "tool 7 on disk" line; the explorer gains the preprocessed signal + a flags strip, and its trace
  colours use the step-6 thresholds (the former `L.DEFAULT_THRESHOLDS` are not tool 7's amplitude defaults).
- Found: SPEC's amplitude-threshold table and `qc_rejected_epochs_lib.DEFAULT_THRESHOLDS` (tool 8's fallback without
  a sidecar) say W 300 / N2-N3 200 µV, tool 7 uses W 250 / N2-N3 300 µV. SPEC corrected; the qc-lib fallback set
  to tool 7's values (user decision, 2026-10-14).
- Tests: step 6 = tool 7 on every participant (marker `tool7`).

### Between phases 5 and 6: user remarks (2026-10-14), **done** (tool-0 tests 15 passed)

Found by remark 4's role check: `edf_io_lib.detect_channel_types` classed **A1 / A2** (auricular references,
common in PSG) as "other", so the synthetic A2 → M2 reference showed "⚠ used as EEG". Tools 1 and 2 already
had A1 / A2 (only the library's copy, taken from the former tool 0, lacked them): added there (user, same day);
only tool 0 uses it, no batch output changes. Follow-up (user): the per-channel inspection is drawn ahead,
right after the overview figures, so it opens at once.

### Phase 6 results (2026-10-14)

- `tools/autoreject_lib.py`: tool 8bis's shared-function cell (decision rule, tables, overview HTML, heatmap),
  verbatim, minus the UI-bound `registry_folder` / `DECISION_DIRNAME` / `S`. Two readers split so tool 0 passes its
  step 6's tables: `build_pair_matrix` → `pair_matrix_from_tables`, `read_event_counts` → `event_counts_from_table`.
  Lifted from the Scan / Run cells: `ordered_stages`, `default_reference_stages`, `event_exclusion_reason`,
  `registry_reason`. Built by script, checked both ways. The Curry twin (verbatim copy) imports it; its
  `test_curry_twins` marker `tool7_deselected` (moved to the library) became `from autoreject_lib import`.
- Tool 0 step 7: source = live 8bis rule (8bis's widgets and defaults, this recording's event counts) or the tool 8
  / 8bis tables on disk (offered when their epoch count equals step 6's); would-be exclusion, "tool 8bis on disk"
  line, 8bis's overview + heatmap. Explorer: decision strip (kept green / rejected red, 8bis's colours) and the
  read-only `box_decision`. Step 6 keeps its event counts (8bis's exclusion threshold).
- Tests: step 7 = tool 8bis on every participant (marker `tool8bis`).
- User decision (2026-10-14): a keep / reject decision is changed in tool 8 only; tool 0 stays read-only (D2 kept).

### Phase 7 results (2026-10-14)

- `tools/spectral_lib.py`: tool 9's setup constants, its computation cell and its plotting / report cell, verbatim
  (the two cells are gone from the notebook). `load_decision` → `decision_from_tables` split. Lifted from the run
  cell: `stage_order`, `prepare_epochs`, `fit_aperiodic`, `usable_bands`, `participant_tables`,
  `participant_report_items` (logs through `warn`, progress through `progress`). Built by script.
- `preprocessing_lib.put_bads_back`: tool 7's save lines (bad channels back in, original order, marked bad); tool 7
  calls it, its Curry twin regenerated (generator unchanged).
- Tool 0 step 8: tool 9's settings with its defaults (bands and peak model at tool 9's defaults), on step 6's
  epochs as tool 7 saves them (bads back + a temporary `.fif` round trip: single precision, so tool 9's numbers to
  the last digit) and step 7's decision (`auto` / `manual`, else `raw` / `raw_tool7flags`); "tool 9 on disk" line,
  exponent per stage, tool 9's report pieces. Step 6 keeps the deselected channels' epochs.
- Tests: step 8 = tool 9 on every participant (marker `tool9`), every per-participant table.

### Phase 8 plan (2026-10-14), **validated: reduced scope + a sync test** (user, same day)

Custom-stage helpers in `edf_io_lib` (user). The copies left in tools 2, 4, 6 and 7 get a quick test
(`tests/test_edf_io_sync.py`) that fails as soon as one differs from `edf_io_lib` (code compared, docstrings
ignored): "keep in sync" becomes an automatic check.

**Inventory** (every function of `edf_io_lib` against its copies, compared as code, docstrings ignored):
- identical: `extract_filter_value` (tools 1, 2), `check_anonymization` (1, 1bis), `drop_suffix_duplicates` /
  `adapt_remap_dict_to_suffixes` (6, 7, `autoscoring_lib`), `get_phys_bounds_uV` (6), `read_edf_start_datetime` /
  `event_companion_paths` / `_TXT_DUR_RE` (4, 7, `autoscoring_lib`), the `_events_df_*` readers and
  `augment_remap_for_manual` (7); `load_custom_stages` (6 libraries + tool 3), `custom_stage_style` (3 libraries),
  `parse_custom_field` (5 libraries; tool 3's differs by a `str()` only).
- different on purpose: `read_edf_header_custom` / `detect_encoding` (the library's = tool 1's + an optional latin-1
  fallback, off by default: same behaviour; tool 1bis reads only the identity fields), `_events_df_from_xml` (the
  library's adds tool 5's `Desaturation` column), tool 4's `load_events` (names only, a list of tuples: another
  interface), `load_json_lenient` in `autoscoring_lib` and tool 7's `load_event_remap` (same behaviour, other code).
- **Curry constraint**: the generators of tools 2, 4, 6 and 7 rewrite exactly these reading blocks (tool 2's header
  scan, tool 4's and tool 7's event chain, tool 6's load block). A twin cannot import `edf_io_lib` as it is (its probe
  looks for `tools_curry/`), and the EDF-only readers (start datetime, companions) have Curry counterparts.

**Proposed scope** (one before-snapshot, the four moves, one verification, one commit each):
- **8a, tool 1** (Voila + both batch scripts) imports the header parser, `extract_filter_value` and
  `check_anonymization` from `edf_io_lib`; tool 1bis imports `check_anonymization` and keeps its own reader of the
  identity fields (it reads exactly what it overwrites). Tool 1's Curry twin (its own reader) is untouched. Check:
  tool 1's known answers + a before / after run of tool 1 (`outside_chain.py` gains a tool-1 runner).
- **8b, `autoscoring_lib`** (tool 3bis) imports its I/O helpers from `edf_io_lib` (`drop_suffix_duplicates`,
  `adapt_remap_dict_to_suffixes`, `read_edf_start_datetime`, `load_json_lenient`, `load_custom_stages`). Check:
  `outside_chain.py 3bis`.
- **8c, custom stages**: one copy of `load_custom_stages` / `parse_custom_field` / `custom_stage_style`, in
  `edf_io_lib` (section *Study configuration*), imported by `macrostructure_lib`, `preprocessing_lib`,
  `qc_rejected_epochs_lib`, `spectral_lib` (its `stage_style` is another function: kept), `autoscoring_lib` and
  tool 3. `quality_lib` keeps its copy (the tool-6 Curry generator inlines that library). Check: the whole chain.
- **8d, tool 3's mid-night `?` rule** → `edf_io_lib.mid_recording_unscored` (tool 3 + its twin through
  `_make_tool3_curry.py`). Check: tool 3's known answers + `test_curry_twins`.
- **Not moved, by decision**: tool 2's parser, tool 4's event loaders, tool 6's load block, tool 7's event chain /
  start datetime / suffix detection: the Curry generators rewrite them. They stay "keep in sync with `edf_io_lib`"
  (CLAUDE.md / SPEC say so), guarded by `test_tool0` (tool 0 = tools 5, 6, 7 on every participant).
- Alternative (full scope): move those too, the format-specific readers passed as arguments and the four generators
  rewritten. Much larger, for code tidiness only.

### Phase 8 results (2026-10-14)

- Built by script (each removed copy checked to be the library's code; tool 1's two header functions once the
  library's optional fallback is taken out). `edf_io_lib` gains *Study configuration* (`BASE_STAGE_COLORS`,
  `CUSTOM_STAGE_PALETTE`, `load_custom_stages`, `parse_custom_field`, `custom_stage_style`), imported by the five
  libraries and tool 3. Tool 1 (Voila + batch scripts, a probe that also works in Spyder without `__file__`), tool
  1bis (`check_anonymization`), tool 3 (+ the `?` rule), `autoscoring_lib` import from it. Tools 1bis / 2 / 7 got
  verbatim copies (`detect_encoding`, the header parser, `load_event_remap`) so the sync test can guard them. Curry
  twins 2, 3, 7 regenerated (generators unchanged). A notebook cell stored as one string stays one string.
- `tests/test_edf_io_sync.py` (quick): 6 files, 25 copies; it fails on a code change, not on a docstring.
- `outside_chain.py` gains tool-1 and tool-3 runners.

Tool 0 only, display only: no computation changes, no batch tool touched.

1. **Step 5, metric figures → report only.** "All electrodes: metric distribution" (pooled box plots) and "All
   electrodes: metric by electrode" leave the tab. They go into the inspection report, inside a **closed**
   `<details>` ("Metrics per electrode: tool 6's figures"), after the tables and the two overview figures already
   kept there. They are drawn only when *Keep an inspection report* is ticked. The tab keeps the butterfly,
   median spectrogram, PSD, amplitude distribution and the averaged-metrics table, plus one line: "The metric
   figures per electrode are in the inspection report".
2. **Step 5, "One channel" → "Per-channel inspection"**, in a closed `Accordion` (hidden until step 5 has run,
   like now). Its figures are drawn **only when it is opened**, then on each channel change while it is open, so
   a step-5 run no longer draws four channel figures nobody looked at. Same content (dropdown with the first
   flagged channel pre-selected, *Show in the explorer*, metrics + histogram + PSD + time series + spectrogram).
3. **Step 1, hypnogram chosen by hand.** `dd_hypno` keeps the detected files (the `.txt` named after the
   recording, tool 7's suffix pre-selected) and gains a last option **"Other file…"**, which shows a file chooser
   (`ipyfilechooser`, opened in the recording's folder, any file anywhere). An orange note under it: the batch
   tools will not find this file (they look for `<recording><suffix>` beside the EDF): rename or copy it for
   them. Loading reads the chosen path the same way (`load_hypno_labels`, tool 7's length rule); "Other file…"
   with nothing chosen = a clear message, nothing loaded. Downstream: step 3's crosstab still compares the other
   detected files with the one in use; step 4 passes the file name where it passed the suffix (only used in tool
   5's "no hypnogram" message); the rescoring is saved beside the chosen file (as now: beside the hypnogram);
   the log says `hypnogram=<name> (chosen by hand)`.
4. **Step 2, channel table split by detected type.** One table per type, in the order EEG, EOG, EMG, ECG,
   other, each titled with its count ("EOG: 2 channels"); EEG / EOG / EMG / ECG say "none detected" when empty
   (that is the thing to see), "other" is omitted when empty. Same columns minus `detected type`. **Role vs
   type**: when a channel's role here disagrees with its detected type (used as EEG by tool 2's remap but
   detected EOG / other, or the EOG-left dropdown on a channel detected as EMG…), its row gets a `check` column
   "⚠ used as EEG" and the Checks list says so. The column only appears when such a disagreement exists.

Docs: SPEC §0 (tabs 1, 2, 5 rows, step-1 paragraph). Tests: `test_tool0.py` (drives `dd_hypno`) still passes
as is; one test added: "Other file…" loads a hypnogram picked by path. Proposed level: standard on `tool0`.

---

## 6. Risks

- **Curry generators (6, 7)** match code by string: the blocks they rewrite stay in place (§3.3), and their
  assertions fail loudly if a match is lost. `test_curry_twins.py` at each phase.
- **Floating-point noise** (tool 9, MKL): the before/after check uses the golden tolerance for floats.
- **Memory**: the whole night stays in RAM (raw + epochs + PSDs). Fine for PSG (3 to 12 EEG channels).
  Above `HD_CHANNEL_THRESHOLD` a warning says it may be slow.
- **GSSC / torch**: heavy import, done only when §3's GSSC is run.
- **Tool 0 must still execute top to bottom without data** (`test_notebooks_execute.py`).

---

## 7. Later / out of scope

- Interactive canvas (click-to-seek, zoom, mouse annotation): the TODO file is kept, trimmed to this item.
- Comparison badges against the batch outputs already on disk ("tool 7: 15 % flagged, here: 15 %").
- A Curry twin of tool 0.
- **Full phase-8 scope for the Curry twins (user, 2026-10-14, TODO)**: tools 2, 4, 6, 7 read the EDF through
  `edf_io_lib`, the format-specific readers passed as arguments, the four Curry generators rewritten; needs a user
  test of the four twins on Curry data. Written in SPEC → *Planned modules* (survives this plan).
- **Batch tools' companion suffixes (user, 2026-10-14, after the redesign)**: tools 3, 3bis, 5, 6, 7 pre-fill an
  editable `Text` with the detected suffix and list the others ("Detected: …"). Give them tool 0's step-1 pattern:
  a dropdown of the detected suffixes (with their counts) + "Other suffix…" revealing the text field. Touches five
  notebooks, their Curry twins / generators and the tests that drive `hypno_suffix` / `txt_suffix` / `csv_suffix`.

## 8. Files touched (whole task)

- **rewritten**: `tools/0_live_explore_1file_voila.ipynb`
- **new**: `tools/edf_io_lib.py`, `autoscoring_lib.py`, `macrostructure_lib.py`, `quality_lib.py`,
  `preprocessing_lib.py`, `autoreject_lib.py`, `spectral_lib.py`, `tests/refactor_check.py`,
  `tests/test_tool0.py`
- **modified** (code moved to a library + import): notebooks 1, 1bis, 2, 3, 3bis, 4, 5, 5bis, 6, 7, 8bis, 9;
  `qc_rejected_epochs_lib.py`; Curry generators 2, 4, 6, 7 (+ 8bis if needed) and their regenerated twins;
  `tests/chaincache.py`; `pytest.ini`; `SPEC.md`; `CLAUDE.md`; `TESTS.md`;
  `tools/TODO_live_explore_interactive_migration.md`
- **new output folder**: `<study>/reports_live_explore/` (inspection log + report, §1.4)
