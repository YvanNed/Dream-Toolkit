# Plan — Tool 9: Spectral features (band power + aperiodic/periodic decomposition)

**Status**: draft for review. Nothing implemented yet. Annotate freely, then confirm.

**Goal**: turn `tools/9_SpectralPower_&_AperiodicFit_PSG.py` (raw-EDF batch script) into
`tools/9_spectral_features_voila.ipynb` — an interactive Voila app that **starts from the clean epoched
data** written by tools 8 / 8bis, lets the user choose every spectral parameter, and produces
per-participant TSV + HTML plus a database-level Excel workbook and report.

This is the **first tool of the feature-extraction family** (10, 11, … will follow the same skeleton:
folder pickers → parameters → participant selection → run loop → per-file TSV/HTML → global
Excel/TSV/HTML).

**Delivery**: Voila only for now (no Jupyter twin, no batch `.py`, no Curry twin — the tool reads only
format-agnostic `.fif`, so a Curry twin would be a verbatim copy and can be produced later if wanted).
All code lives **in the notebook** (no shared `spectral_features_lib.py`).

Files touched:
- **new** `tools/9_spectral_features_voila.ipynb`
- **moved** `tools/9_SpectralPower_&_AperiodicFit_PSG.py` → `tools/old/`
- `environment.yml` (+ `openpyxl`)
- `SPEC.md` (new §9 + a *Cross-cutting* note), `CLAUDE.md` (one short reminder)

---

## 1. Section 1 — Folders, discovery, provenance

### 1.1 Pickers

| Widget | Role |
|---|---|
| `fc_data` — **Data folder** (optional) | Convenience, like tools 8/8bis: on selection it `reset(path=…)`s `fc_clean` to `<data>/derivatives`. Also the anchor for `config_param/custom_stages.json` and for the sidecar search (§1.3). |
| `fc_clean` — **Clean-epochs folder** (required) | Points at `derivatives/clean_epo_manual/`, `clean_epo_auto/`, **or** `raw_epo/`. Not one root scanned recursively — versioned runs stay apart (same rationale as 8bis). |
| `fc_subjinfo` — **Participant info** (optional) | `.csv`/`.tsv` joined into the aggregated outputs (e.g. group, night, treatment). A `Join column` dropdown is populated from its header once loaded; matching is `os.path.normcase`-wrapped on `file_id`. Non-fatal when absent/unmatched (unmatched rows keep NaN + a ⚠ count). |

### 1.2 Source auto-detection (`_clean-epo` vs `_all-epo`)

`Scan` looks for `*_clean-epo.fif` recursively; if none is found it falls back to `*_all-epo.fif` and
shows an **amber warning**: *"raw (unrejected) epochs — no manual/automatic rejection was applied"*.
In that case a checkbox **`Exclude epochs flagged by tool 7`** (ON) filters on the `.fif`
`reject_flag` metadata column. The chosen source (`clean_manual` / `clean_auto` / `raw`) is recorded in
every output as a `source` column so a database mixing both is still readable. Source label is derived
from the folder name (`clean_epo_manual` / `clean_epo_auto` / else `raw`).

### 1.3 Tool-7 sidecar (provenance) — auto-located, non-fatal

`{file_id}_preprocessing_params.json` lives in `raw_epo/`, **not** next to the clean `.fif`. Searched in
order, first hit wins:
1. beside the `.fif`;
2. `<clean_root>.parent/raw_epo*/**` (sibling of the clean folder — the standard tool-7 layout);
3. `<data_root>/derivatives/**`.

Found → read with the same logic as `qc_rejected_epochs_lib.load_params` (filters, resample, notch,
`epoch_length_s`, `psd_smoothing`, `methods_run`, thresholds) purely as **provenance** for the report and
the Excel `parameters` sheet. Not found → **⚠ warning per file** (surfaced in the scan info line and in
the run log), and `epoch_length_s` is derived from the `.fif` itself (`len(epochs.times)/sfreq`). Nothing
blocks.

### 1.4 Custom stages

`config_param/custom_stages.json` read from the data folder (same helper as tools 6/7/8) into an editable
`Custom stages` field; the stage checkboxes and every per-stage table/figure iterate
`['W','N1','N2','N3','R'] + custom_stages`.

### 1.5 Scan output

- participant list (`file_id`, subtree, n epochs, n channels, sfreq, stages present);
- **"N / M already processed"** line (skip gate = `_spectral_stage.tsv` **and** `_spectral_report.html`
  both on disk, per the toolkit's interruption-safe convention);
- **channel-coverage preview**: the union of channel names with the number of participants carrying each
  (this is where a 8/8bis channel drop becomes visible *before* running);
- warnings: missing sidecar, heterogeneous sfreq, heterogeneous epoch length, stage sets.

---

## 2. Section 2 — Parameters

Layout follows tool 7: headed sub-sections, a grey one-line description per parameter, and every box
revealed by a checkbox derives its **initial** `display` from that checkbox's value
(*checkbox-revealed-widgets* rule).

### 2a. PSD

| Widget | Default | Note |
|---|---|---|
| `dd_psd_method` | `welch` | `welch` / `multitaper` (MNE). |
| `txt_fmin` / `txt_fmax` | 0.5 / 45 Hz | Capped at `sf/2 − 0.5` per file (a ⚠ is logged when capped). |
| `txt_win_s` | 4 s | Welch window → 0.25 Hz resolution. Guard: `win_s ≤ epoch_length_s`, else clipped + ⚠. |
| `txt_overlap_pct` | 50 % | Welch only. |
| `dd_window` | `hann` | Welch only. |
| `txt_bandwidth` | 2 Hz | Multitaper only (revealed by the method dropdown). |

`average='mean'`, `n_fft = n_per_seg` — not exposed (per your answer).

### 2b. Aperiodic / periodic fit (specparam) — **per epoch × channel**

| Widget | Default | Note shown in the UI |
|---|---|---|
| `txt_fit_fmin` / `txt_fit_fmax` | 2 / 45 Hz | Independent of the PSD range; ≥ 2 Hz limits slow-wave influence (tool-7 convention). |
| `dd_aperiodic_mode` | `fixed` | **fixed** = straight line in log-log; use **knee** only when the log-log spectrum visibly *bends* (typically when fitting above ~45–50 Hz or from very low frequencies). On a 2–45 Hz sleep PSD `fixed` is almost always right; `knee` adds a parameter that is unstable when there is no real bend. Figure 4 of the individual report (log-log PSD + fit) is the check. |
| `txt_peak_width_min` / `max` | 0.5 / 20 Hz | Min must be > 2 × freq resolution (guarded). |
| `txt_min_peak_height` | 0.3 | |
| `txt_max_n_peaks` | 8 | Bounds the peak search (same rationale as tool 7: keeps noisy short-epoch fits fast). |
| `txt_peak_threshold` | 2.0 | In SD of the flattened spectrum. |
| `txt_r2_min` / `txt_mae_max` | 0.90 / 0.15 | **Quality flag only** — a row below them keeps its values with `fit_ok = False` (never dropped silently); counted per file in the report. |
| `cb_exclude_bad_fit` | **OFF** | *Exclude poorly-fitted epochs from the aggregates* — when ON, epochs with `fit_ok = False` are left out of the (stage × channel) means; the excluded count is always reported either way. OFF by default: tool 7 already flagged and rejected epochs on the very same 1/f criteria, so what reaches this tool should be fine — the checkbox is there for the residual cases. Affects **only** the aggregated tables, never the epoch-level ones. |
| `cb_parallel` + `txt_n_jobs` | ON, `-1` | joblib over epochs; exact (same numbers as serial), falls back to the serial loop on any failure. Cost reference: ~11 ms/fit → ~1.2 min for 1000 epochs × 6 ch, ~6.5 min × 32 ch. |

**No epoch subsampling** here (unlike tool 8's QC distributions): a feature table with holes is not usable.

**Smoothing of the fit copy only** (`cb_smooth_fit`, **ON**; `cb_smooth_median` ON + `txt_median_span` 3.0;
`txt_lowess_span` 2.0) — reuses tool 7's `smooth_psd_median` → `smooth_psd_lowess` (median first, then
LOWESS, on the linear PSD). The smoothed array feeds **only** specparam; **band powers are always computed
on the unsmoothed PSD**. Stated explicitly in the UI and in the report's parameter table.

### 2c. Frequency bands

An editable list of rows (`name`, `fmin`, `fmax`) with **Add band** / **Remove** buttons. Defaults:

```
delta 0.5–4 | theta 4–8 | alpha 8–12 | sigma 12–16 | beta 16–30 | gamma 30–45
```

Guards: `fmin < fmax`, band inside `[PSD fmin, PSD fmax]` (else the band is skipped for that file with a
⚠), unique names. Bands are allowed to overlap (some designs want `slow 0.5–1` inside delta) — no check.
Bin selection is `freqs >= fmin & freqs < fmax` (half-open, so adjacent bands never double-count a bin).
No JSON save/load button (dropped per your answer).

### 2d. Measures to export — one checkbox each, with a one-line description

| Column | Default | Description shown |
|---|---|---|
| `power_mean_uV2_Hz` | ON | Mean of the PSD over the band's bins — **power density** (µV²/Hz). Independent of band width. |
| `power_abs_uV2` | ON | Integral of the PSD over the band (`np.trapz`) — **absolute power** (µV²). Scales with band width. |
| `power_db` | ON | `10·log10(power_mean_uV2_Hz)` — the usual log scale for stats. |
| `power_rel` | ON | `power_abs_uV2` of the band ÷ absolute power over the whole PSD range — dimensionless, all bands sum to ~1. **Fixes the old script**, which divided a *mean* by a *sum*. |
| `power_ap_removed_db` | ON | Mean, over the band, of `PSD_dB − aperiodic_fit_dB`: how many dB the spectrum sits **above the 1/f floor**. Removes inter-subject/stage offset differences. |
| `power_ratio_over_aperiodic` | ON | Mean, over the band, of `PSD_linear / aperiodic_fit_linear` — the same idea in linear space (a fold-change: 2 = twice the 1/f floor). Not a re-expression of the previous one (mean of a ratio ≠ ratio of means), and it stays strictly positive. |
| `power_rel_periodic` | OFF | The true *relative* analogue after 1/f removal: residual power `max(PSD_lin − AP_lin, 0)` integrated over the band ÷ the same integrated over the whole range. Off by default because the clipping at 0 makes it fragile when the fit is poor. |

**Aperiodic model → dB curve** (works for both modes, extrapolated over the *whole* PSD range, including
below the fit's lower bound):

```
fixed :  ap_log10(f) = offset − exponent · log10(f)
knee  :  ap_log10(f) = offset − log10(knee + f**exponent)
ap_dB  = 10 · ap_log10        ap_linear = 10**ap_log10
```

### 2e. Stages, aggregation and coverage

- **Stage checkboxes**: built at Scan from the stages actually present in the `.fif` files (AASM + custom),
  **all ON by default** — `W` included, it is a legitimate comparison baseline; untick it if not wanted.
  A stage absent from every file simply gets no checkbox.
- **Derived groups** (checkboxes, additive rows in the aggregated tables): `NREM = N2+N3`,
  `sleep = N1+N2+N3+R`, `all`.
- **`txt_min_epochs`** (default **20**): an aggregated (stage × channel) cell computed on fewer clean
  epochs is written as **NaN** with `n_epochs` kept, so the shortage is visible rather than silent. The
  epoch-level tables are never filtered.
- **`cb_thirds`** (default OFF, *draft feature*): adds a `third` column and a companion aggregated file.
  Definition: the sleep period = first to last epoch whose stage ≠ `W`, split into three equal spans of
  `epoch_idx`; epochs outside are `pre` / `post`. **NREM–REM cycle detection is deliberately not
  attempted** in v1 (it needs the full, unrejected hypnogram and a rule set of its own) — noted as a
  follow-up.

### 2f. Averaging space (epoch → stage)

Both are exported, per your answer, because they answer different questions and the old script had an
open TODO on it:
- `*_from_lin` — average in **linear** space (µV²/Hz), then convert to dB. Preserves total power; a few
  high-power epochs dominate.
- `*_from_log` — average the per-epoch dB values. Closer to what a per-epoch statistical model sees;
  robust to single loud epochs.

The report states the difference in one sentence and the individual figures use `_from_lin` (with the
`_from_log` curve dashed behind) so the gap is visible.

---

## 3. Section 3 — Participant selection

Checkbox list (one per participant, `file_id` + subtree + n epochs), **Select all** / **None**, plus the
**`Skip already processed`** checkbox (ON) and the "N / M already done" line. A ⚠ marks any participant
whose report and data files disagree (interrupted run) — it is reprocessed.

---

## 4. Section 4 — Run

Two progress bars, like tools 6/7: participants (`i/N`) and an inner pipeline bar with a 5-segment legend
— **Load / PSD / 1-f fit / Bands / Report** — advanced from a `set_phase()` helper, the 1/f segment
animated per epoch batch.

Per participant (each step in its own `try/except`; fatal → `failed` list + `continue`, non-fatal → ⚠):

1. `mne.read_epochs` (preload), read `metadata` (`epoch_idx`, `stage`, and the tool-8/8bis decision
   columns when present); apply the `raw`-source filter if relevant.
2. Welch/multitaper PSD → `psds_uV2` (n_ep, n_ch, n_freqs), µV²/Hz.
3. Smoothed copy → specparam per (epoch × channel), joblib-parallel: `offset`, `exponent`, `knee`
   (knee mode), `r2`, `mae`, `fit_ok`, `n_peaks` + one row per peak (CF, PW, BW).
4. Aperiodic curve extrapolated to the full PSD range → the aperiodic-removed measures.
5. Band powers (unsmoothed PSD) for every selected measure.
6. Aggregation to (stage × channel × band) — and (stage × channel) for the aperiodic — with `n_epochs`,
   `mean`, `median`, `sd`, in both averaging spaces; `< min_epochs` → NaN.
7. Write the per-file TSVs (**data before report**), then the HTML report.

Then, once over all files: rebuild every global table by **globbing the per-file TSVs on disk**
(interruption-safe — never from an in-memory list), write the Excel workbook and the database report.

---

## 5. Outputs

### 5.1 Per participant — data → `derivatives/features_spectral/<subtree>/`

| File | Content |
|---|---|
| `{file_id}_bandpower_epoch.tsv` | Long format: `file_id, source, epoch_idx, stage, third?, channel, band, <selected measures>`. |
| `{file_id}_aperiodic_epoch.tsv` | `file_id, epoch_idx, stage, third?, channel, offset, exponent, knee, r2, mae, fit_ok, n_peaks`. |
| `{file_id}_periodic_peaks.tsv` | One row per specparam peak: `…, peak_idx, n_peaks, center_freq_hz, peak_power, bandwidth_hz`. |
| `{file_id}_spectral_stage.tsv` | **Master table for stats** — `file_id, stage, channel, band, n_epochs, n_fit_ok`, then `<measure>_mean_from_lin/_from_log`, `_median_…`, `_sd_…`. Also the **skip-gate data marker**. |
| `{file_id}_aperiodic_stage.tsv` | Aggregated aperiodic: `file_id, stage, channel, n_epochs, n_fit_ok, offset/exponent/r2/mae × mean/median/sd`. |
| `{file_id}_psd_stage.tsv` | Mean PSD per `stage × channel × freq` (`psd_uV2_Hz_from_lin`, `psd_db_from_lin`, `psd_db_from_log`, `psd_ap_removed_db`, `n_epochs`) — lets the group figures be redrawn without recomputing. |
| `{file_id}_spectral_stage_third.tsv` | Only when `cb_thirds` is ON: same schema as `_spectral_stage.tsv` + a `third` column. |

### 5.2 Per participant — report → `reports_features_spectral/<subtree>/{file_id}_spectral_report.html`

An `mne.Report` with (per your C4 list):
1. **Hypnogram of the kept epochs** — stage vs `epoch_idx`, with the epochs removed by tool 7/8/8bis
   drawn as gaps + a "clean epochs per stage" bar.
2. **Mean PSD per stage**, one subplot per channel (dB/Hz ± SEM across epochs).
3. **Mean PSD per stage after removing the aperiodic component**, same layout.
4. **log-log PSD + aperiodic fit per stage** (fit superposed, `R²`/exponent annotated) — the visual check
   for `fixed` vs `knee`.
5. **Exponent & offset across the night** — one point per epoch, coloured by stage, one subplot per
   channel (+ per-stage mean lines).
6. **Band power per stage** — boxplots over epochs, grouped by band, one subplot per channel.
7. **Parameters + coverage tables** — every widget value used, the sidecar provenance (filters, epoch
   length, tool-7 PSD smoothing), `n_epochs`/`n_fit_ok` per stage × channel, and the `min_epochs`
   NaN count.

### 5.3 Database level — `reports_features_spectral/` root

| File | Content |
|---|---|
| `global_spectral_stage.tsv` | Concatenation of every `_spectral_stage.tsv` found by `rglob` (+ the participant-info join). |
| `global_aperiodic_stage.tsv` | Same for `_aperiodic_stage.tsv`. |
| `global_psd_stage.tsv` | Same for `_psd_stage.tsv` (feeds the group figures). |
| **`spectral_features_database.xlsx`** | Multi-sheet workbook (needs `openpyxl`): `bandpower_stage`, `aperiodic_stage`, `psd_stage`, `coverage`, `parameters`, `failed`. Built from the globbed TSVs — **the TSVs stay the source of truth**; a sheet exceeding Excel's ~1 048 576-row limit is truncated with a bold warning row and a pointer to the TSV. |
| `spectral_database_report.html` | Group figures: mean PSD per stage across participants (± SEM between subjects, one subplot per channel), band-power boxplots per stage, exponent/offset distributions per stage, clean-epoch count per participant, channel-coverage table. |
| `spectral_features_failed.tsv` | One row per failed participant + reason. |

### 5.4 Heterogeneous channels (your A3)

Confirmed: at the **epoch level** a channel dropped by tool 8/8bis simply has no rows (it does not exist
in that `.fif`). At the **aggregated/group level** the tables are reindexed on the **union** of channels
across participants, so a participant missing `C4` gets `C4` rows with `n_epochs = 0` and **NaN** values —
nothing crashes, and the `coverage` sheet says how many participants actually carry each channel.

---

## 6. Conventions inherited (no deviation)

- **Skip + cumulative merge**, data written before report, globals rebuilt by globbing the disk.
- `os.path.normcase` on **both** sides of every path/id comparison.
- Per-item `try/except` with fatal vs non-fatal; every button callback wrapped; errors surfaced in a widget.
- All user-facing strings in **English**.
- Data → `derivatives/`, reports → `reports_*` (split precomputed at scan time).
- Subtree mirroring: `<subtree>` = the participant's path relative to the chosen clean-epochs folder.

---

## 7. Doc & environment updates (proposed, to confirm at the end)

- `environment.yml`: add `openpyxl` (conda-forge) — the only new dependency.
- `SPEC.md`: rewrite **§9** for the new tool (inputs, parameters, outputs, figures), add the new tool to
  the *Repository structure* + *How to run* lists, and add one *Cross-cutting* bullet: **"Feature tools
  start from `clean_epo_*`; the tool-7 sidecar is provenance only and its absence is never fatal."**
- `CLAUDE.md`: one short reminder + pointer (feature-tool family, clean-epo input, smoothing applies to
  the fit copy only, band powers always on the raw PSD).
- Delete this plan file once the task is confirmed done.

---

## 8. Resolved decisions (answered on review)

1. **Stages** — every stage present in the `.fif` gets a checkbox, **all ON by default**, `W` included (§2e).
2. **`power_rel_periodic`** — kept as an opt-in checkbox, **default OFF** (§2d).
3. **Stage-level specparam fit** (fitting the per-stage *average* PSD for cleaner periodic peaks) —
   **left out of v1**. Noted here as the obvious follow-up if the per-epoch peaks turn out too noisy.
4. **`fit_ok = False` epochs** — kept in the aggregates by default, with a `cb_exclude_bad_fit` checkbox
   (**OFF**) to drop them; the count is reported either way (§2b). Tool 7 already rejects on the same 1/f
   criteria, so this should stay a rare case.

## 8bis. Deltas decided during implementation (all validated on the test data)

Small design calls made while building, each flagged here because they refine what §2-§5 said:

1. **Band integral = rectangular sum (`sum(PSD) x df`), not `np.trapz`.** Over the band's own bins a
   trapezoidal rule drops the half-bin at each edge, so the six bands summed to **0.971** instead of 1.
   A rectangular sum gives each band its exact width and the relative powers now sum to **0.99997**
   (the residual is the single 45 Hz bin excluded by the half-open gamma band).
2. **No 1/f-fit subsampling, and parallelism is auto-disabled on small jobs.** Measured 2.8 ms/fit;
   spawning joblib workers costs ~3 s on Windows, so below **5000 (epoch x channel) fits** the serial
   loop wins and the parallel path is skipped (`PARALLEL_MIN_FITS`). Serial and parallel were verified
   to return **bit-identical** offsets, exponents, R2, MAE and peaks.
3. **`cb_do_fit` master checkbox** (ON): unticking skips specparam entirely and disables the three
   1/f-based measures - the main speed lever, mirroring tool 7's `cb_1f`.
4. **Bands partly outside the PSD range**: fewer than 2 frequency bins -> the band is skipped with a
   warning; partly covered -> kept with a warning that its power is not comparable to a fully covered
   band. (Found by setting `fmax = 30` with a 30-45 Hz gamma band: it used to be computed from a
   single bin.)
5. **`file_id` forced to `str`** when the per-file tables are re-read: a numeric id (`73`) came back as
   an int, giving a mixed int/str column that breaks matching and sorting.
6. **Database outputs are non-fatal**: a failure of the workbook or the database report no longer
   swallows the run summary (found via a real Windows 260-character path overflow). The scan now also
   warns when the longest output path exceeds 245 characters - relevant on this machine's OneDrive tree.
7. **Channel-union padding applies to `bandpower_stage` and `aperiodic_stage` only**, not to
   `psd_stage` (one row per frequency bin would multiply the padding by ~180).
8. **`knee` mode degeneracy warning**: on real sleep data with no bend, specparam returns
   `knee ~ 5e28`; the run log now says so explicitly and points at `fixed`.

## 8ter. Changes from the first review round (after your test run)

1. **Widget disappearance bug — fixed.** `_IND_A` / `_IND_B` were single `Layout` objects **shared**
   by several boxes; ipywidgets keeps the reference, so setting `display='none'` on one hid every
   widget sharing it (switching the PSD method twice, or unticking the fit, wiped the section). Both
   are now **factory functions** (`_ind(px)` / `_w(px)`) returning a fresh `Layout` per widget — the
   exact rule CLAUDE.md already states for tool 7 ("the shared `Layout` objects are never mutated").
   Verified by replaying your sequence: at every step only the intended box is hidden, and every
   `Layout` object is unique.
2. **Pre-fit smoothing tightened** to two rows, tool-7 style: the LOWESS span sits on the same row as
   its master checkbox, the median span next to its own; both revealed by the master. Same for the
   parallel workers field.
3. **Fold-change over the 1/f fit** now **OFF** by default; the *relative periodic power* description
   states that it rests on the specparam fit and says to check the `fit_ok` rate first.
4. **Derived groups (NREM / sleep / all) removed** — one groupby away in the exported tables, and they
   doubled the table size.
5. **Data-table previews removed from the individual report.** It now carries only the parameters used
   and the coverage table, pointing at `{file_id}_spectral_stage.tsv` for the numbers.
6. **Averaging space is now a choice** (§2e): `log space` **ON** by default, `linear space` off; at
   least one required. Only the selected suffixes are exported (`_from_log` / `_from_lin`), which
   halves the column count in the default setup. It also fixed a **real bug**: the PSD figure's title
   promised a dashed linear-space curve that was never drawn — both curves now appear when both spaces
   are selected (measured gap on the test data: 0.3 dB in delta, 1.4 dB in gamma).
7. **"Where to look" section** closing the database report: a question → file → column table over the
   seven kinds of table produced, naming the exact column for the selected averaging space, plus the
   caveats to check before reading a number (`n_epochs`, the `n_epochs = 0` padding rows, the 1/f
   extrapolation under delta).

## 8quater. Changes from the second review round

1. **Smoothing layout**: the two spans now sit on their own lines, one indent deeper than the master
   checkbox; `Running-median pre-step` renamed to `Median span (Hz)` with its value field beside it,
   symmetric with the LOWESS line.
2. **Why some spectra still fail the quality bounds** (your question) is now answered *in the report*,
   next to the `fit_ok` rate, and the new `fit_settings_diff()` names the 1/f settings that differ
   from the tool-7 run when its sidecar is available (on the test data: *"R2 bound 0.95 in tool 7,
   0.9 here"*). Two independent reasons: tool 7 flags **(epoch, channel) pairs** while 8/8bis decide
   **per epoch**, and the fit is **recomputed here** with this notebook's own settings.
3. **`_sd` columns added to `_psd_stage.tsv`** (`psd_uV2_Hz_sd`, `psd_db_from_*_sd`,
   `psd_ap_removed_db_sd`). This was the one thing that could **not** be redrawn from the tables: the
   +/- SEM band. With `sem = sd / sqrt(n_epochs)` every report figure is now reproducible from the
   TSVs. Cost measured at 0.01 s per participant (0.1 % of the pipeline).
4. **Report guidance extended** with two sections: *Redrawing a figure yourself* (which table and
   which columns per figure, including the analytical aperiodic curve) and *Grouping participants
   differently* (a 5-line pandas snippet; it states whether the participant-info table was joined).
5. **Report ordering**: `save_report_html` now takes an ORDERED list of `('fig'|'html', title,
   payload)` items instead of "all figures then all text", so *Clean epochs per participant* sits
   below *Summary*, with *Where to look* closing the database report.
6. **Closing markdown cell + an end-of-run banner** pointing at the database report as the first
   thing to open.

## 8quinquies. Third review round — the smoothing sub-section

`LOWESS` gained its **own checkbox**, so each filter can now be used alone (median only, LOWESS only,
both, or neither); `smooth_for_fit` honours both flags in the oscip order (median first). Both labels
end with `:` and both checkboxes share a 200 px width, so the two value fields line up vertically.
Unticking **both** filters sets `enabled = False`, so the report never claims a smoothing that did not
run. Verified on a spectrum with an injected spike: raw 9.62 → median+LOWESS 1.66, median alone 1.67,
LOWESS alone 3.34 — i.e. the median pre-step is what actually removes narrow spikes, as documented.

## 8sexies. Bug found by comparing against `Thomas4Carl_extract_features.py`

**`normalization='full'` is now passed explicitly to the multitaper PSD.** MNE's
`compute_psd(method='multitaper')` defaults to `normalization='length'`, which is **not** a power
density: it is larger by a factor of `sfreq` (+23.5 dB measured at 256 Hz) and the factor follows the
sampling rate, so recordings at 256 and 512 Hz would not have been comparable. Welch was never
affected (MNE returns a proper density there) and it is the default, so no run made so far is wrong —
but the multitaper option was.

After the fix, Welch and multitaper agree to within 0.2 dB on every band above 4 Hz (delta differs by
0.7 dB, a genuine resolution effect: 0.033 Hz vs 0.25 Hz on a very steep slope) and the aperiodic
exponent matches to 0.01 (2.46 vs 2.47 on N3/C3).

*Related, in the other script only*: `np.log10(psd + 1e-10)` on a PSD in V²/Hz is a **100 µV²/Hz
floor**, above most of the sleep spectrum; it flattens the 1/f slope (−0.47 exponent) and distorts
band ratios. Tool 9 has no additive epsilon anywhere — every log is guarded by
`np.where(psd > 0, psd, np.nan)`, so a null bin propagates as NaN instead of an invented value.

## 9. Implementation order

1. Notebook skeleton + Section 1 (pickers, scan, source detection, sidecar lookup, coverage preview).
2. Section 2 (all parameter widgets + guards + the revealed-box rule).
3. Section 3 (participant checkboxes, skip gate).
4. Section 4 core: load → PSD → specparam → aperiodic removal → band powers → per-file TSVs.
5. Individual HTML report (7 figures).
6. Global tables + Excel workbook + database report.
7. Move the old script to `tools/old/`, update `environment.yml`, propose the SPEC.md / CLAUDE.md edits.

Validated at each step on a real dataset (`tools/test_data/` + a real `clean_epo_*` folder if available).

---

## 10. Proposed doc edits (paste-ready — awaiting your confirmation)

`environment.yml` is already updated (`openpyxl=3.1.*`, the only new dependency). The edits below are
**proposed, not applied**, per the CLAUDE.md rule.

### 10.1 SPEC.md — *Repository structure*, replace the tool-9 line

    │   ├── 9_spectral_features_voila.ipynb        # Spectral features: band power + aperiodic fit (Voila GUI)

and note in the `old/` line that `9_SpectralPower_&_AperiodicFit_PSG.py` (the raw-EDF batch ancestor)
now lives there.

### 10.2 SPEC.md — *Conda environment*, add

    - **openpyxl 3.1** — writes tool 9's `.xlsx` workbook (via `pandas.ExcelWriter`)

### 10.3 SPEC.md — *How to run → Interactive Voila apps*, add

    voila tools/9_spectral_features_voila.ipynb

### 10.4 SPEC.md — *Cross-cutting procedures*, new bullet

- **Feature tools start from the clean epochs (`clean_epo_*`)**: every feature-extraction tool
  (9, then 10, 11 …) reads `{file_id}_clean-epo.fif` from `derivatives/clean_epo_manual/` or
  `clean_epo_auto/`, chosen through **one explicit folder picker** (never a recursively-scanned root —
  versioned runs must stay apart, as in tools 8/8bis). When that folder holds no `*_clean-epo.fif` the
  tool **falls back to tool 7's `*_all-epo.fif`** with an amber warning and an option to drop the epochs
  tool 7 flagged; the origin is recorded in a `source` column of every output (`clean_manual` /
  `clean_auto` / `raw`). Tool 7's `{file_id}_preprocessing_params.json` lives in `raw_epo/`, not beside
  the clean `.fif`: it is **auto-located** (beside the fif, then under the clean folder's parent, then
  under the data folder), used for **provenance only**, and its absence is **never fatal** (the epoch
  length then comes from the `.fif` itself). Outputs follow the data/reports split: data →
  `derivatives/features_<name>/<subtree>/`, reports → `reports_features_<name>/<subtree>/`. Because
  tools 8/8bis drop channels per participant, the database-level tables are **padded to the union of
  channels** (`n_epochs = 0` + NaN) so a channel missing for one participant is visible rather than
  silently absent.

### 10.5 SPEC.md — replace the whole of *§9. Spectral Analysis* with

### 9. Spectral features — band power + aperiodic/periodic decomposition (`9_spectral_features_voila.ipynb`)

First tool of the **feature-extraction family** (10, 11 … will follow the same skeleton). Voila only for
now (no Jupyter twin, no batch `.py`, no Curry twin — it reads only format-agnostic `.fif`, so a Curry
twin would be a verbatim copy). All code lives in the notebook. Replaces the raw-EDF batch script
`9_SpectralPower_&_AperiodicFit_PSG.py`, now archived in `tools/old/`; unlike its ancestor it **starts
from the clean epochs** (tools 8/8bis) and never re-does the preprocessing.

**Section 1 — Folders and scan**: an optional **Data folder** (pre-points the others, supplies
`config_param/custom_stages.json`), the required **Clean-epochs folder**, and an optional **Participant
info** `.csv`/`.tsv` joined on a chosen column (`os.path.normcase` on both sides). `Scan` reads each
`.fif` header + metadata **without loading the signal**, shows the participant table (epochs, channels,
sfreq, epoch length, stages, sidecar found or not), the **N / M already processed** count, the
**channel-coverage table**, and warns about: a missing tool-7 sidecar, mixed sampling rates or epoch
lengths, a data/report mismatch (interrupted run → reprocessed), and an output path close to the
**260-character Windows limit**. See the cross-cutting bullet above for the source detection, the
sidecar lookup and the output layout.

**Section 2 — Parameters** (all echoed in the report and in the workbook's `parameters` sheet):
- **PSD**: `welch` (default) or `multitaper`; `fmin` 0.5 / `fmax` 45 Hz (capped at Nyquist per file);
  Welch window 4 s (→ 0.25 Hz resolution), 50 % overlap, `hann` (window clipped to the epoch length when
  shorter); multitaper `bandwidth` 2 Hz. `average='mean'` is not exposed.
- **Aperiodic/periodic fit (specparam), per epoch × channel**, behind a master checkbox `cb_do_fit`
  (ON; unticking skips specparam entirely and disables the three 1/f-based measures): fit range 2–45 Hz
  **independent of the PSD range**, `aperiodic_mode` `fixed` (default) or `knee`, `peak_width_limits`
  [0.5, 20], `min_peak_height` 0.3, `max_n_peaks` 8, `peak_threshold` 2.0. Quality bounds `R² ≥ 0.90` /
  `MAE ≤ 0.15` set `fit_ok` — poorly-fitted epochs are **kept** and counted, with an opt-in
  `cb_exclude_bad_fit` (OFF) to drop them from the aggregates only. **Smoothing (median 3 Hz then
  LOWESS 2 Hz, tool 7's helpers, ON) is applied to the copy of the spectrum that feeds the fit ONLY —
  the band powers always use the raw PSD.** `cb_parallel` (joblib, exact — verified bit-identical to the
  serial loop) is auto-skipped below `PARALLEL_MIN_FITS` (5000 fits), where worker start-up costs more
  than it saves. **No epoch subsampling** (unlike tool 8's QC distributions: a feature table with holes
  is unusable).
- **Bands**: editable rows (name, fmin, fmax) + *Add band*; defaults `delta 0.5–4 / theta 4–8 /
  alpha 8–12 / sigma 12–16 / beta 16–30 / gamma 30–45`. Half-open `[fmin, fmax)`, so adjacent bands
  never share a bin and tiling bands give relative powers summing to 1. A band with < 2 bins inside the
  PSD range is skipped (warning); a partly covered band is kept with a warning.
- **Measures** (checkbox + one-line description each): `power_mean_uV2_Hz`, `power_abs_uV2`,
  `power_db`, `power_rel`, `power_ap_removed_db`, `power_ratio_over_aperiodic` (all ON) and
  `power_rel_periodic` (OFF). **`power_abs_uV2` is a rectangular sum `Σ PSD × df`**, not a trapezoidal
  rule over the band's own bins — the latter drops the half-bin at each edge and made the six bands sum
  to 0.971 instead of 1. `power_rel` = band absolute power ÷ absolute power of the whole PSD range
  (**fixing the old script**, which divided a *mean* by a *sum*). `power_ap_removed_db` and
  `power_ratio_over_aperiodic` are **not** re-expressions of each other (the mean of a ratio is not the
  ratio of the means). The aperiodic model is evaluated **analytically over the whole PSD axis**
  (`fixed`: `offset − exponent·log10(f)`; `knee`: `offset − log10(knee + f^exponent)`), so it is
  extrapolated below the fit's lower bound — which is why delta sits *below* the fit and its
  1/f-corrected value must be read with care.
- **Stages and aggregation**: one checkbox per stage actually present in the files (AASM + custom, all
  ON); derived groups `NREM = N2+N3`, `sleep = N1+N2+N3+R`, `all`; `Min epochs per cell` (default 20 →
  the cell is NaN but `n_epochs` is kept, so a shortage is visible); optional **night thirds** (draft,
  OFF): the sleep period (first to last non-W epoch) split into three equal spans of the **original**
  epoch index, written to a companion table — NREM–REM cycle detection is deliberately not attempted.
- **Two averaging spaces** (the old script's open TODO): the dB measures are aggregated **twice** —
  `_from_log` (mean of the per-epoch dB, what a per-epoch model sees, robust to one loud epoch) and
  `_from_lin` (dB of the linear mean, preserves total power but a few high-power epochs dominate). The
  linear measures get a plain mean/median/sd. The individual PSD figure draws the `_from_log` curve
  solid and the `_from_lin` one dashed, so the gap is visible.

**Section 3 — Participants**: one checkbox each (+ *Select all* / *None*), `Skip already processed`
(ON). Skip gate = **both** `{file_id}_spectral_stage.tsv` and `{file_id}_spectral_report.html` on disk.

**Section 4 — Run**: two progress bars (participants + a 5-segment pipeline bar *Load / PSD / 1-f fit /
Bands / Report*). Per participant, each step in its own `try/except`; **data is written before the
report**, with `_spectral_stage.tsv` written **last** as the skip marker. The database tables are then
rebuilt by **globbing the per-file tables from disk**; the workbook and the database report are
**non-fatal** (a locked file or an over-long path leaves every per-participant output intact and still
prints the run summary).

**Outputs — data → `derivatives/features_spectral/<subtree>/`**: `{file_id}_bandpower_epoch.tsv` (one
row per epoch × channel × band), `{file_id}_aperiodic_epoch.tsv` (offset, knee, exponent, R², MAE,
`fit_ok`, `n_peaks`), `{file_id}_periodic_peaks.tsv` (one row per specparam peak: centre frequency,
power, bandwidth), **`{file_id}_spectral_stage.tsv`** (the master table for statistics: stage × channel
× band with `n_epochs`, `n_fit_ok`, mean/median/sd of every measure), `{file_id}_aperiodic_stage.tsv`,
`{file_id}_psd_stage.tsv` (mean PSD per stage × channel × frequency, so the figures can be redrawn
without recomputing), and `{file_id}_spectral_stage_third.tsv` when the thirds option is on.

**Outputs — reports → `reports_features_spectral/<subtree>/{file_id}_spectral_report.html`**: clean
epochs against the original epoch index (removed epochs appear as gaps) + count per stage; mean PSD per
stage per channel (± SEM); the same after removing the aperiodic component; **log-log PSD with the
aperiodic fit dashed** (the visual check of `fixed` vs `knee`); **exponent and offset across the
night**, one point per epoch coloured by stage (the sleep-cycle structure is plainly visible); band
power boxplots per stage; and the parameters + coverage tables.

**Database level** (at the `reports_features_spectral/` root): `global_spectral_stage.tsv`,
`global_aperiodic_stage.tsv`, `global_psd_stage.tsv` (all rebuilt by globbing, with the participant info
joined and — for the first two — **padded to the channel union**),
**`spectral_features_database.xlsx`** (sheets `bandpower_stage`, `aperiodic_stage`, `psd_stage`,
`coverage`, `parameters`, `failed`; a sheet over Excel's row limit is truncated with a pointer to the
TSV, which stays the source of truth), `spectral_database_report.html` (group PSD per stage ± SEM
between subjects, the same after 1/f removal, aperiodic parameters per stage, band power across
participants, clean epochs per participant) and `spectral_features_failed.tsv`.

**Validated** on `tools/test_data` (8_N1, 73, 100 — 702 / 625 / 359 clean epochs, 3 / 2 / 4 channels to
exercise the heterogeneous-channel path): 3 participants in ~34 s, exponent N3 2.72 > N2 2.29 > R 2.25 >
W 0.80 and a sigma fold-change of ~6 over the 1/f floor in N2, i.e. the expected physiology.

### 10.6 CLAUDE.md — new short reminder under *Per-tool*

- **Spectral features (`9_spectral_features_voila`)**: first **feature tool** — starts from
  `clean_epo_*` (tool 8/8bis), never from the raw; tool 7's params sidecar is **provenance only**
  (auto-located, absence non-fatal). Two invariants: the **PSD smoothing feeds the specparam fit only**
  (band powers always use the raw PSD), and **band integrals are rectangular sums** (`Σ PSD × df`, so
  tiling bands give relative powers summing to 1). Outputs split data → `derivatives/features_spectral/`,
  reports → `reports_features_spectral/`; database tables are globbed from disk and **padded to the
  channel union** because 8/8bis drop channels per participant. → SPEC §9.
