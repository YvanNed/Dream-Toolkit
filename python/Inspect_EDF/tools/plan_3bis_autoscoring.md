# Plan: tool 3bis — automatic sleep scoring with GSSC (`3bis_autoscoring_voila.ipynb`)

## Context
Some recordings have no expert hypnogram, and where one exists the user wants to compare it with an
automatic scoring and get the per-epoch stage probabilities (hypnodensity). First version = **simple**:
GSSC native consensus hypnogram + hypnodensity + confusion matrix / agreement vs the expert hypnogram.
More outputs will come in later updates (keep the code easy to extend).

Decisions taken with the user:
- Channels come from **tool 2's `remap_reref_persubject.json`** (EEG = remapped channels, re-referenced per
  `ref_channels`; EOG = `context_channels.eog_left/eog_right`), pre-filled but editable.
- Consensus = **GSSC native** `loudest_vote`: `mne_infer` scores every combination (each chosen EEG or
  none × each EOG or none, minus "nothing"), and **per epoch keeps the most confident combination**
  (lowest entropy), not a majority vote nor a mean. The hypnodensity = that winning combination's
  softmax, so stage and probabilities follow the same rule (unlike Thomas's script, which averages the
  probabilities: not reproduced). Other consensus rules (mean of probabilities, majority vote) = future
  update (needs GSSC's private API).
- **Scoring mode** chosen by the user: **one electrode (fast)** = 1 EEG + EOG (3 combinations), or
  **several electrodes (consensus)** = N EEG + EOG (≈ 2N+1 combinations, time grows accordingly).
  Default: one electrode, **C3** (C4 if C3 absent).
- Hypnogram written **beside the EDF** (tool-3 format, usable by tools 5–9 through the suffix) + data in
  `derivatives/`, reports in `reports_*/`.
- **Voila only** (like tools 5–9), no Jupyter twin, no Curry twin, no batch `.py` for now.

## Dependencies (answer to the user's question)
- **gssc is pip-only** (not on conda-forge; latest 0.0.9, 2023, deps: mne, torch, pandas, matplotlib).
- **torch is cleaner from conda-forge** (the env is all conda-forge): you install it yourself, I update
  `environment.yml`:
  ```powershell
  # from a miniforge prompt
  conda install -n inspect_edf -c conda-forge pytorch-cpu
  & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" -m pip install gssc==0.0.9
  ```
  pip then sees torch already installed and does not reinstall it. If conda hits the ICM SSL error (as
  for pytest), fallback: `python -m pip install torch --index-url https://download.pytorch.org/whl/cpu`.
- `environment.yml`: add `- pytorch-cpu` (conda-forge deps) and `- gssc==0.0.9` under `pip:`.
- Known trap (seen in both example scripts): torch ≥ 2.6 defaults `torch.load(weights_only=True)`, which
  refuses GSSC's packaged weights → set `os.environ['TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD'] = '1'` before
  `EEGInfer()` (Nico's approach; simpler than monkey-patching `torch.load`), with a comment.

## Notebook layout (`tools/3bis_autoscoring_voila.ipynb`)
Pattern copied from tool 7 (folder chooser + hypnogram-suffix detection + load/re-reference block) and
the `ParticipantSelector` from `tools/participant_selection_lib.py` (called with `excluded=None`: like
tool 5, tool 3bis knows nothing of exclusions).

**Imports cell**: gssc/torch imported in `try/except`, a red box with the install commands if missing
(the rest of the notebook still renders).

**Section 1 — Data folder & inputs**
- `FileChooser` data folder → recursive `*.edf` scan (as tools 6/7).
- Tool-2 JSON auto-located at `<data>/config_param/remap_reref_persubject.json` (lenient loader). Info line
  `N / M EDFs configured in tool 2`, `K with EOG declared (context_channels)`. EDFs absent from the JSON
  are listed and not run.
- `Expert hypnogram suffix` (auto-detected, tools 6/7 rule: longest among ≥ 50 %, `event` excluded; the
  GSSC output suffix is also excluded from the candidates). Empty / missing file → scoring only, no
  comparison (non-fatal, stated in the report).
- `Output hypnogram suffix` text, default `_Hypnogram_gssc.txt`.
- `ParticipantSelector` + `Skip already processed` + "N / M already done" line.

**Section 2 — Channels**
- `Scoring mode` `ToggleButtons`: **`One electrode (fast)`** (default) / `Several electrodes (consensus)`,
  with a one-line explanation of each and of the loudest-vote rule.
- EEG choice over the **harmonized** (post-remap) EEG names minus the reference channels (common across
  configs, so one global choice is enough). One-electrode mode → `Dropdown`, default `C3`, else `C4`, else
  the first central/any channel. Consensus mode → `SelectMultiple`, default = every EEG channel. The
  displayed widget follows the mode (initial `display` derived from the mode value, per the
  checkbox-revealed-widgets rule).
- Participant lacking a chosen electrode: one-electrode mode → ⚠ listed, **not run** (no silent
  substitution); consensus mode → run on the chosen channels present (⚠ if some are missing, not run if
  none). Channels actually used recorded per participant (`eeg_used`).
- EOG mode `Dropdown`: **`L − R bipolar` (default)** / `L and R separately` / `none`. Participant without
  `context_channels` → EEG only, ⚠ in the log and report.
- A short note: GSSC filters 0.3–30 Hz and resamples itself; inference on CPU takes ~1–3 min/night.

**Section 3 — Run** (per participant, wrapped in `try/except`, fatal → `failed` list + `continue`)
1. `mne.io.read_raw_edf(..., include=remap keys + EOG raw names, preload=True)` →
   `drop_suffix_duplicates` → `adapt_remap_dict_to_suffixes` → rename → channel types → re-reference per
   `ref_channels` (copy of tool 7's block: list / `'average'`). EOG L−R built as a new channel
   (`add_channels` of a `RawArray`).
2. `EEGInfer(use_cuda=torch.cuda.is_available())`, `mne_infer(raw, eeg=[...], eog=[...], eeg_drop=True,
   eog_drop=True, filter=True)`. Accept 3 return values (`stages, times, probs`); if only 2 are returned
   (older gssc) → clear error naming the required version (no use of GSSC's private API). Rows of `probs`
   checked to sum to 1 (softmax applied otherwise, with a comment).
3. Alignment: GSSC cuts 30 s epochs from t = 0 and drops the tail (`cut='back'`) → epoch *i* = expert line
   *i* (same convention as tool 7, `#N = index + 1`). Compare on `min(n_gssc, n_expert)`; ⚠ if they differ
   by > 1 epoch.
4. Comparison on epochs whose expert label ∈ {W, N1, N2, N3, R} (MT / custom / `?` excluded and counted).
   `sklearn`: accuracy, Cohen's κ, per-stage F1, confusion matrix (counts + row %).
5. Write data **before** the report (interruption-safe), then the report.

**Section 4 — Database summary** (rebuilt by globbing per-file TSVs on disk, never from memory).

## Outputs
Beside the EDF:
- `{id}_Hypnogram_gssc.txt`: one label per line (`W/N1/N2/N3/R`), exactly the tool-3 format.

`<data>/derivatives/autoscoring_gssc/<subtree>/`:
- `{id}_hypnodensity.tsv`: `epoch_idx, onset_s, stage_gssc, p_W, p_N1, p_N2, p_N3, p_R, confidence`
  (max prob), `stage_expert` (empty if none), `agree`.
- `{id}_scoring_params.json`: gssc/torch versions, scoring mode, consensus rule (`loudest_vote`),
  n combinations, EEG/EOG channels used, EOG mode, reference, date.

`<data>/reports_autoscoring_gssc/<subtree>/`:
- `{id}_autoscoring_comparison.tsv`: one row: `file_id, n_epochs_gssc, n_epochs_expert, n_compared,
  n_excluded_expert, accuracy, kappa, f1_W … f1_R, scoring_mode, eeg_used, eog_used, status`.
- `{id}_confusion.tsv`: long format `expert, gssc, n`.
- `{id}_autoscoring_report.html` (base64 PNGs): hypnograms (expert black / GSSC red, stacked), hypnodensity
  (stacked area, fixed stage colours), confusion matrix (counts annotated, row-% colour), metrics table.

`<data>/reports_autoscoring_gssc/` root:
- `global_autoscoring_summary.tsv` (glob of the per-file comparison TSVs).
- `global_confusion.tsv` (sum of every `_confusion.tsv`) + `global_autoscoring_report.html`: pooled
  confusion matrix, κ / accuracy per participant (strip plot), table.
- `failed_autoscoring.tsv` if any.

Skip gate: `_Hypnogram_gssc.txt` **and** `_hypnodensity.tsv` **and** report present (mismatch → ⚠ + rerun).
All path comparisons through `os.path.normcase`.

## Files touched
- **New** `tools/3bis_autoscoring_voila.ipynb`.
- `environment.yml` (pytorch-cpu + gssc).
- `SPEC.md` (proposed for confirmation): repository tree, *How to run* list, new §3bis section (inputs,
  channels, outputs, alignment, comparison rules), `Conda environment` (torch, gssc).
- `CLAUDE.md`: one short per-tool reminder (GSSC native consensus, `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD`,
  hypnogram beside the EDF, Voila only).
- `TESTS.md` + `tests/`: the notebook is picked up automatically by `test_notebooks_execute.py` (quick).
  Add a `tool3bis` marker to `pytest.ini`. No golden/chain change (3bis is outside the 5→9 chain).
- During implementation: copy this plan to `tools/plan_3bis_autoscoring.md` (project rule), deleted at the end.

## Reuse
- Tool 7 (`tools/7_preprocessing_voila.ipynb`): `_detect_hypno_suffixes`, EDF load/`drop_suffix_duplicates`/
  `adapt_remap_dict_to_suffixes`, re-reference block (~cell lines 2304–2381), lenient JSON loader.
- `tools/participant_selection_lib.py`: `ParticipantSelector`, `id_key`.
- Example scripts `tools/autoscoring_example_scripts/Thomas4Carl_…_gssc.py` (EOG L−R, figure layout,
  confusion/κ/F1) and `Nico_NapNest_GSSC.py` (3-value `mne_infer`, env var).

## Verification
1. After you install torch + gssc: `python -c "import gssc, torch"` in the env.
2. Run the Voila on `tools/test_data` (real PSG 731/734/736/737 + 73/100): check the hypnogram TXT has one
   label per expert line (±1), κ in a plausible range on clean files (literature for GSSC ≈ 0.7–0.8),
   lower on the defective ones (dead channel, clipping), reports render, skip works on a second run,
   global summary rebuilt from disk.
3. Tests: proposed level **quick** (`-m quick`: the new notebook executes; nothing in the chain changes).
