# Plan: rejection decisions cover every epoch (stage choice moves out of tools 8 / 8bis decisions)

Status: **final plan, validated by the user 2026-10-05 (incl. the design choices below), in progress.**

## Why
A rejection decision should answer "is this epoch clean?", not "will I analyse it?". Today tools 8 and 8bis
write `in_scope = False` for the stages not selected, so a review done for an N2/N3 study must be redone for
a later REM study. The analysis choice of stages already lives downstream (tool 9 stage checkboxes,
statistics filtering on `stage`). Same spirit as "compute everything, flag, choose at output time".

## Decisions (user, 2026-10-05)
1. **Tool 8 stage checkboxes = review focus.** They filter the navigator and the Section-2 report only. The
   saved decision covers every epoch. A flagged epoch outside the focus (never reviewed) keeps tool 7's
   automatic decision (rejected).
2. **Tool 8 `out of scope (view only)` mode is removed.** `Show` = `flagged` / `kept` / `all`, all three
   restricted to the focus stages. The grey "out of scope" banner and the grey strip cells go away.
3. **Tool 8bis reference stages.** The stage selection (default N2/N3/R, unchanged) becomes the
   **reference stages**: channel badness and participant exclusion are judged on them. The epoch rule is
   then applied to **every** epoch. Exclusion: `all channels dropped` (unchanged, badness is already
   reference-based) and `all reference-stage epochs rejected` (replaces `all in-scope epochs rejected`).
4. **Counts: both denominators + per stage**, in both tools. `pct_rejected` (8) / `pct_epochs_rejected`
   (8bis) now count **every** epoch. Additive `n_epochs_ref` / `n_rejected_ref` / `pct_rejected_ref` on the
   reference (8bis) / focus (8) stages. `n_<stage>` / `n_rejected_<stage>` for **every** stage present (tool 8
   had them for the selected stages only; 8bis had none).
5. **`in_scope` column kept, always True** in new decision tables (my default, as announced): tool 9 needs no
   change and old decision tables stay readable (its "out of scope" counts and message remain for them).

## Design details (my choices, flag any you disagree with)
- **Tool 8, one name per meaning.** In the notebook `S['in_scope']` becomes `S['focus']` (stage ∈ ticked
  stages) so the code never confuses "reviewed here" with "analysed downstream". The focus is also the
  **reference for the channel-triage badness** (`channel_badness(..., focus)`), exactly like 8bis's
  reference stages (W movement artefacts would otherwise suggest dropping healthy channels). Section-1 label:
  `Stages to review (also the reference for channel badness): the decision covers every stage`.
- **Changing only the focus keeps the manual overrides.** Today "Apply selection" resets every override.
  The decision no longer depends on the stages, so Apply resets the overrides **only** when methods, event
  types, dropped channels or the epoch rule changed. Otherwise widening the focus to W after reviewing N2/N3
  would silently throw the N2/N3 review away. The status line says which case happened.
- **`recompute_reject` (shared lib):** `base_reject = any_hit` (no `& in_scope` any more); the third return
  value keeps its position but is documented as the focus mask. With every method / event type ticked it
  still reproduces tool 7's `reject_flag` exactly, now **whatever the stage ticks**. The batch twin does not
  call it (checked), so it is unaffected.
- **Tool 8 saved tables:** `in_scope = True`, `rejected = final_reject` on every epoch, `reject_reason`
  filled for every rejected epoch. `n_in_scope` / `n_out_of_scope` stay in the record (now `n_epochs` / 0) for
  schema continuity. `stages_used` keeps its name and now means "focus stages". `n_seen` / `pct_seen` keep
  their denominator = focus ∧ flagged (the set the navigator walks).
- **Tool 8 `_print_decision`:** `==> REJECTED x / N epochs (all stages) = p %` then
  `in focus stages (N2, N3): y / M = q %`.
- **Tool 8 Section 2 report** (PSD grid, distributions, stage × method table): focus stages only, as today
  (`fit_mask = focus`). **Section-4 review strip**: every epoch drawn in its decision colour, no grey
  (`plot_review_strip(..., in_scope=None)`; the lib parameter stays for compatibility).
- **8bis `auto_reject_decision`:** badness over the reference epochs (unchanged), epoch rule over every epoch
  with the remaining good channels. No reference stage present in a file → still `failed` (badness is
  undefined), unchanged. `stages_of_interest` keeps its column name (old global summaries keep concatenating)
  and is labelled **Reference stages** in the UI and the report.
- **8bis outputs:** `build_epoch_decision` → `in_scope = True`, `rejected` on every epoch, reason
  `k/n good channels flagged`. Record / global summary: `n_epochs`, `n_epochs_rejected`, `n_epochs_kept`,
  `pct_epochs_rejected` over every epoch + the additive `_ref` and per-stage columns. Compact table shows both
  percentages. The **bar plot uses the reference-stage %** (what the quality is judged on, comparable across
  participants whatever their amount of W). The heatmap's rejected-epoch strip and the per-stage table cover
  every stage.
- **Tool 9:** no code change expected, verified by the tests. Its `stages requested but not analysed`
  message keeps its all-rejected / absent cases. Its out-of-scope case now only arises on old decision tables.

## Tool 8 Section-3 display (added 2026-10-05, display only, no output changes)
Goal: montage + PSD + metric table on **one screen**, the epoch navigator right above the plots.
- **Smaller detail panel**: `plot_epoch_detail` figure height 4.8 → ~3.2 in (sparse 13 in wide, HD 18 in),
  as one constant `DETAIL_FIG_HEIGHT_IN` in the lib; fonts/legend checked on a rendered PNG (synthetic + a
  dense montage) so the table stays readable. Only the navigator calls it (checked): no saved figure changes.
- **New widget order** (`display(widgets.VBox([...]))` of the Section-3 cell):
  1. `Start navigator` + Show / Stage / Method / Contains event, help line
  2. `◀ Prev` / Index slider / `Next ▶` / epoch n° offset, then position line + `reviewed n/N` strip
  3. **montage** (`out_nav`) + **PSD/table** (`out_detail`)
  4. **Decision** row (Keep / Reject / Apply to all shown / Cancel) + its confirm line + the status banner
  5. Annotate-this-epoch box
  6. Display settings: `± context` + `Show EOG/EMG context`, `Channels` + `Table rows`, `Scale (µV/row)`
  (`± context` / `Show EOG/EMG context` move down too: they are display settings, and leaving them above
  would separate the navigator from the plots.)

## Golden snapshots: proposal NOT to regenerate (deviation from the draft)
The synthetic golden predates the decision tables: it still holds the former `clean_epo_*` copies, which
`test_exclusion_tools.py` / `test_tool9.py` use as **the reference** for the kept epochs. Regenerating now
would delete that reference (tools 8/8bis no longer write clean-epo). And what changes is well delimited:
- **Tool 8**: the chain ticks every stage (default), so its outputs are already "every epoch in scope":
  **unchanged** except additive columns (allowed by the comparison). No intended change.
- **Tool 8bis**: the chain runs N2/N3/R. The decision on those stages is unchanged (same badness, same rule),
  W/N1 epochs gain a decision. → add `reports_rejection_auto/` to `INTENDED_CHANGES['8bis']`
  (`n_epochs` 420 → total, etc.), with the dedicated tests below.
- **Tool 9**: already in `INTENDED_CHANGES`.
So the golden stays as is. (If you prefer a fresh golden anyway, it must come with moving the clean-epo
reference elsewhere: more work, no gain I can see.)

## Steps
0. Baseline: run `-m "tool8 or tool8bis or tool9"` (standard) on the current code, to start from green.
1. `qc_rejected_epochs_lib.py`: `recompute_reject` (base over every epoch), `build_manual_decision_row`
   (counts over every epoch + `_ref` + per-stage for every stage), docstrings.
2. `8_reject_manually_voila.ipynb` (rename-to-`.txt` method): `S['focus']`, Section-1 label + help line,
   `_recompute_reject` override-preserving logic, `_print_decision`, `dd_show` without `out`, `_status_html`
   / `_progress_html` / Method-dropdown `out` branches removed, `_build_decision_tables`, review strip, all
   markdown. Re-run `tools_curry/_make_tool8_curry.py` (asserts 3/3).
3. `8bis_reject_automatically_voila.ipynb`: `auto_reject_decision`, `build_epoch_decision`,
   `decision_summary_row`, exclusion reason, `per_stage_table_html`, `decision_overview_html`, heatmap strip,
   compact table + bar plot, labels/markdown. Re-run `tools_curry/_make_tool8bis_curry.py`.
4. Tests:
   - `test_pipeline.py`: `INTENDED_CHANGES['8bis'] = ['reports_rejection_auto/']`.
   - `test_exclusion_tools.py::test_tool8bis_decision_tables_match_the_clean_epochs` → **dedicated test of
     the intended change**: on the reference stages the kept epochs equal the golden clean-epo; every epoch
     is `in_scope`; a W/N1 epoch is rejected iff more than 20 % of the good channels flag it (recomputed from
     the reports TSV); `_ref` counts equal the old golden `n_epochs` / `n_epochs_rejected`.
   - `test_tool8_decision_tables`: unchanged (all stages ticked); add a check that `in_scope` is all True.
   - `test_variants.py::test_tool8_out_of_scope_is_view_only` → `test_tool8_unticked_stage_stays_in_the_decision`:
     W unticked → `Show` has no `out` option, the navigator lists no W epoch, the saved table has every W epoch
     `in_scope` True with tool 7's flagged W epochs rejected; and a Keep override survives an Apply that only
     changes the focus (but not one that changes a method).
   - `test_tool9_names_a_stage_out_of_scope_in_tool8` → kept as the **legacy-table** test: it writes an
     old-style decision table (N3 `in_scope` False) on a copy, so tool 9's backward-compatible path stays
     covered.
   - Proposed level: **full** on `tool8`, `tool8bis`, `tool9` (+ `test_curry_twins.py`), real data included.
5. Docs: SPEC §8, §8bis, §9, *Cross-cutting → Rejected-data policy* (table row: `n_epochs_out_of_scope` = 0
   with new decisions) and *Feature tools read tool 7's epochs + a decision*; CLAUDE.md tool-8 bullet
   ("out-of-scope stages = `in_scope False`" → focus only) and 8bis bullet ("computed over the stages of
   interest" → reference stages, decision on every epoch). SPEC edits proposed for confirmation as usual.
6. Real-data check on one participant (tools 8bis → 8 → 9), then suggest a commit message; delete this plan
   and update the memory once you confirm.

## Backward compatibility
- Old decision tables (with `in_scope` False rows) stay valid input for tool 9: nothing removed.
- Old `_autoreject_decision.tsv` / `_manualreject_decision.tsv` keep concatenating in the global summaries
  (same column names, new columns NaN for them). Note: on old rows `pct_epochs_rejected` was over the
  reference stages; the global summary does not distinguish them (re-run 8bis to homogenise, it is cheap).
