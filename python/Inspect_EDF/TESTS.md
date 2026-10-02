# Tests: what is checked, and why the outputs can be trusted

The toolkit's notebooks are long and are edited often. A change made for one purpose can silently move a
number somewhere else: a sleep efficiency, a flagged-epoch rate, a band power. The test suite in `tests/`
catches that before it reaches a dataset. This page explains, without code, **what the tests check and how**.
The technical details (commands, file layout, how to add a test) are in SPEC.md, *How to run → Tests*.

## In one paragraph

The tests run the **real notebooks** (the same code you run in Voila) on two small test databases. They
click through the tools as a user would: **5 → 5bis → 6 → 7 → 8bis → 8 → 9**, the full chain from sleep
metrics to spectral features. Then they check two things: the results are **the same numbers as a reference
run** taken before the code changed, and each new behaviour does **exactly what it claims**. No part of the
tools is copied or simulated: if a test passes, the notebook itself produced the result.

## How it works

- **The notebooks are driven, not copied.** Each tool runs in its own Jupyter kernel. A small driver picks
  the folders in the file choosers, sets the widgets (thresholds, criteria, subsets) and clicks the buttons.
  The results are read back from the files the tool wrote on disk.
- **Two test databases**, with the same layout and the same four defects:

  | | Synthetic | Real |
  |---|---|---|
  | What | 4 nights of 4 h, generated from fixed random seeds | 4 full nights (7.8 h) from `tools/test_data` |
  | Signal | artificial EEG: 1/f background, alpha in wake, theta in N1/REM, spindles in N2, slow activity in N3 | a real person's PSG, with defects injected |
  | Defects | clipped Fp1, dead C3, 50 Hz on O1, movement bursts on Fp1 | the same |
  | Reference | M2 (mastoid), as in a clinical montage | the same |
  | Shared? | yes: generated on any machine, its reference results are versioned | **no**: real participant data, never leaves this machine |

  The synthetic nights were tuned until every tool sees them as plausible sleep EEG: tool 6 flags the same
  channels, tool 7 flags a similar share of epochs (0 to 14 % per stage), and the aperiodic exponent sits
  around 2.2. A test on a signal the tools find absurd would prove little.
- **A reference run** ("golden snapshot") is taken once, on code known to be good, before a change. Each
  later run is compared to it, table by table and value by value.

Running everything takes about 18 minutes; the synthetic database alone, about 6 minutes.

## The main checks

### 1. Nothing changed by accident

For every tool, every table it writes (TSV) and every settings file (JSON) is compared with the reference
run. A changed value, a removed column or a missing file fails the test, and the message names the first
value that differs (for example `column 'tst_min': 1 value differs, first at row 2: 265.5 -> 266.0`).
Adding a column or a file is allowed, so a new feature does not break the comparison as long as it leaves
the existing numbers alone. Two runs of the chain on the same data were checked to give byte-identical
results, so any difference is real, not noise.

When a change is *meant* to move numbers (for example the per-stage summary of tool 7 now leaves the
excluded participants out), that output is listed as an intended change and gets its own dedicated test
instead.

### 2. Tool 9 gives the same features as before the exclusion refactor

Tool 9 used to read a cleaned copy of the epochs (only the kept ones). It now reads every epoch and joins
the keep/reject decision of tool 8 / 8bis. The test takes the new results, keeps only the kept epochs of the
kept channels, and compares them with the reference run made with the old method: every per-epoch table,
every per-stage table and every database table must give **the same values**. This is the proof that the
new way of working changed the organisation of the data, not the science.

(One deliberate difference is excluded from the comparison: the night thirds now span the whole night
instead of the kept epochs only.)

### 3. A deselected channel is referenced like the others

Tool 7 now keeps a channel you deselect (for example a dead C3 found by tool 6), marked "bad", so that its
features can still be computed later. MNE, the signal library, leaves bad channels in their original
reference, so tool 7 re-references them itself. Two tests check that this is exactly right:

- **M2 reference**: the clipping participant is processed twice, once with Fp1 deselected and once with Fp1
  kept. Every channel, Fp1 included, must be identical in the two runs: the deselected channel received the
  same reference and the same filtering, and the good channels did not depend on it.
- **Average reference** (used for high-density montages): with Fp1 deselected, the good channels must
  average to zero without Fp1 (the dead channel never entered the reference). And the bipolar difference
  Fp1 − C3 must be the same whether Fp1 was deselected or not: any common reference cancels in a
  difference, so this proves Fp1 carries the same reference as the others.

The flagging tables of tool 7 are also unchanged by keeping these channels (check 1).

### 4. "Re-aggregate only" gives the same result as a full computation

Tool 9 can rebuild its stage averages and reports after a decision change without recomputing any spectrum
(a few seconds instead of minutes). The test rejects 25 more epochs in a decision table, runs
*Re-aggregate only*, then runs the full computation with the same decision: every table must match. This
test found two real bugs during development.

### 5. Exclusion decisions behave as documented

A participant can be excluded from several places: tool 5bis (criteria, or by hand), tool 7, tool 8, and
tool 8bis automatically. The tests check the rules of the shared exclusion registry:

- each tool only replaces its **own** decisions, never another tool's;
- a **forced inclusion** (5bis, with a mandatory comment) overrides every exclusion;
- an excluded participant is **still processed** everywhere: present in every table, flagged
  `participant_excluded`, left out of the group statistics and figures, with a red banner on its own report;
- tool 5 knows nothing of exclusions: it extracts the raw sleep metrics, 5bis takes the decisions;
- a decision taken after a tool ran is picked up when the tool is run again, and the tool warns at scan time
  that its tables are out of date;
- an absent, empty or damaged registry never blocks a run (a damaged one is reported and never overwritten).

### 6. Running a subset only touches that subset

Re-running one participant must leave the others alone. The test runs a subset with *Skip already
processed* ticked (nothing must run, and the orange warning must name the skipped participants), then
unticked (only the chosen participant's files change, checked through their modification dates), and checks
that the database tables still cover everyone.

## Smaller checks

| Area | What is checked |
|---|---|
| Every tool | No Python error appears in any tool's output during the chain |
| Every tool | Every reference file is owned by a tool (none escapes the comparison) |
| Tool 5bis | Validating writes the criteria exclusions to the registry; a new validation replaces them but keeps the manual decisions; a forced inclusion is refused without a comment; the exclusion columns are never mistaken for sleep metrics |
| Tool 6 | The quality tables flag the excluded participants; `dataset_overview.html` counts only the others and names the excluded one; the excluded participant's report shows the banner |
| Tool 7 | Deselected channels are in the `.fif` (marked bad) and in the settings file, absent from the flagging tables; the per-stage summary is pooled over the non-excluded participants only; a manual exclusion needs a reason; a re-run with nobody to process still refreshes the tables |
| Tool 8bis | The decision tables keep exactly the epochs and channels the former cleaned copy held; tool-7 deselected channels are recorded as dropped; an event-count threshold excludes a participant, a re-run without it lifts the exclusion |
| Tool 8 | The decision tables match the former cleaned copy; the *Exclude* / *Remove my exclusion* button writes and deletes the decision |
| Tool 9 | Every epoch × channel is in the tables with the right decision flags; dropped channels are still computed; the per-epoch spectrum file exists; night-third counts (rejected, out of scope, % rejected) are those of each stage × third cell |
| Participant selector | Skip applies in both modes; a subset with already-processed participants shows the warning; the add box only offers participants not yet selected; pasted and clicked ids switch to subset mode; unknown ids are refused; a rescan drops vanished ids |
| Registry details | Ids match regardless of case; no temporary file is left behind; rewriting the same decisions does not touch the file (so it never looks "changed") |
| Curry versions | The Curry twins of tools 6, 7, 8 and 8bis are valid Python and contain the new features of their EDF original (they are generated, not run, see below) |

## What the tests do not check

- **The look of the notebooks and reports.** Widgets are driven and HTML reports are written, but nobody
  looks at them: layout, colours and figure readability still need a human eye in Voila.
- **Figures.** Only the tables behind them are compared.
- **Curry recordings.** The Curry twins are checked for consistency with their EDF original, not run (the
  Curry test data is 17 GB).
- **Tools 0 to 4 and 1bis** are not part of the chain.
- **The science of the methods.** The tests prove the tools compute what they computed before, or what a
  change intends. Whether a threshold or a method is the right choice for a study remains a scientific
  decision.

## When a test fails

Read the message: it names the tool, the file, the column and the first value that differs. Either the
change broke something (fix the code), or the change is meant to move that value (then it must be declared
as an intended change, with its own test explaining what the new value should be). The reference run is
regenerated **only** on code known to be good, before starting a change, never to make a failing test pass.
