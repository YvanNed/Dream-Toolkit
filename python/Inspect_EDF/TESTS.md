# Tests: what is checked, and why the outputs can be trusted

The toolkit's notebooks are long and are edited often. A change made for one purpose can silently move a
number somewhere else: a sleep efficiency, a flagged-epoch rate, a band power. The test suite in `tests/`
catches that before it reaches a dataset. This page explains, without code, **what the tests check and how**.
The technical details (commands, file layout, how to add a test) are in SPEC.md, *How to run → Tests*.

## In one paragraph

The tests run the **real notebooks** (the same code you run in Voila) on two small test databases. They
click through the tools as a user would: **5 → 5bis → 6 → 7 → 8bis → 8 → 9**, the full chain from sleep
metrics to spectral features. Then they check three things: the results are **the same numbers as a reference
run** taken before the code changed, each new behaviour does **exactly what it claims**, and on artificial
nights whose content is known, the tools **find what was put in** (the injected defects, the generated sleep
rhythms). No part of the tools is copied or simulated: if a test passes, the notebook itself produced the result.

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

## Three levels: run what the change needs

Not every change needs every test. The suite has three levels, from the fastest to the most thorough:

| Level | What it runs | Duration | When |
|---|---|---|---|
| **Quick** | small checks of the shared code, and every Voila notebook (EDF tools and Curry twins) opened and run from top to bottom without clicking anything | about 2 min | a change that touches no tool's results (a test, a text, a Curry generator) |
| **Standard** | quick + the full chain on the synthetic nights (comparison with the reference run, every check below) + the known answers + tools 1, 3, 4 | about 9 min (5 when the chain is reused) | a change to a tool |
| **Full** | standard + the chain on the real nights + the parameter variants | about 24 min | before merging a branch, or after a change to code shared by several tools |

Two things keep the runs short:

- **Only the tools concerned.** Each test is labelled with the tool it checks. A change to tool 7 runs the tests
  of tool 7 and of the tools that use its outputs (8bis, 8, 9), not those of tools 5 or 6. A change to code
  shared by several tools runs the tests of every tool that uses it.
- **The chain is not rerun when nothing changed.** The test run remembers which version of each tool produced
  the current results. If only tool 9 changed, tools 5 to 8 keep their results and only tool 9 is rerun (about
  2 minutes instead of 4.5). A change to the test machinery itself reruns everything.

At each change, the level to run is proposed with a recommendation, and you choose.

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

### 7. The tools find what was put in the synthetic nights (known answers)

Comparing with a reference run proves that nothing **changed**, not that the result was right in the first
place: a mistake present when the reference was taken would pass forever. On the synthetic nights the truth is
known, because the test machinery generated it, so these tests recompute the expected answer from the
generator itself:

- **Tool 5**: total sleep time, sleep efficiency, sleep onset latency, WASO, time and percentage in each stage,
  stage latencies: all recomputed from the generated hypnogram, and **exactly** equal.
- **Tool 6**: exactly two channels flagged, the clipped Fp1 and the dead C3. Nothing else.
- **Tool 7**: the 10 movement bursts are re-created with the same random seed on a silent signal, which gives
  their exact position. Every epoch holding at least 1 s of burst is flagged on Fp1, and no other Fp1 epoch.
- **Tool 9**, the generated physiology:
  - occipital **alpha** stronger in wake than in N2 (about 15 dB), and stronger on O1 than on Fp1 in wake
    (about 11 dB);
  - a **spindle peak** in N2 on C3: 12–14 Hz above both 9–11 Hz (about 3.7 dB) and 16–20 Hz (about 10 dB).
    The dead C3 of one participant is skipped, as it should be;
  - a steeper **aperiodic slope** (higher exponent) in N3 than in wake.

  The test thresholds are set well below these measured differences (for example 6 dB for a 15 dB effect),
  so a genuine change of method fails while ordinary variation passes.

### 8. The preparation tools 1, 3 and 4

The chain starts from data already prepared by tools 1 to 4. Three short runs check them on a fresh
synthetic database:

- **Tool 1** (header inspection): 4 files of 4 channels at 256 Hz, the Fp1 recorded with a ±50 µV range
  reported as a too-narrow dynamic range, every header anonymous.
- **Tool 3** (hypnogram remapping): a raw hypnogram in the Compumedics style (stages coded `0` to `5`, N3
  sometimes coded `4` as in the old R&K stage 4, unscored `?` epochs) is remapped. The result must give back
  the generated AASM stages epoch by epoch, and the `?` in the middle of the night must be listed for review.
- **Tool 4** (event labels): the raw labels are listed with their files; the harmonised file gets the
  canonical names, a French label (`Hypopnée`) is recognised, a label ticked "ignore" is saved as ignored, and
  a mapping saved earlier is kept.

Tools 0 and 2 are too interactive for a useful run: they are only checked to open and run (quick level).
Tool 1bis (anonymisation) checks its own work on every file (the signal part of the file must stay
byte-identical) and is not part of the routine tests.

### 9. Options other than the defaults (parameter variants)

The chain runs every tool with its default settings. The full level also runs the main options, each on a
copy of the database:

| Tool | Option | What must happen |
|---|---|---|
| 7 | epochs of **10 s** | 3 times more epochs, each inheriting the stage of its 30 s epoch; tools 8bis and 9 then work on them |
| 7 | **resampling** to 128 Hz | the saved epochs and the settings file are at 128 Hz |
| 7 | **notch off** | the 50 Hz line noise of one participant stays more than 10 dB above the notched run |
| 7 | a **custom stage** `N4` | N4 gets its own amplitude threshold, the saved epochs read back whole, N4 reaches tools 8bis and 9 |
| 6 | **resampling + high-pass** | the dead channel stays flagged, no healthy channel is flagged |
| 9 | **multitaper** spectrum | band powers within 3 dB of the default (Welch) |
| 9 | **knee** aperiodic mode | the knee is estimated and the fits are good |
| 8 | **rescue** an epoch in the navigator | it is saved as kept, marked as a manual rescue |
| 8 → 7 | a **manual annotation** | added in tool 8, re-read by tool 7 when its option is ticked, and it flags its epoch |

The custom-stage variant found a real bug when it was written: an epoch with a stage outside W/N1/N2/N3/R made
the saved epochs file of tool 7 unreadable for tools 8 and 9. It is fixed, and this test keeps it fixed.

The tool 6 variant also documents a deliberate choice: tool 6 measures the signal **after** its optional
filtering (a high-pass can rescue a channel worth keeping). The filtering smooths the flat tops of a clipped
channel, so with these options the clipped Fp1 of the synthetic set is no longer flagged.

## Smaller checks

| Area | What is checked |
|---|---|
| Every Voila notebook | Opens and runs from top to bottom without error, EDF tools and Curry twins (quick level) |
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
- **Curry recordings.** The Curry twins are checked for consistency with their EDF original and opened in a
  kernel (quick level), but never run on a recording (the Curry test data is 17 GB).
- **Tools 0 and 2** are only opened (too interactive for a useful run), **tool 1bis** checks itself.
- **Every option.** The variants cover the main options of tools 6 to 9, not every combination.
- **The science of the methods.** The tests prove the tools compute what they computed before, or what a
  change intends. Whether a threshold or a method is the right choice for a study remains a scientific
  decision.

## When a test fails

Read the message: it names the tool, the file, the column and the first value that differs. Either the
change broke something (fix the code), or the change is meant to move that value (then it must be declared
as an intended change, with its own test explaining what the new value should be). The reference run is
regenerated **only** on code known to be good, before starting a change, never to make a failing test pass.
