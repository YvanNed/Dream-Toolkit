# TODO — `8_live_explore_1file` Section 3: migrate to a fully interactive canvas (v2)

**Status**: deferred. v1 ships with static matplotlib PNG rendering (chosen for speed + reliability under Voila + easy HTML-report embedding).

## Why v1 is static PNG
- Each navigation step (prev / next / slider) re-renders the epoch figure to a PNG shown in the browser (~100–300 ms — snappy).
- No extra dependency; predictable layout; figures embed directly into the optional saved HTML snapshot.
- Limitations accepted for v1: no mouse zoom/pan, no click-to-seek on the navigator strip, no drawing annotations with the mouse — every interaction is a full redraw.

## v2 goal
Migrate **Section 3 (epoch explorer)** only to an interactive canvas to unlock:
1. **Click-to-seek** on the navigator strip (click the hypnogram/spectrogram to jump to that epoch).
2. **Mouse zoom/pan** on the epoch time series without a full redraw.
3. **Direct annotation writing** — draw/label events (e.g. artifact spans, custom markers) on the trace with the mouse and persist them (candidate format: append to / mirror the `*_event_xml.csv` schema `Name, Start, Duration`).

## Candidate backends (evaluate at v2 time)
- **`ipympl` / `%matplotlib widget`**: least code change from the static matplotlib figures; supports click events and pan/zoom. Cons: extra dependency, occasionally flaky under Voila, click events add round-trip latency, harder to snapshot to a static report.
- **Plotly**: fully interactive + hover tooltips; supports shape drawing for annotations. Cons: large dependency, slower first render for dense multichannel traces.

## Constraints to respect in v2
- Must stay usable under Voila (code hidden) **and** in the Jupyter twin.
- Keep navigation fast — interactive mode must not regress the snappy feel of the static version; consider keeping a static-PNG fallback toggle.
- Annotation persistence must never overwrite the original Compumedics `*_event_xml.csv`; write a separate file (e.g. `<stem>_manual_events.csv`).
- Add the chosen backend (`ipympl` or `plotly`) to `environment.yml` when implemented.

---

# TODO — `8_live_explore_1file`: automatic sleep staging (YASA) add-on (deferred)

**Status**: deferred idea (from the old `plan_inspect_file_live` "Out of scope (v1)"). Not committed — recorded so it isn't lost.

Section 3 already supports **manual** rescoring (stage buttons / keyboard shortcuts → `<hypno_stem>_rescored.txt` + a rescore log, never overwriting the original). A natural follow-up is to **pre-fill / suggest** the staging automatically with **YASA's auto-scorer** (`yasa.SleepStaging`), which the user then reviews and corrects — turning Section 3 into a scoring-*review* tool when no expert hypnogram exists.

Notes for whoever picks this up:
- **YASA is already a dependency** (used for the hypnospectrogram), so no new package — but the auto-scorer needs a central EEG (± EOG/EMG) channel and expects standard referencing; surface the channel choice in the UI.
- Treat the auto-scored hypnogram as a **suggestion**: load it as the starting point for the existing manual-override flow, keep the same "never overwrite the original" rule, and write to a distinct file (e.g. `<stem>_autostaged.txt`) with a provenance note that it is YASA-generated.
- Show YASA's per-epoch **confidence** (it exposes probabilities) to guide where the reviewer should look first.
