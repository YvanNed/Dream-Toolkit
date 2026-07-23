"""
_make_tool4_curry.py — Generate 4_remap_events_curry_voila.ipynb from the EDF original.

Adaptations from the EDF version:
- File discovery: .edf -> .cdt (bare suffix)
- Event source: Compumedics CSV / *.edf.XML companions -> Curry French text export
  (*_ScoredEvents_Export.txt), parsed via curry_io.load_events_curry using the .cdt.dpo
  recording-start datetime. load_events() keeps its (path, suffix) -> (events_list, source)
  signature (events_list = list of (name, start, duration)) so the scan/harmonize/verify
  logic is untouched.
- Section 1bis (CSV vs XML consistency) is REMOVED: Curry has a single event source, so the
  cross-source check is meaningless. That drops the XML/CSV parser helpers and _events_multiset.
- Curry shared-module imports (curry_header, curry_io); the now-unused xml.etree import is dropped.

Everything else (config grouping, harmonize UI, merge/save of config_param/event_remap.json,
verification) is format-agnostic and kept as-is.

Run from Inspect_EDF root:
    & "$env:LOCALAPPDATA\\miniforge3\\envs\\inspect_edf\\python.exe" tools_curry/_make_tool4_curry.py
"""

import json, os, sys, ast

SRC = "tools/4_remap_events_edf_voila.ipynb"
DST = "tools_curry/4_remap_events_curry_voila.ipynb"

with open(SRC, encoding="utf-8") as f:
    nb = json.load(f)

errors = []


def get_src(cell):
    s = cell["source"]
    return "".join(s) if isinstance(s, list) else s


def set_src(cell, text):
    orig = cell["source"]
    if isinstance(orig, list):
        cell["source"] = text.splitlines(keepends=True)
    else:
        cell["source"] = text


def replace_in_cell(cell, old, new, label):
    src = get_src(cell)
    if old not in src:
        errors.append(label)
        print(f"  X NOT FOUND: {label}")
        return False
    if src.count(old) > 1:
        print(f"  ! MULTIPLE matches: {label} — replacing all")
    set_src(cell, src.replace(old, new))
    print(f"  ok {label}")
    return True


def remove_between(cell, start_marker, end_marker, label):
    """Delete text from start_marker up to (but not including) end_marker."""
    src = get_src(cell)
    i = src.find(start_marker)
    j = src.find(end_marker)
    if i == -1 or j == -1 or j <= i:
        errors.append(label)
        print(f"  X NOT FOUND / bad range: {label}")
        return False
    set_src(cell, src[:i] + src[j:])
    print(f"  ok {label}")
    return True


md0  = nb["cells"][0]
code = nb["cells"][1]

# ===========================================================================
print("=== Markdown cell 0 ===")
NEW_MD = (
    "# Harmonize scored-event labels of a Curry 9 (.cdt) EEG/PSG database\n"
    "\n"
    "This notebook lets you **visualize** the scored-event configurations present in a Curry 9\n"
    "database (events annotated during sleep scoring and exported by Curry as a text file:\n"
    "arousals, apnea, hypopnea, limb movements, desaturations…) and **harmonize** their raw\n"
    "labels to a single canonical vocabulary.\n"
    "\n"
    "It returns a JSON file `config_param/event_remap.json` (a flat python dict\n"
    "`{raw_label: canonical_label}`, with `null` for labels you choose to ignore) that downstream\n"
    "tools (epoch rejection in `6_preprocessing_curry`) can use.\n"
    "\n"
    "---\n"
    "**To use this notebook, interact with the widgets and read the output below. You first have to\n"
    "select a database to make the widgets appear.**\n"
    "\n"
    "Sections:\n"
    "1. Select your study folder and scan the events\n"
    "2. Event configurations found\n"
    "3. Harmonize the labels\n"
    "4. Preview & save the JSON\n"
    "5. Verify\n"
    "\n"
    "The events are read from the Curry text export `*_ScoredEvents_Export.txt` next to each `.cdt`.\n"
    "The raw labels are the export strings (often in French, e.g. `Microéveil 1 ARO SPONT`); the\n"
    "canonical suggestions below are English defaults, so most raw labels will need a manual choice.\n"
)
set_src(md0, NEW_MD)
print("  ok markdown rewritten")

# ===========================================================================
print("\n=== Cell 1: imports ===")
# Drop the now-unused XML parser import
replace_in_cell(code, "    import xml.etree.ElementTree as ET\n", "", "remove ET import")
# Add curry shared-module imports inside the try block
replace_in_cell(code,
    "    from IPython.display import display, HTML, clear_output\n"
    "except ImportError as e:",
    "    from IPython.display import display, HTML, clear_output\n"
    "    import sys as _sys\n"
    "    # curry shared modules — located next to this notebook\n"
    "    _here = os.path.dirname(os.path.abspath('__file__'))\n"
    "    if _here not in _sys.path:\n"
    "        _sys.path.insert(0, _here)\n"
    "    from curry_header import read_curry_header\n"
    "    from curry_io import load_events_curry, rec_start_from_header\n"
    "except ImportError as e:",
    "curry imports")

# ===========================================================================
print("\n=== Cell 1: event loader ===")
OLD_EVENTS = (
    '# ---- event companion loading: CSV first, XML (<ScoredEvents>) fallback ----\n'
    'def event_companion_paths(edf_path, csv_suffix="_event_xml.csv"):\n'
    '    """Return (csv_path_or_None, xml_path_or_None) for an EDF stem.\n'
    '    csv_suffix is the configurable Compumedics event-CSV suffix (default \'_event_xml.csv\')."""\n'
    '    edf_path = Path(edf_path)\n'
    '    csv = edf_path.with_name(f"{edf_path.stem}{csv_suffix}")\n'
    '    csv = csv if csv.exists() else None\n'
    '    xml = None\n'
    '    for cand in (f"{edf_path.name}.XML", f"{edf_path.name}.xml"):\n'
    '        p = edf_path.with_name(cand)\n'
    '        if p.exists():\n'
    '            xml = p\n'
    '            break\n'
    '    return csv, xml\n'
    '\n\n'
    'def load_events_from_csv(csv_path):\n'
    '    """Parse a Compumedics *_event_xml.csv -> list of (name, start, duration)."""\n'
    '    df = pd.read_csv(csv_path)\n'
    '    for col in ("Name", "Start", "Duration"):\n'
    '        if col not in df.columns:\n'
    '            raise ValueError(f"missing column \'{col}\' in {Path(csv_path).name}")\n'
    '    events = []\n'
    '    for _, row in df.iterrows():\n'
    '        events.append((str(row["Name"]).strip(), float(row["Start"]), float(row["Duration"])))\n'
    '    return events\n'
    '\n\n'
    'def load_events_from_xml(xml_path):\n'
    '    """Parse the <ScoredEvents> of a Profusion CMPStudyConfig .edf.XML\n'
    '    -> list of (name, start, duration). <Input> is ignored (absent from the CSV)."""\n'
    '    root = ET.parse(xml_path).getroot()\n'
    '    events = []\n'
    '    for se in root.iter("ScoredEvent"):\n'
    '        name_el = se.find("Name")\n'
    '        if name_el is None or name_el.text is None:\n'
    '            continue\n'
    '        start_el = se.find("Start")\n'
    '        dur_el = se.find("Duration")\n'
    '        start = float(start_el.text) if (start_el is not None and start_el.text) else float("nan")\n'
    '        dur = float(dur_el.text) if (dur_el is not None and dur_el.text) else float("nan")\n'
    '        events.append((name_el.text.strip(), start, dur))\n'
    '    return events\n'
    '\n\n'
    'def load_events(edf_path, csv_suffix="_event_xml.csv"):\n'
    '    """CSV-first, XML-fallback event loader.\n'
    '    Returns (events_list, source) with source in {\'csv\', \'xml\'} or (None, None) when\n'
    '    neither companion is usable. events_list = list of (name, start, duration)."""\n'
    '    csv, xml = event_companion_paths(edf_path, csv_suffix)\n'
    '    if csv is not None:\n'
    '        try:\n'
    '            return load_events_from_csv(csv), "csv"\n'
    '        except Exception:\n'
    '            pass  # fall through to the XML\n'
    '    if xml is not None:\n'
    '        try:\n'
    '            return load_events_from_xml(xml), "xml"\n'
    '        except Exception:\n'
    '            pass\n'
    '    return None, None\n'
    '\n\n'
)
NEW_EVENTS = (
    '# ---- event loading: Curry French text export (*_ScoredEvents_Export.txt) ----\n'
    'def load_events(cdt_path, event_suffix="_ScoredEvents_Export.txt"):\n'
    '    """Load Curry scored events next to the .cdt file (French text export).\n'
    '    Returns (events_list, \'txt\') with events_list = list of (name, start, duration) in\n'
    '    seconds, or (None, None) when the export is absent/unreadable. The .cdt.dpo header gives\n'
    '    the recording-start datetime used to convert the export clock times (load_events_curry)."""\n'
    '    cdt_path = Path(cdt_path)\n'
    '    ev_path = cdt_path.with_name(f"{cdt_path.stem}{event_suffix}")\n'
    '    if not ev_path.exists():\n'
    '        return None, None\n'
    '    try:\n'
    '        hdr = read_curry_header(str(cdt_path))\n'
    '        rec_start = rec_start_from_header(hdr)\n'
    '        if rec_start is None:\n'
    '            return None, None\n'
    '        df = load_events_curry(str(ev_path), rec_start)\n'
    '        if df is None or df.empty:\n'
    '            return None, None\n'
    '        events = [(str(n).strip(), float(s), float(d))\n'
    '                  for n, s, d in zip(df["Name"], df["Start"], df["Duration"])]\n'
    '        return events, "txt"\n'
    '    except Exception:\n'
    '        return None, None\n'
    '\n\n'
)
replace_in_cell(code, OLD_EVENTS, NEW_EVENTS, "event loader")

# Remove _events_multiset (only used by the CSV-vs-XML check)
replace_in_cell(code,
    'def _events_multiset(events):\n'
    '    """Multiset of (name_lower, round(start,3), round(dur,3)) for CSV/XML comparison."""\n'
    '    return Counter((n.strip().lower(), round(s, 3), round(d, 3)) for (n, s, d) in events)\n'
    '\n\n',
    "",
    "remove _events_multiset")

# ===========================================================================
print("\n=== Cell 1: Section 1 banner + widgets ===")
# section1 banner (CSV/XML wording -> Curry text export)
replace_in_cell(code,
    'section1 = widgets.HTML("""\n'
    '<hr style="height:4px; background-color:black; border:none;">\n'
    '<h2>1. Select your study folder and scan the events</h2>\n'
    '<p>Pick the folder of your .edf database. Each .edf is expected to have a Compumedics/Profusion\n'
    'event companion next to it: a CSV (read first; its suffix is set in the <b>CSV suffix</b>\n'
    'field below — default <code>_event_xml.csv</code>) or, as a fallback, the\n'
    '<code>&lt;ScoredEvents&gt;</code> of the <code>*.edf.XML</code>.\n'
    '<br>&#x2022; Selecting the folder auto-detects the CSV suffix and refreshes the info line below.\n'
    '<br>&#x2022; Click <b>Run scan</b> to read the events and list the configurations.\n'
    '<br>&#x2022; "Skip labels already mapped" hides labels already present in an existing\n'
    '<code>event_remap.json</code> (incremental harmonization when you add a new cohort).</p>\n'
    '""")',
    'section1 = widgets.HTML("""\n'
    '<hr style="height:4px; background-color:black; border:none;">\n'
    '<h2>1. Select your study folder and scan the events</h2>\n'
    '<p>Pick the folder of your Curry (.cdt) database. Each .cdt is expected to have a Curry event\n'
    'text export next to it, whose suffix is set in the <b>Event export suffix</b> field below\n'
    '(default <code>_ScoredEvents_Export.txt</code>).\n'
    '<br>&#x2022; Selecting the folder auto-detects the suffix and refreshes the info line below.\n'
    '<br>&#x2022; Click <b>Run scan</b> to read the events and list the configurations.\n'
    '<br>&#x2022; "Skip labels already mapped" hides labels already present in an existing\n'
    '<code>event_remap.json</code> (incremental harmonization when you add a new cohort).</p>\n'
    '""")',
    "section1 banner")

# Remove section1bis banner
replace_in_cell(code,
    'section1bis = widgets.HTML("""\n'
    '<hr style="height:4px; background-color:black; border:none;">\n'
    '<h2>1bis. (Optional) Check CSV vs XML consistency</h2>\n'
    '<p>Verifies that, for every file having <b>both</b> companions, the <code>*_event_xml.csv</code>\n'
    'and the <code>&lt;ScoredEvents&gt;</code> of the <code>*.edf.XML</code> describe the same events\n'
    '(name + start + duration, tolerance 1e-3 s). Writes <code>config_param/event_source_mismatch.tsv</code>.\n'
    'This is opt-in because it forces reading both files for every EDF.</p>\n'
    '""")\n\n',
    "",
    "remove section1bis banner")

# csv_suffix widget default + label
replace_in_cell(code,
    'csv_suffix = widgets.Text(value="_event_xml.csv", description="CSV suffix:",',
    'csv_suffix = widgets.Text(value="_ScoredEvents_Export.txt", description="Event export suffix:",',
    "csv_suffix widget")

# Remove section1bis widgets (run_check_button + out_check)
replace_in_cell(code,
    '# Section 1bis\n'
    'run_check_button = widgets.Button(description="Run CSV vs XML check", button_style="info", icon="check")\n'
    'out_check = widgets.Output()\n\n',
    "",
    "remove section1bis widgets")

# ===========================================================================
print("\n=== Cell 1: _update_info + run_scan ===")
# _update_info discovery + event detection
replace_in_cell(code,
    '        folder = Path(chooser.selected_path)\n'
    '        edfs = [f for f in folder.rglob("*") if f.suffix.lower() == ".edf" and not f.name.startswith("._")]\n'
    '        if not edfs:\n'
    '            existing_info.value = \'<small style="color:#888;">No EDF files found in selected folder.</small>\'\n'
    '            csv_suffix_info.value = ""\n'
    '            return\n'
    '        existing = load_existing_mapping(folder)\n'
    '        msg = f"<small>{len(edfs)} EDF file(s) found. "\n'
    '        msg += (f"{len(existing)} label(s) already mapped in event_remap.json."\n'
    '                if existing else "No existing event_remap.json yet.")\n'
    '        existing_info.value = msg + "</small>"\n'
    '        # --- Event-CSV suffix auto-detection (mirrors 5_quality_overview hypno-suffix block) ---\n'
    '        all_csv = [f for f in folder.rglob("*") if f.suffix.lower() == ".csv"]\n'
    '        suffix_counts = {}\n'
    '        for edf in edfs:\n'
    '            for csvf in all_csv:\n'
    '                if os.path.normcase(csvf.name).startswith(os.path.normcase(edf.stem)):\n'
    '                    suf = csvf.name[len(edf.stem):]\n'
    '                    suffix_counts[suf] = suffix_counts.get(suf, 0) + 1\n'
    '        if not suffix_counts:\n'
    '            csv_suffix_info.value = (\n'
    '                \'<small style="color:#e67e00;">No event CSV detected next to the EDFs \'\n'
    '                \'— set the suffix manually or rely on the XML fallback.</small>\')\n'
    '        else:\n'
    '            # Events have no "more specific remapped" variant (unlike hypnograms), so prefer the\n'
    '            # MOST FREQUENT suffix (shortest on ties).\n'
    '            best_suffix, best_count = max(\n'
    '                suffix_counts.items(), key=lambda x: (x[1], -len(x[0])))\n'
    '            csv_suffix.value = best_suffix\n'
    '            parts = [f\'<b>{s}</b>&nbsp;(×{c})\'\n'
    '                     for s, c in sorted(suffix_counts.items(), key=lambda x: -x[1])]\n'
    '            color = \'#2e7d32\' if best_count == len(edfs) else \'#e67e00\'\n'
    '            csv_suffix_info.value = (\n'
    '                f\'<small style="color:{color};">Detected:&nbsp;\'\n'
    '                f\'{"&nbsp;·&nbsp;".join(parts)}&nbsp;— \'\n'
    '                f\'{best_count}/{len(edfs)} EDF file(s) matched</small>\')\n',
    '        folder = Path(chooser.selected_path)\n'
    '        cdts = [f for f in folder.rglob("*") if f.suffix == ".cdt" and not f.name.startswith("._")]\n'
    '        if not cdts:\n'
    '            existing_info.value = \'<small style="color:#888;">No .cdt files found in selected folder.</small>\'\n'
    '            csv_suffix_info.value = ""\n'
    '            return\n'
    '        existing = load_existing_mapping(folder)\n'
    '        msg = f"<small>{len(cdts)} .cdt file(s) found. "\n'
    '        msg += (f"{len(existing)} label(s) already mapped in event_remap.json."\n'
    '                if existing else "No existing event_remap.json yet.")\n'
    '        existing_info.value = msg + "</small>"\n'
    '        # --- Event-export suffix auto-detection (Curry: *_ScoredEvents_Export.txt) ---\n'
    '        # Scan .txt files whose name contains \'event\' so hypnogram .txt files are not counted.\n'
    '        all_evt = [f for f in folder.rglob("*")\n'
    '                   if f.suffix.lower() == ".txt" and "event" in f.name.lower()]\n'
    '        suffix_counts = {}\n'
    '        for cdt in cdts:\n'
    '            for evtf in all_evt:\n'
    '                if os.path.normcase(evtf.name).startswith(os.path.normcase(cdt.stem)):\n'
    '                    suf = evtf.name[len(cdt.stem):]\n'
    '                    suffix_counts[suf] = suffix_counts.get(suf, 0) + 1\n'
    '        if not suffix_counts:\n'
    '            csv_suffix_info.value = (\n'
    '                \'<small style="color:#e67e00;">No event export detected next to the .cdt files \'\n'
    '                \'— set the suffix manually.</small>\')\n'
    '        else:\n'
    '            # Prefer the MOST FREQUENT suffix (shortest on ties).\n'
    '            best_suffix, best_count = max(\n'
    '                suffix_counts.items(), key=lambda x: (x[1], -len(x[0])))\n'
    '            csv_suffix.value = best_suffix\n'
    '            parts = [f\'<b>{s}</b>&nbsp;(×{c})\'\n'
    '                     for s, c in sorted(suffix_counts.items(), key=lambda x: -x[1])]\n'
    '            color = \'#2e7d32\' if best_count == len(cdts) else \'#e67e00\'\n'
    '            csv_suffix_info.value = (\n'
    '                f\'<small style="color:{color};">Detected:&nbsp;\'\n'
    '                f\'{"&nbsp;·&nbsp;".join(parts)}&nbsp;— \'\n'
    '                f\'{best_count}/{len(cdts)} .cdt file(s) matched</small>\')\n',
    "_update_info block")

# run_scan discovery
replace_in_cell(code,
    '        folder = Path(chooser.selected_path)\n'
    '        edfs = sorted(f for f in folder.rglob("*")\n'
    '                      if f.suffix.lower() == ".edf" and not f.name.startswith("._"))\n'
    '        if not edfs:\n'
    '            print("⚠️ No EDF files found in the selected folder.")\n'
    '            return\n',
    '        folder = Path(chooser.selected_path)\n'
    '        cdts = sorted(f for f in folder.rglob("*")\n'
    '                      if f.suffix == ".cdt" and not f.name.startswith("._"))\n'
    '        if not cdts:\n'
    '            print("⚠️ No .cdt files found in the selected folder.")\n'
    '            return\n',
    "run_scan discovery")

# run_scan loop head + failed reason
replace_in_cell(code,
    '        for edf in edfs:\n'
    '            fid = edf.stem\n'
    '            try:\n'
    '                events, source = load_events(edf, csv_suffix.value)\n'
    '            except Exception as e:\n'
    '                failed.append((fid, f"{type(e).__name__}: {e}"))\n'
    '                continue\n'
    '            if events is None:\n'
    '                failed.append((fid, "no readable event companion (.csv or .edf.XML)"))\n'
    '                continue\n',
    '        for cdt in cdts:\n'
    '            fid = cdt.stem\n'
    '            try:\n'
    '                events, source = load_events(cdt, csv_suffix.value)\n'
    '            except Exception as e:\n'
    '                failed.append((fid, f"{type(e).__name__}: {e}"))\n'
    '                continue\n'
    '            if events is None:\n'
    '                failed.append((fid, "no readable event export (*_ScoredEvents_Export.txt)"))\n'
    '                continue\n',
    "run_scan loop")

# run_scan "Scanned" summary line
replace_in_cell(code,
    'print(f"✅ Scanned {len(edfs)} EDF file(s): {n_with} with events, {len(failed)} failed."',
    'print(f"✅ Scanned {len(cdts)} .cdt file(s): {n_with} with events, {len(failed)} failed."',
    "run_scan summary")

# Remove the CSV/XML source-count print (single source in Curry)
replace_in_cell(code,
    '        n_csv = sum(1 for s in source_by_file.values() if s == "csv")\n'
    '        n_xml = sum(1 for s in source_by_file.values() if s == "xml")\n'
    '        print(f"   event source: {n_csv} from CSV, {n_xml} from XML fallback.")\n',
    "",
    "remove source-count print")

# ===========================================================================
print("\n=== Cell 1: remove run_csvxml_check + wiring + layout ===")
remove_between(code,
    "def run_csvxml_check(_=None):",
    "# ========================= Wiring & layout =========================",
    "remove run_csvxml_check fn")

replace_in_cell(code,
    "run_scan_button.on_click(run_scan)\n"
    "run_check_button.on_click(run_csvxml_check)\n"
    "preview_save_button.on_click(on_preview_save)",
    "run_scan_button.on_click(run_scan)\n"
    "preview_save_button.on_click(on_preview_save)",
    "remove run_check wiring")

replace_in_cell(code,
    "    section1bis, run_check_button, out_check,\n",
    "",
    "remove section1bis from layout")

# ===========================================================================
print("\n=== Validation ===")
try:
    out_str = json.dumps(nb, ensure_ascii=False, indent=1)
    json.loads(out_str)
    print(f"JSON valid. n_cells={len(nb['cells'])}")
except Exception as e:
    print(f"X JSON invalid: {e}")
    sys.exit(1)

for i, cell in enumerate(nb["cells"]):
    if cell["cell_type"] == "code":
        try:
            ast.parse(get_src(cell))
        except SyntaxError as e:
            print(f"X SyntaxError in code cell {i}: {e}")
            errors.append(f"syntax cell {i}")

if errors:
    print(f"\n! {len(errors)} issue(s): {errors}")

with open(DST, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write("\n")
print(f"Written: {DST}")
