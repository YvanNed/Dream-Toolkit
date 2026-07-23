"""
One-shot script to generate 3_remap_hypno_curry_voila.ipynb from the EDF original.
Run from the Inspect_EDF root:
    & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" tools_curry/_make_tool3_curry.py
"""
import json
import os

SRC  = "tools/3_remap_hypno_voila.ipynb"
DST  = "tools_curry/3_remap_hypno_curry_voila.ipynb"

replacements = [
    # File discovery anchor: bare .edf → bare .cdt
    (
        "if f.suffix.lower() == '.edf' and not f.name.startswith('._')",
        "if f.suffix == '.cdt' and not f.name.startswith('._')"
    ),
    # The message "No EDF files" → Curry.
    # NB: unlike _make_tool{5,6}_curry.py, this script string-replaces the RAW notebook JSON
    # (see the open() below), where a newline is the two characters \n and a quote is \" — so a
    # pattern must never contain a real newline or a bare ". Keep every pattern here a plain,
    # escape-free substring. Occurs once in the source.
    (
        "No EDF files found in selected folder",
        "No .cdt files found in selected folder"
    ),
    # Suffix match: txt.name starts with edf.stem → cdt.stem
    # The stem for .cdt files: Path('y_S005.cdt').stem == 'y_S005', which is correct.
    # No code change needed — stem() already works on .cdt paths.

    # Title markdown
    (
        "# Remap hypnogram labels",
        "# Remap hypnogram labels — Curry 9 (.cdt)"
    ),
    # Match message "EDF files matched" → ".cdt files matched"
    (
        "f'&nbsp;— {best_count}/{n_total} EDF files matched</small>'",
        "f'&nbsp;— {best_count}/{n_total} .cdt files matched</small>'"
    ),
    # chooser title hint
    (
        "Choose your data folder",
        "Choose your data folder (containing .cdt files)"
    ),
]

with open(SRC, encoding="utf-8") as f:
    content = f.read()

n_applied = 0
for old, new in replacements:
    if old in content:
        content = content.replace(old, new, 1)
        n_applied += 1
        print(f"  ✓ replaced: {repr(old[:60])}")
    else:
        print(f"  ⚠ NOT FOUND (may already be changed or pattern mismatch): {repr(old[:60])}")

# Validate JSON after replacements
try:
    nb = json.loads(content)
    print(f"\nJSON valid. {n_applied}/{len(replacements)} replacements applied.")
except json.JSONDecodeError as e:
    print(f"\n✗ JSON invalid after replacements: {e}")
    raise

with open(DST, "w", encoding="utf-8") as f:
    f.write(content)

print(f"Written: {DST}")
