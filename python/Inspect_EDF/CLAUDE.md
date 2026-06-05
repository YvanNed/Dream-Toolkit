# Dream-Toolkit — Inspect_EDF

For full project description, tool details, and planned modules, see [SPEC.md](SPEC.md).

## Key design decisions

- **Custom EDF binary parser for headers**: EDF headers are read using a hand-written binary parser rather than pyedflib. pyedflib was insufficiently robust for the heterogeneous EDF files encountered in real-world datasets (encoding edge cases, non-standard headers). This parser is used **only for header inspection** — it never loads signal data.
- **MNE for signal loading**: Whenever actual signal data needs to be loaded (e.g. for spectral analysis or future preprocessing), `mne.io.read_raw_edf()` is used. Do not use pyedflib for signal loading.
- **edfio for EDF export**: When writing EDF files (e.g. in `generate_test_data.py`), use `edfio.EdfSignal()` and `edfio.Edf()` directly — not `mne.export.export_raw()`. Direct edfio usage allows per-channel physical range overrides that MNE's export layer does not expose.
- **Dual delivery (Jupyter + Voila)**: Every user-facing tool exists in both forms. The Voila version hides all code cells; the Jupyter version is kept in sync for debugging. Always maintain both versions when modifying a tool.
- **TSV + HTML outputs**: TSV files are machine-readable for scripting; HTML reports are human-readable for non-programmers.
- **MNE duplicate channel suffixes**: MNE ≥ 1.8 renames duplicate channel names by appending `-0`, `-1`, etc. Use `drop_suffix_duplicates()` to keep only the `-0` variants, and `adapt_remap_dict_to_suffixes()` to match base channel names in remap dicts to MNE's suffixed names. Apply this pattern whenever loading an EDF with MNE and using a channel remap dict.
- **Physical bounds from MNE internals**: MNE 1.9 stores EDF physical range as `physical_max + offset` in `_raw_extras`, not as explicit `physical_min`/`physical_max`. Use `get_phys_bounds_uV()` (defined in `quality_overview_voila.ipynb`) to reconstruct these bounds when needed.

## Development context

- **Developer background**: The primary developer is a sleep research engineer with domain expertise in PSG/EEG but limited software engineering experience. Prioritize readability and correctness over abstraction or cleverness.
- **No unnecessary abstraction**: Keep functions explicit and readable. Three clear lines beat a clever one-liner. Do not introduce helper functions or classes unless the complexity clearly justifies it.
- **Error handling in notebooks and new features**: When implementing a new feature or creating a notebook, add `try/except` blocks proactively. Distinguish fatal vs non-fatal failures: a fatal step (e.g. file I/O, epoching) should add the item to a `failed` list and `continue` to the next iteration; a non-fatal step (e.g. re-referencing, report generation) should print a `⚠` warning and continue processing. In Voila notebooks specifically, every button callback and every per-item processing loop must be wrapped so that a single failure never silently crashes the whole run or freezes the UI. Always surface error messages to the user via a widget (e.g. `HTML`, `Output`) or `print()` inside a `widgets.Output` context.
- **Comments**: Add them only where the EEG/PSG domain logic would not be obvious to a reader unfamiliar with sleep research (e.g. why a 500 µV threshold flags clipping, or why stage 4 is remapped to N3).
- **No formal test suite**: Validate tools against real EDF files from the dataset. Do not add pytest or unittest infrastructure unless explicitly requested.
- **OS**: Windows 10 with Conda. PowerShell is the default shell. Use forward slashes in Python path strings (MNE and pathlib handle them correctly on Windows).
- **Running Python scripts from the terminal**: `conda` is NOT in the PATH. miniforge3 is installed at `$env:LOCALAPPDATA\miniforge3` (not `$env:USERPROFILE\miniforge3`). To run a script, always use the explicit executable path:
  ```powershell
  & "$env:LOCALAPPDATA\miniforge3\envs\inspect_edf\python.exe" tools/my_script.py
  ```
  Do NOT use `conda run`, `conda activate`, or bare `python` — they will all fail.
- **Updating environment.yml**: When adding a new `import` to any tool, add the corresponding package to `environment.yml` as part of the task — conda-forge section if available, pip section otherwise. Do not rely on transitive dependencies for packages that are directly imported.
- **Updating SPEC.md**: At the end of a task, propose updates to SPEC.md for any new feature, workflow change, constant, output file, or module detail — but do not write without explicit confirmation.
- **Updating CLAUDE.md**: If a design decision or development rule changes, update this file directly as part of the task.
- **Read SPEC.md at the start of each new task**: Before implementing anything, read `SPEC.md` once to get an up-to-date picture of the project structure, tool inventory, and planned modules. This avoids re-deriving context that is already documented there.
- **Ask clarifying questions before starting**: Before implementing any non-trivial task, ask the user clarifying questions to resolve ambiguities (format, scope, backward compatibility, etc.), then present a plan for validation before touching any file.
- **Language in tools and documentation**: All user-facing strings in notebooks (widget labels, print statements, HTML banners, section descriptions) must be in English. SPEC.md must also be written in English.

