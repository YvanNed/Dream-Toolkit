"""Drive a Voila notebook headlessly, the way a user would.

A real Jupyter kernel runs every cell of the notebook (so the code under test is exactly the Voila
code), then the test injects small "driver" cells that set widget values, pick files in the
FileChoosers and click buttons. Results are read back from the files the tool wrote on disk, and the
text printed by the tool (cell outputs + Output widgets) is collected so a test can assert on it.

Usage:
    with NotebookSession('tools/9_spectral_features_voila.ipynb') as nb:
        nb.pick('fc_clean', some_folder)
        nb.click('btn_scan')
        n = nb.value('len(S["parts"])')
"""
import json
import os
import re
from contextlib import ExitStack
from pathlib import Path

import nbformat
from nbclient import NotebookClient

REPO_ROOT = Path(__file__).resolve().parent.parent      # Inspect_EDF/ = the canonical launch cwd

# Helpers defined inside the kernel once the notebook has run.
_PRELUDE = r'''
import os as _dtk_os, json as _dtk_json, re as _dtk_re

def _dtk_pick(fc, path):
    """Select a file or a folder in an ipyfilechooser FileChooser, then fire its callback,
    exactly like opening the dialog and clicking its Select button."""
    path = _dtk_os.path.abspath(str(path))
    if _dtk_os.path.isdir(path):
        folder, name = path, ''
    else:
        folder, name = _dtk_os.path.split(path)
    fc.reset(path=folder, filename=name)
    fc._show_dialog()
    fc._on_select_click(None)        # dialog open -> apply the selection + run the callback

def _dtk_widget_text():
    """Text of every Output widget (where the tools print their run logs)."""
    import ipywidgets as _w
    chunks = []
    for w in list(_w.Widget.widgets.values()):
        if isinstance(w, _w.Output):
            for o in w.outputs:
                if 'text' in o:
                    chunks.append(o['text'])
                elif 'data' in o and 'text/html' in o['data']:
                    # the tools log a lot through display(HTML(...)): keep the text, drop the tags
                    chunks.append(_dtk_re.sub(r'<[^>]+>', ' ', o['data']['text/html']))
                elif 'data' in o and 'text/plain' in o['data']:
                    chunks.append(o['data']['text/plain'])
                elif o.get('output_type') == 'error':
                    chunks.append('\n'.join(o.get('traceback', [])))
    return '\n'.join(chunks)
'''

_MARK = '@@DTK_JSON@@'


class DriverError(RuntimeError):
    """A driver cell raised, or a widget callback printed a traceback."""


class NotebookSession:
    def __init__(self, notebook, cwd=REPO_ROOT, run_all=True):
        self.path = (REPO_ROOT / notebook) if not Path(notebook).is_absolute() else Path(notebook)
        self.cwd = Path(cwd)
        self.run_all_cells = run_all
        self.log = []                         # every text output, in order (for debugging a failure)

    # -- lifecycle ---------------------------------------------------------
    def __enter__(self):
        self.nb = nbformat.read(str(self.path), as_version=4)
        self.client = NotebookClient(self.nb, timeout=None, kernel_name='python3',
                                     resources={'metadata': {'path': str(self.cwd)}})
        self._stack = ExitStack()
        self._stack.enter_context(self.client.setup_kernel(cwd=str(self.cwd)))
        if self.run_all_cells:
            for i, cell in enumerate(list(self.nb.cells)):
                if cell.cell_type == 'code':
                    self._exec_cell(cell, i)
        self.run(_PRELUDE)
        return self

    def __exit__(self, *exc):
        self._stack.close()
        return False

    # -- execution ---------------------------------------------------------
    def _exec_cell(self, cell, index):
        try:
            self.client.execute_cell(cell, index)
        except Exception as e:                       # CellExecutionError: re-raise with the log
            raise DriverError(f'{self.path.name} cell {index} failed:\n{e}') from None
        text = _outputs_text(cell.get('outputs', []))
        self.log.append(text)
        return text

    def run(self, code):
        """Execute a driver cell; return its printed text. Raises when the cell itself errors, or
        when a widget callback it triggered printed a Python traceback (ipywidgets swallows callback
        exceptions and only prints them)."""
        cell = nbformat.v4.new_code_cell(code)
        self.nb.cells.append(cell)
        text = self._exec_cell(cell, len(self.nb.cells) - 1)
        if 'Traceback (most recent call last)' in text:
            raise DriverError(f'callback traceback while running:\n{code}\n---\n{text}')
        return text

    # -- convenience -------------------------------------------------------
    def pick(self, chooser, path):
        return self.run(f'_dtk_pick({chooser}, r"{Path(path)}")')

    def set(self, widget, value):
        return self.run(f'{widget}.value = {value!r}')

    def click(self, button):
        return self.run(f'{button}.click()')

    def value(self, expr):
        """Evaluate an expression in the kernel; JSON round-trip (non-JSON objects as str)."""
        out = self.run(f'print({_MARK!r} + _dtk_json.dumps({expr}, default=str))')
        line = [l for l in out.splitlines() if l.startswith(_MARK)][-1]
        return json.loads(line[len(_MARK):])

    def widget_text(self):
        """Everything the tool printed into its Output widgets so far."""
        return self.value('_dtk_widget_text()')


def _outputs_text(outputs):
    chunks = []
    for o in outputs:
        if o.get('output_type') == 'stream':
            chunks.append(o.get('text', ''))
        elif o.get('output_type') == 'error':
            chunks.append('\n'.join(o.get('traceback', [])))
        elif 'data' in o and 'text/html' in o['data']:
            chunks.append(re.sub(r'<[^>]+>', ' ', o['data']['text/html']))
        elif 'data' in o and 'text/plain' in o['data']:
            chunks.append(o['data']['text/plain'])
    return ''.join(c if c.endswith('\n') else c + '\n' for c in chunks)
