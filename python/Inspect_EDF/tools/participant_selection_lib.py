"""Shared participant handling: the exclusion registry and the participant selector widget.

A deliberate exception to the toolkit's "duplicate small helpers" rule (like qc_rejected_epochs_lib):
this code is used by tools 5, 5bis, 6, 7, 8, 8bis, 9 and their Curry twins, and seven copies would drift.

1. Exclusion registry: `<data_folder>/config_param/participant_exclusions.tsv`
   One row per decision, written by the tool that took it:
       file_id | source | decision | reason | comment | decided_at
   source   = 5bis_criteria / 5bis_manual / 7_manual / 8_manual / 8bis_auto
   decision = exclude, or include (a FORCED inclusion, written by 5bis's manual editor only)
   Each tool rewrites only its own rows. The effective status of a participant is:
       a forced inclusion exists         -> included (it overrides every exclusion)
       else at least one exclusion       -> excluded, reasons joined with their source
       else                              -> included
   Excluded participants are still processed by every tool; they are only flagged in the tables
   (`participant_excluded`, `participant_exclude_reason`) and left out of group statistics / figures.

2. ParticipantSelector: the "All participants / Subset only" block shared by the processing tools.
"""
import datetime
import os
import re
from pathlib import Path

import pandas as pd

REGISTRY_NAME = 'participant_exclusions.tsv'
REGISTRY_COLUMNS = ['file_id', 'source', 'decision', 'reason', 'comment', 'decided_at']
SOURCES = ('5bis_criteria', '5bis_manual', '7_manual', '8_manual', '8bis_auto')
EXCLUDE, INCLUDE = 'exclude', 'include'
COL_EXCLUDED = 'participant_excluded'
COL_REASON = 'participant_exclude_reason'


def id_key(file_id):
    """Comparison key for a participant id (normcase: Windows ids differ only by case / slashes)."""
    return os.path.normcase(str(file_id).strip())


# ---------------------------------------------------------------------------
# Locating the data folder
# ---------------------------------------------------------------------------
def registry_path(data_folder):
    return Path(data_folder) / 'config_param' / REGISTRY_NAME


def find_data_folder(path):
    """The data folder a path belongs to: the first ancestor (or the path itself) that holds a
    `config_param/` folder, else the parent of the `derivatives/` or `reports_*` folder the path lies
    in, else None. Lets a tool that only knows e.g. `<data>/derivatives/raw_epo` find the registry."""
    if path is None:
        return None
    p = Path(path)
    if p.is_file():
        p = p.parent
    for anc in [p] + list(p.parents):
        if (anc / 'config_param').is_dir():
            return anc
    for anc in [p] + list(p.parents):
        if anc.name == 'derivatives' or anc.name.startswith('reports_'):
            return anc.parent
    return None


# ---------------------------------------------------------------------------
# Reading / writing the registry
# ---------------------------------------------------------------------------
def empty_registry():
    return pd.DataFrame(columns=REGISTRY_COLUMNS)


def load_registry(data_folder):
    """(registry DataFrame, warning or None). Absent file -> empty registry, no warning (nobody
    excluded yet). Unreadable file -> empty registry + a warning: every participant is then treated
    as included, which the calling tool must show to the user."""
    if data_folder is None:
        return empty_registry(), None
    path = registry_path(data_folder)
    if not path.is_file():
        return empty_registry(), None
    try:
        df = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    except Exception as e:
        return empty_registry(), f'exclusion registry unreadable ({path}): {e}. Everyone treated as included.'
    missing = [c for c in ('file_id', 'source', 'decision') if c not in df.columns]
    if missing:
        return empty_registry(), (f'exclusion registry {path} has no {", ".join(missing)} column: '
                                  f'everyone treated as included.')
    for c in REGISTRY_COLUMNS:
        if c not in df.columns:
            df[c] = ''
    df = df[REGISTRY_COLUMNS + [c for c in df.columns if c not in REGISTRY_COLUMNS]]
    df['decision'] = df['decision'].str.strip().str.lower()
    return df.reset_index(drop=True), None


def _write_atomic(df, path):
    """Write to a temporary file, then rename it over the registry: an interruption never leaves a
    half-written registry behind (os.replace is atomic on the same volume)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    df.to_csv(tmp, sep='\t', index=False)
    os.replace(tmp, path)


def write_registry_rows(data_folder, source, rows, replace_ids=None):
    """Replace this tool's rows in the registry and return the registry written.

    source      : who decides (one of SOURCES)
    rows        : list of dicts with file_id, decision ('exclude' / 'include'), reason, comment
    replace_ids : the participants whose `source` rows are replaced by `rows`; None = ALL the rows
                  of `source` (e.g. a new 5bis validation replaces the previous one entirely).
    Rows of the other sources are never touched. Raises when the existing registry is unreadable
    (refusing to overwrite decisions we could not read)."""
    if source not in SOURCES:
        raise ValueError(f'unknown registry source {source!r}')
    path = registry_path(data_folder)
    reg, warn = load_registry(data_folder)
    if warn:
        raise RuntimeError(warn.replace('Everyone treated as included.', 'Not overwritten.'))
    mine = reg['source'] == source
    if replace_ids is not None:
        ids = {id_key(i) for i in replace_ids}
        mine &= reg['file_id'].map(id_key).isin(ids)
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    new = pd.DataFrame([{'file_id': str(r['file_id']), 'source': source,
                         'decision': r.get('decision', EXCLUDE), 'reason': r.get('reason', ''),
                         'comment': r.get('comment', ''), 'decided_at': r.get('decided_at') or now}
                        for r in rows], columns=REGISTRY_COLUMNS)
    for r in new.itertuples():
        if r.decision not in (EXCLUDE, INCLUDE):
            raise ValueError(f'decision must be exclude or include, got {r.decision!r}')
        if r.decision == INCLUDE and source != '5bis_manual':
            raise ValueError('a forced inclusion can only be written by 5bis_manual')
    out = pd.concat([reg[~mine], new], ignore_index=True)
    out = out.sort_values(['file_id', 'source', 'decided_at'], kind='stable').reset_index(drop=True)
    _write_atomic(out, path)
    return out


def remove_registry_rows(data_folder, source, file_ids):
    """Delete this source's rows for these participants (e.g. 'Remove exclusion' in tool 8)."""
    return write_registry_rows(data_folder, source, [], replace_ids=file_ids)


# ---------------------------------------------------------------------------
# Effective status
# ---------------------------------------------------------------------------
def effective_status(registry):
    """{normcase(file_id): {'file_id', 'excluded', 'reason', 'forced_include'}} for every participant
    that has at least one row. A participant with no row is included (use `status_of`)."""
    out = {}
    if registry is None or len(registry) == 0:
        return out
    for key, grp in registry.groupby(registry['file_id'].map(id_key), sort=False):
        forced = grp[grp['decision'] == INCLUDE]
        excl = grp[grp['decision'] == EXCLUDE]
        reasons = ' | '.join(f'{r.source}: {r.reason}' if r.reason else r.source for r in excl.itertuples())
        out[key] = {'file_id': grp['file_id'].iloc[0],
                    'excluded': len(excl) > 0 and len(forced) == 0,
                    'reason': reasons if len(forced) == 0 else '',
                    'forced_include': len(forced) > 0,
                    'overridden_reason': reasons if len(forced) > 0 else ''}
    return out


def status_of(status, file_id):
    """(excluded, reason) of one participant; (False, '') when the registry says nothing about it."""
    s = status.get(id_key(file_id))
    return (s['excluded'], s['reason']) if s else (False, '')


def excluded_ids(status):
    return sorted(s['file_id'] for s in status.values() if s['excluded'])


def add_exclusion_columns(df, status, id_col='file_id'):
    """Copy of `df` with participant_excluded (bool) and participant_exclude_reason (str) appended
    (replaced when already present, so a rebuilt table always reflects the CURRENT registry)."""
    out = df.drop(columns=[c for c in (COL_EXCLUDED, COL_REASON) if c in df.columns]).copy()
    pairs = [status_of(status, fid) for fid in out[id_col]]
    out[COL_EXCLUDED] = [p[0] for p in pairs]
    out[COL_REASON] = [p[1] for p in pairs]
    return out


def excluded_banner_html(reason):
    """Red banner for the top of an excluded participant's own report."""
    import html as _html
    return ('<div style="background:#fde8e8;border:2px solid #d03b3b;color:#8b1a1a;padding:8px 12px;'
            'margin:6px 0;border-radius:4px;font-weight:bold">EXCLUDED participant: '
            f'{_html.escape(reason)}<br><span style="font-weight:normal;font-size:90%">Still computed so '
            'the decision can be revisited; left out of the group statistics and figures.</span></div>')


# ---------------------------------------------------------------------------
# Participant selector (All / Subset + skip), shared by the processing tools
# ---------------------------------------------------------------------------
def parse_id_list(text):
    """'165, 176;180\\n181' -> ['165', '176', '180', '181'] (comma / semicolon / space / newline)."""
    return [t for t in re.split(r'[,;\s]+', text or '') if t]


class ParticipantSelector:
    """'All participants / Subset only' + 'Skip already processed', with a summary line and a
    prominent warning box. Usage in a tool:

        sel = ParticipantSelector()
        display(sel.widget)
        sel.set_participants(ids, done={...}, excluded={id: reason})   # after the scan
        to_run = sel.selected_ids()                                      # at Run
    """

    def __init__(self, skip_default=True):
        import ipywidgets as w
        self._w = w
        self._ids = []                     # all participant ids, display order
        self._done = set()                 # normcase keys of the already-processed ones
        self._excluded = {}                # normcase key -> reason
        self.dd_mode = w.ToggleButtons(options=['All participants', 'Subset only'],
                                       value='All participants', style={'button_width': '150px'})
        self.tags = w.TagsInput(value=[], allowed_tags=[], allow_duplicates=False,
                                layout=w.Layout(width='620px'))
        self.txt_paste = w.Text(placeholder='or paste a list: 165, 176 ...',
                                layout=w.Layout(width='300px'))
        self.btn_paste = w.Button(description='Add list', icon='plus', layout=w.Layout(width='100px'))
        self.btn_clear = w.Button(description='Clear', icon='times', layout=w.Layout(width='80px'))
        self.tgl_list = w.ToggleButton(value=False, description='Show participant list', icon='list',
                                       layout=w.Layout(width='190px'))
        self.lst = w.Select(options=[], rows=10, layout=w.Layout(width='620px', display='none'))
        self.cb_skip = w.Checkbox(value=skip_default, indent=False,
                                  description='Skip already processed participants',
                                  layout=w.Layout(width='420px'))
        self.lbl_summary = w.HTML('<i>Scan first.</i>')
        self.html_warn = w.HTML('')
        self.lbl_paste = w.HTML('')
        self.box_subset = w.VBox([w.HBox([self.tags]),
                                  w.HBox([self.txt_paste, self.btn_paste, self.btn_clear, self.lbl_paste])],
                                 layout=w.Layout(display='none'))
        self.widget = w.VBox([w.HBox([self.dd_mode, self.tgl_list]), self.box_subset, self.lst,
                              self.cb_skip, self.lbl_summary, self.html_warn])
        self.dd_mode.observe(self._on_mode, names='value')
        self.tags.observe(self._refresh, names='value')
        self.cb_skip.observe(self._refresh, names='value')
        self.tgl_list.observe(self._on_toggle_list, names='value')
        self.lst.observe(self._on_list_pick, names='value')
        self.btn_paste.on_click(self._on_paste)
        self.btn_clear.on_click(lambda _: setattr(self.tags, 'value', []))

    # -- data -----------------------------------------------------------------
    def set_participants(self, ids, done=(), excluded=None):
        """ids: every participant found by the scan; done: the already-processed ones;
        excluded: {file_id: reason} from the registry (they stay selectable: still processed)."""
        self._ids = [str(i) for i in ids]
        self._done = {id_key(i) for i in done}
        self._excluded = {id_key(k): v for k, v in (excluded or {}).items()}
        known = {id_key(i) for i in self._ids}
        self.tags.allowed_tags = list(self._ids)
        self.tags.value = [t for t in self.tags.value if id_key(t) in known]
        self._fill_list()
        self._refresh()

    def _status_text(self, fid):
        k = id_key(fid)
        parts = ['processed' if k in self._done else 'not processed']
        if k in self._excluded:
            parts.append(f'EXCLUDED ({self._excluded[k]})')
        return ' - '.join(parts)

    def _fill_list(self):
        width = max([len(i) for i in self._ids] + [8])
        self.lst.options = [(f'{fid.ljust(width)}   {self._status_text(fid)}', fid) for fid in self._ids]
        self.lst.value = None

    # -- selection ------------------------------------------------------------
    @property
    def subset_mode(self):
        return self.dd_mode.value == 'Subset only'

    def chosen_ids(self):
        """The participants chosen, before the skip rule (all of them, or the subset in id order)."""
        if not self.subset_mode:
            return list(self._ids)
        chosen = {id_key(t) for t in self.tags.value}
        return [fid for fid in self._ids if id_key(fid) in chosen]

    def skipped_ids(self):
        """Chosen participants that the skip rule leaves out (already processed, skip ticked)."""
        if not self.cb_skip.value:
            return []
        return [fid for fid in self.chosen_ids() if id_key(fid) in self._done]

    def selected_ids(self):
        """The participants to run now."""
        skipped = {id_key(f) for f in self.skipped_ids()}
        return [fid for fid in self.chosen_ids() if id_key(fid) not in skipped]

    def warning_html(self):
        """The prominent box shown when a subset contains already-processed participants that the
        skip rule will leave out ('' otherwise). Also printed at the top of the run log."""
        if not self.subset_mode:
            return ''
        skipped = self.skipped_ids()
        if not skipped:
            return ''
        names = ', '.join(skipped[:20]) + (f' (+{len(skipped) - 20} more)' if len(skipped) > 20 else '')
        return ('<div style="background:#fff4e5;border:2px solid #ef6c00;color:#7a3e00;padding:8px 12px;'
                'margin:6px 0;border-radius:4px"><b>&#9888; '
                f'{len(skipped)} of the {len(self.chosen_ids())} selected participant(s) will NOT be '
                f'processed: already processed</b> ({names}).<br>Untick <b>Skip already processed '
                'participants</b> to reprocess them.</div>')

    # -- callbacks ----------------------------------------------------------------
    def _on_mode(self, change=None):
        self.box_subset.layout.display = '' if self.subset_mode else 'none'
        self._refresh()

    def _on_toggle_list(self, change):
        self.lst.layout.display = '' if change['new'] else 'none'
        self.tgl_list.description = 'Hide participant list' if change['new'] else 'Show participant list'

    def _on_list_pick(self, change):
        fid = change['new']
        if fid is None:
            return
        if fid not in self.tags.value:
            self.tags.value = list(self.tags.value) + [fid]
        if not self.subset_mode:
            self.dd_mode.value = 'Subset only'
        self.lst.value = None

    def _on_paste(self, _=None):
        known = {id_key(i): i for i in self._ids}
        wanted = parse_id_list(self.txt_paste.value)
        found = [known[id_key(t)] for t in wanted if id_key(t) in known]
        unknown = [t for t in wanted if id_key(t) not in known]
        self.tags.value = list(self.tags.value) + [f for f in found if f not in self.tags.value]
        if found and not self.subset_mode:
            self.dd_mode.value = 'Subset only'
        self.txt_paste.value = ''
        self.lbl_paste.value = (f'<span style="color:#c62828">unknown id(s): {", ".join(unknown)}</span>'
                                if unknown else '')

    def _refresh(self, change=None):
        n = len(self._ids)
        if not n:
            self.lbl_summary.value = '<i>Scan first.</i>'
            self.html_warn.value = ''
            return
        n_done = sum(1 for i in self._ids if id_key(i) in self._done)
        n_excl = sum(1 for i in self._ids if id_key(i) in self._excluded)
        chosen = self.chosen_ids()
        run = self.selected_ids()
        excl_txt = (f', <b>{n_excl}</b> excluded (still computed, flagged in the tables)' if n_excl else '')
        if self.subset_mode:
            head = f'Subset: <b>{len(chosen)}</b> of {n} participants'
        else:
            head = f'<b>{n}</b> participants, <b>{n_done}</b> already processed{excl_txt}'
        self.lbl_summary.value = f'{head} &rarr; <b>{len(run)}</b> will run.'
        self.html_warn.value = self.warning_html()
