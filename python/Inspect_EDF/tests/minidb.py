"""Build the small test databases used by the pipeline tests. Two datasets, same layout:

'synthetic' (versioned golden, runs anywhere): four 4-h nights generated from fixed seeds by
    tests/synthetic.py: no real participant's data, nothing large stored in git.

'real' (local only, golden git-ignored): four full nights (4 channels at 256 Hz, ~934 epochs) copied
    from tools/test_data, derived from a real recording (73.edf). Pseudonymised headers, but the
    signal is a real person's PSG: these files and their golden snapshot never leave this machine.

Both carry the same four defects (injected by tools/generate_test_data.py):
    clipping          Fp1 clipped at 75 uV           -> many amplitude flags, Fp1 excluded by tool 6
    dead_channel      C3 scaled x0.01                -> dead channel, excluded by tool 6
    line_noise        50 Hz added on O1              -> notch / 1/f flags
    movement_bursts   high-amplitude bursts on Fp1   -> epoch rejection

The originals are never touched: everything is written into `root/data/`. The channel config (tool
2's remap_reref_persubject.json) is written here, so the tests do not depend on the state of
tools/test_data/config_param. Every recording uses the M2 (mastoid) reference of a classical clinical PSG
montage; the average reference (meant for high-density montages) gets a dedicated tool-7 test.
"""
import json
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
REAL_SRC = REPO_ROOT / 'tools' / 'test_data'

PARTICIPANTS = {
    'real': ['731_clipping', '734_dead_channel', '736_line_noise', '737_movement_bursts'],
    'synthetic': ['sim01_clipping', 'sim02_dead_channel', 'sim03_line_noise', 'sim04_movement_bursts'],
}
REAL_COMPANIONS = ['.edf', '.edf.XML', '_Hypnogram_Export.txt', '_Hypnogram_remapped.txt',
                   '_event_xml.csv']
REMAP = {'A2': 'M2', 'C3': 'C3', 'Fp1': 'Fp1', 'O1': 'O1'}
# tool 4's raw -> canonical event labels (the Compumedics labels of both datasets)
EVENT_REMAP = {
    'Arousal (ARO Limb)': 'arousal_limb', 'Arousal (ARO PLM)': 'arousal_limb',
    'Arousal (ARO RES)': 'arousal_respiratory', 'Arousal (ARO SPONT)': 'arousal_spontaneous',
    'Central Apnea': 'apnea_central', 'Hypopnea': 'hypopnea',
    'Limb Movement (Left)': 'limb_movement', 'Limb Movement (Right)': 'limb_movement',
    'Obstructive Apnea': 'apnea_obstructive', 'PLM (Left)': 'plm', 'PLM (Right)': 'plm',
    'SpO2 desaturation': 'spo2_desaturation',
}
AGES = [25, 70, 40, 55]          # the 2nd participant (dead channel) is the one '5bis age > 60' excludes


def manual_review_participant(dataset):
    """The participant tool 8 (manual review) is driven on."""
    return PARTICIPANTS[dataset][3]


def available(dataset):
    return dataset == 'synthetic' or (REAL_SRC / '731_clipping.edf').exists()


def build(root, dataset='synthetic'):
    """Create `root/data/` (recordings + config_param/ + participants.tsv), return its path."""
    data = Path(root) / 'data'
    cfg_dir = data / 'config_param'
    cfg_dir.mkdir(parents=True, exist_ok=True)
    pids = PARTICIPANTS[dataset]
    if dataset == 'real':
        for pid in pids:
            for suf in REAL_COMPANIONS:
                src = REAL_SRC / f'{pid}{suf}'
                if src.exists() and not (data / src.name).exists():
                    shutil.copy2(src, data / src.name)
    else:
        import synthetic
        for pid in pids:
            if not (data / f'{pid}.edf').exists():
                synthetic.write_recording(data, pid)
    config = {pid: {'config': 'config. 1', 'remap': dict(REMAP), 'ref_channels': ['M2']} for pid in pids}
    (cfg_dir / 'remap_reref_persubject.json').write_text(json.dumps(config, indent=2), encoding='utf-8')
    # A participant table: the real recordings all derive from the same night (identical sleep
    # metrics), so `age` is what lets a tool-5bis criterion exclude exactly one participant.
    (data / 'participants.tsv').write_text(
        'participant_id\tage\n' + ''.join(f'{pid}\t{age}\n' for pid, age in zip(pids, AGES)),
        encoding='utf-8')
    (cfg_dir / 'event_remap.json').write_text(json.dumps(EVENT_REMAP, indent=2), encoding='utf-8')
    return data
