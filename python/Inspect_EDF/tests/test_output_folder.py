"""Tools 5 and 7: the output folder is optional and defaults to the study folder; a different output
folder is allowed but warned about, since config_param/ (the exclusion registry) would be split in two.
Tool 7 also renames the former derivatives/raw_epo/ folder to derivatives/preproc_epo/ once, at Run.

No processing is run: the folders are empty, only the pickers and what they resolve to are checked.
"""
from pathlib import Path

import pytest

from nbdriver import NotebookSession
from pipeline import T5, T7


@pytest.fixture
def folders(tmp_path):
    study, elsewhere = tmp_path / 'my_study', tmp_path / 'elsewhere'
    study.mkdir()
    elsewhere.mkdir()
    return study, elsewhere


@pytest.mark.tool7
def test_tool7_output_folder_defaults_to_the_study_folder(folders):
    study, elsewhere = folders
    with NotebookSession(T7) as nb:
        nb.pick('fc_edf', study)
        assert Path(nb.value('str(out_root_path())')) == study
        assert Path(nb.value('str(registry_folder())')) == study
        assert nb.value('out_warn.value') == ''

        nb.pick('fc_out', elsewhere)
        assert Path(nb.value('str(out_root_path())')) == elsewhere
        assert 'not the study folder' in nb.value('out_warn.value')

        nb.pick('fc_out', study)
        assert nb.value('out_warn.value') == ''


@pytest.mark.tool7
def test_tool7_renames_the_former_raw_epo_folder(folders):
    """The epochs folder was called derivatives/raw_epo/ before: Run renames it once to preproc_epo/, so a
    study never has its participants split across two folders. When both exist, raw_epo/ is left alone."""
    study, _ = folders
    old, new = study / 'derivatives' / 'raw_epo', study / 'derivatives' / 'preproc_epo'
    old.mkdir(parents=True)
    (old / 'sub1_all-epo.fif').write_text('witness', encoding='utf-8')
    config = study / 'config.json'
    config.write_text('{}', encoding='utf-8')
    with NotebookSession(T7) as nb:
        nb.pick('fc_edf', study)
        nb.pick('fc_config', config)
        nb.click('btn_run')                        # empty study: nothing to process after the rename
        assert not old.exists()
        assert (new / 'sub1_all-epo.fif').read_text(encoding='utf-8') == 'witness'
        assert 'renamed to derivatives/preproc_epo/' in nb.widget_text()

        old.mkdir()
        nb.click('btn_run')
        assert old.is_dir() and new.is_dir()


@pytest.mark.tool5
def test_tool5_warns_when_the_output_folder_is_not_the_study_folder(folders):
    study, elsewhere = folders
    with NotebookSession(T5) as nb:
        nb.pick('fc_data', study)
        assert nb.value('lbl_out_warn.value') == ''

        nb.pick('fc_out', elsewhere)
        assert 'not the study' in nb.value('lbl_out_warn.value')

        nb.pick('fc_data', study)            # picking the study folder again resets the output picker
        assert nb.value('lbl_out_warn.value') == ''
