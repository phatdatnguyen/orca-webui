"""Regression tests for the Gradio event handlers.

These modules import Gradio (and, for the working directory, nglview) at import time,
so the whole module is skipped when those optional UI deps are absent. Run them in the
full application environment (the venv) to exercise the fixed bugs end-to-end.
"""
import os

import pytest


class FakeProgress:
    """Stand-in for gr.Progress(): .tqdm() just passes the iterable through."""
    def tqdm(self, iterable, total=None, desc=None):
        return iterable


# --- working_directory.on_file_list_change: chronological sort ---------------

def test_file_list_sorted_by_mtime_newest_first(tmp_path):
    pytest.importorskip("gradio")
    pytest.importorskip("nglview")
    from orca_webui import working_directory

    # Create three files and stamp explicit, increasing modification times.
    names = ["old.xyz", "middle.inp", "new.log"]
    for i, name in enumerate(names):
        p = tmp_path / name
        p.write_text("data")
        mtime = 1_000_000 + i * 1000  # old < middle < new
        os.utime(p, (mtime, mtime))

    df = working_directory.on_file_list_change(str(tmp_path))
    ordered = list(df["File"])
    # Newest first: the ctime-string sort bug would have ordered by weekday name.
    assert ordered == ["new.log", "middle.inp", "old.xyz"]


# --- conformer_generation.on_generate_conformers: fixed count/invalid SMILES -

def test_generate_conformers_invalid_smiles_returns_error(tmp_path):
    pytest.importorskip("gradio")
    pytest.importorskip("gradio_molecule2d")
    from orca_webui import conformer_generation

    status, _files = conformer_generation.on_generate_conformers(
        str(tmp_path), "not_a_valid_smiles", 0, 1, 3, "conformer", "xyz",
        progress=FakeProgress())
    assert "invalid SMILES" in status
    assert "color:red" in status
    assert not list(tmp_path.glob("conformer*"))


def test_generate_conformers_writes_expected_files(tmp_path):
    pytest.importorskip("gradio")
    pytest.importorskip("gradio_molecule2d")
    from orca_webui import conformer_generation

    status, _files = conformer_generation.on_generate_conformers(
        str(tmp_path), "CCO", 0, 1, 3, "conformer", "xyz",
        progress=FakeProgress())
    assert "color:green" in status
    produced = sorted(p.name for p in tmp_path.glob("conformer_*.xyz"))
    # RDKit may embed fewer than requested; every produced file is named
    # sequentially and no GetConformer(bad_id) crash occurs (the fixed bug).
    assert len(produced) >= 1
    assert produced[0] == "conformer_1.xyz"
    for name in produced:
        assert (tmp_path / name).read_text().strip() != ""
