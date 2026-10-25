"""Regressions for workspace data preservation and conformer generation."""
import os
from pathlib import Path
import subprocess
import sys

import pytest

pytest.importorskip("gradio")
pytest.importorskip("nglview")

from orca_webui import conformer_generation, working_directory


class Progress:
    def __call__(self, fraction, desc=None):
        pass

    def tqdm(self, iterable, total=None, desc=None):
        return iterable


def generate(directory, **overrides):
    values = dict(working_directory_path=str(directory) if directory else None,
                  input_smiles="CCO", charge=0, multiplicity=1, num_confs=1,
                  energy_threshold=0.1, rms_threshold=0.5, file_name="conformer",
                  file_type="xyz", progress=Progress())
    values.update(overrides)
    return conformer_generation.on_generate_conformers(**values)


def test_text_editor_saves_loaded_file_and_rejects_directory_switch(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    (first / "a.inp").write_text("original a")
    (first / "b.inp").write_text("original b")
    (second / "a.inp").write_text("other directory")

    _text, _save, loaded_path = working_directory.on_view_text_file(str(first), "a.inp")
    working_directory.on_save_viewed_text_file(str(first), loaded_path, "edited a")
    assert (first / "a.inp").read_text() == "edited a"
    assert (first / "b.inp").read_text() == "original b"

    working_directory.on_save_viewed_text_file(str(second), loaded_path, "stale editor")
    assert (second / "a.inp").read_text() == "other directory"


def test_text_save_failure_does_not_truncate_existing_file(tmp_path):
    target = tmp_path / "job.inp"
    target.write_text("precious input")
    working_directory.on_save_text_file(str(tmp_path), "job.inp", None)
    assert target.read_text() == "precious input"
    assert [p.name for p in tmp_path.iterdir()] == ["job.inp"]


def test_upload_does_not_replace_existing_structure(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    (workspace / "molecule.xyz").write_text("existing geometry")
    uploaded = tmp_path / "molecule.xyz"
    uploaded.write_text("new geometry")

    files = working_directory.on_upload_file(str(workspace), str(uploaded))
    assert files == ["molecule.xyz"]
    assert (workspace / "molecule.xyz").read_text() == "existing geometry"


@pytest.mark.parametrize("directory_name", [None, "", "../outside", "/outside", "..", r"..\outside"])
def test_invalid_directory_open_preserves_current_state(tmp_path, directory_name):
    updates = working_directory.on_open_working_directory(directory_name, base_path=tmp_path / "data")
    assert len(updates) == 5
    assert all(update == {"__type__": "update"} for update in updates)
    assert not list(tmp_path.iterdir())


def test_workspace_open_uses_requested_data_root(tmp_path):
    data_root = tmp_path / "alternate" / "data"
    dropdown, directory, files, _upload, _clean = working_directory.on_open_working_directory("experiment", base_path=data_root)
    assert directory == str(data_root / "experiment")
    assert Path(directory).is_dir()
    assert dropdown["choices"] == ["experiment"]
    assert files == []


def test_cannot_delete_or_save_outside_workspace(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    outside = tmp_path / "outside.inp"
    outside.write_text("keep")
    working_directory.on_delete_file(str(workspace), "../outside.inp")
    working_directory.on_save_text_file(str(workspace), "../outside.inp", "overwrite")
    assert outside.read_text() == "keep"


def test_clean_without_directory_does_not_clean_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    scratch = tmp_path / "important.tmp"
    scratch.write_text("keep")
    assert working_directory.on_clean_working_directory(None) == []
    assert scratch.read_text() == "keep"


def test_startup_cleanup_preserves_logs_and_external_files(tmp_path):
    from orca_webui.app import _cleanup_transient_files

    (tmp_path / "calculation.log").write_text("ORCA result")
    static = tmp_path / "static"
    static.mkdir()
    (static / "structure.html").write_text("old viewer")
    outside = tmp_path / "outside.html"
    outside.write_text("keep")
    (static / "linked.html").symlink_to(outside)

    _cleanup_transient_files(tmp_path)

    assert (tmp_path / "calculation.log").read_text() == "ORCA result"
    assert not (static / "structure.html").exists()
    assert outside.read_text() == "keep"


def test_structure_viewers_do_not_overwrite_each_other_or_follow_old_symlink(tmp_path, monkeypatch):
    static = tmp_path / "static"
    static.mkdir()
    outside = tmp_path / "outside.html"
    outside.write_text("keep")
    (static / "structure.html").symlink_to(outside)
    structure = tmp_path / "hydrogen.xyz"
    structure.write_text("0 1\nH 0 0 0\nH 0 0 0.75\n")
    monkeypatch.setattr(working_directory.nglview, "write_html",
                        lambda filename, views: Path(filename).write_text("viewer"))

    first = working_directory.on_view_structure_file(str(tmp_path), structure.name, static_dir=static)
    second = working_directory.on_view_structure_file(str(tmp_path), structure.name, static_dir=static)

    assert first and second and first != second
    assert len(list(static.glob("structure_*.html"))) == 2
    assert outside.read_text() == "keep"


def test_app_builds_with_bundled_editor_and_alternate_root(tmp_path):
    # Build from the source checkout, with runtime files outside the launch directory.
    script = """
import sys
sys.path.insert(0, sys.argv[2])
from pathlib import Path
from fastapi.staticfiles import StaticFiles
from orca_webui import app, molecular_editor
original = app.working_directory_blocks
def checked_blocks(*, base_path, static_dir):
    assert base_path == Path(sys.argv[1]) / 'data'
    assert static_dir == Path(sys.argv[1]) / 'static'
    return original(base_path=base_path, static_dir=static_dir)
app.working_directory_blocks = checked_blocks
application = app.build_app(Path(sys.argv[1]))
editor_mount = next(route for route in application.routes if route.path == '/molecule-editor')
assert isinstance(editor_mount.app, StaticFiles)
assert Path(editor_mount.app.directory) == molecular_editor.KETCHER_DIR
assert (molecular_editor.KETCHER_DIR / 'index.html').is_file()
gradio_mount = application.routes[-1]
assert gradio_mount.path == ''
assert application.routes.index(editor_mount) < application.routes.index(gradio_mount)
blocks = gradio_mount.app.get_blocks()
assert app._STYLES_PATH.read_text() in blocks.css
editor = next(block for block in blocks.blocks.values()
              if getattr(block, 'elem_id', None) == 'molecular-editor')
assert editor.value == {'smiles': '', 'load_id': None}
assert callable(molecular_editor.smiles_from_editor)
"""
    run_root = tmp_path / "run"
    source_root = Path(__file__).resolve().parents[1] / "src"
    completed = subprocess.run([sys.executable, "-c", script, str(run_root), str(source_root)],
                               capture_output=True, text=True, timeout=60,
                               cwd=tmp_path,
                               env={**os.environ, "GRADIO_ANALYTICS_ENABLED": "False"})
    assert completed.returncode == 0, completed.stderr
    assert (run_root / "data").is_dir()
    assert (run_root / "static").is_dir()
    assert not (tmp_path / "data").exists()


def test_conformer_generation_preserves_existing_files_before_embedding(tmp_path, monkeypatch):
    previous = tmp_path / "conformer_1.xyz"
    previous.write_text("saved geometry")
    monkeypatch.setattr(conformer_generation, "generate_unique_conformers",
                        lambda *args: pytest.fail("Should reject the name before embedding"))
    status, files, table = generate(tmp_path)
    assert "color:red" in status and "already exists" in status
    assert previous.read_text() == "saved geometry"
    assert files == [previous.name]
    assert table.empty


@pytest.mark.parametrize("overrides", [
    {"input_smiles": ""}, {"input_smiles": None}, {"num_confs": 0},
    {"num_confs": 1.5}, {"file_name": "../outside"}, {"file_name": ""},
    {"rms_threshold": float("nan")}, {"energy_threshold": -1}, {"file_type": "txt"},
])
def test_invalid_conformer_inputs_return_errors_without_files(tmp_path, overrides):
    status, files, table = generate(tmp_path, **overrides)
    assert "color:red" in status
    assert files == []
    assert table.empty


def test_conformer_generation_requires_directory_before_embedding(monkeypatch):
    monkeypatch.setattr(conformer_generation, "generate_unique_conformers",
                        lambda *args: pytest.fail("Should require a directory before embedding"))
    status, files, table = generate(None)
    assert "color:red" in status
    assert files == []
    assert table.empty


def test_failed_embedding_reports_error(tmp_path, monkeypatch):
    monkeypatch.setattr(conformer_generation, "generate_unique_conformers",
                        lambda mol, *args: (mol, [], 0))
    status, files, table = generate(tmp_path)
    assert "color:red" in status and "No conformers" in status
    assert files == []
    assert table.empty


def test_unconverged_conformers_are_not_written(tmp_path, monkeypatch):
    monkeypatch.setattr(conformer_generation.AllChem, "MMFFOptimizeMoleculeConfs",
                        lambda mol, **kwargs: [(1, -1.0)] * mol.GetNumConformers())
    status, files, table = generate(tmp_path)
    assert "color:red" in status and "No conformers" in status
    assert files == []
    assert table.empty


def test_empty_editor_is_safe():
    from orca_webui.molecular_editor import smiles_from_editor

    assert smiles_from_editor({"smiles": ""}) == ""
