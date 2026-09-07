"""Run-handler regressions using tiny Python programs in place of ORCA."""
import sys
import threading
from concurrent.futures import ThreadPoolExecutor

import pytest

pytest.importorskip("gradio")
from orca_webui import calculation


def test_run_uses_selected_directory_for_relative_input_and_output(tmp_path, monkeypatch):
    directory = tmp_path / "calculation work"
    directory.mkdir()
    (directory / "geometry.xyz").write_text("geometry from selected directory")
    (directory / "run.inp").write_text(
        "from pathlib import Path\n"
        "print(Path('geometry.xyz').read_text())\n"
        "Path('scratch.dat').write_text('calculation artifact')\n"
    )
    monkeypatch.chdir(tmp_path)

    status, files = calculation.on_run_calculation("calculation work", "run.inp", sys.executable)

    assert "color:green" in status
    assert "geometry from selected directory" in (directory / "run.log").read_text()
    assert "scratch.dat" in files
    assert (directory / "scratch.dat").is_file()
    assert not (tmp_path / "scratch.dat").exists()


def test_missing_input_preserves_both_previous_logs(tmp_path):
    (tmp_path / "run.log").write_text("previous run")
    (tmp_path / "run.log.bak").write_text("older run")

    status, _ = calculation.on_run_calculation(str(tmp_path), "run.inp", sys.executable)

    assert "color:red" in status
    assert (tmp_path / "run.log").read_text() == "previous run"
    assert (tmp_path / "run.log.bak").read_text() == "older run"


def test_failed_process_launch_restores_previous_log(tmp_path, monkeypatch):
    (tmp_path / "run.inp").write_text("input")
    (tmp_path / "run.log").write_text("previous run")

    def fail_launch(*args, **kwargs):
        raise OSError("invalid executable format")

    monkeypatch.setattr(calculation.subprocess, "Popen", fail_launch)
    status, _ = calculation.on_run_calculation(str(tmp_path), "run.inp", sys.executable)

    assert "invalid executable format" in status
    assert (tmp_path / "run.log").read_text() == "previous run"
    assert calculation._orca_process is None


def test_concurrent_run_cannot_replace_active_process_or_log(tmp_path, monkeypatch):
    (tmp_path / "run.inp").write_text("input")
    launched = threading.Event()
    finish = threading.Event()

    class Process:
        def __init__(self, *args, stdout, **kwargs):
            stdout.write("active calculation output")

        def wait(self):
            launched.set()
            assert finish.wait(5), "test did not release the running calculation"
            return 0

    monkeypatch.setattr(calculation.subprocess, "Popen", Process)
    with ThreadPoolExecutor(max_workers=1) as executor:
        first = executor.submit(
            calculation.on_run_calculation, str(tmp_path), "run.inp", sys.executable)
        try:
            assert launched.wait(5)
            running_process = calculation._orca_process
            status, _ = calculation.on_run_calculation(str(tmp_path), "run.inp", sys.executable)
            assert "already running" in status
            assert calculation._orca_process is running_process
            assert (tmp_path / "run.log").read_text() == "active calculation output"
            assert not (tmp_path / "run.log.bak").exists()
        finally:
            finish.set()
        assert "color:green" in first.result(timeout=5)[0]
    assert calculation._orca_process is None


def test_stop_reports_process_permission_error(monkeypatch):
    class Process:
        pid = 123

        def poll(self):
            return None

    def inaccessible_process(pid):
        raise calculation.psutil.AccessDenied(pid)

    monkeypatch.setattr(calculation, "_orca_process", Process())
    monkeypatch.setattr(calculation, "_orca_stopped_by_user", False)
    monkeypatch.setattr(calculation.psutil, "Process", inaccessible_process)

    status = calculation.on_stop_calculation()

    assert "Unable to stop calculation" in status
    assert calculation._orca_stopped_by_user is False


@pytest.mark.parametrize("name", ["../outside.inp", "run.log"])
def test_run_rejects_invalid_input_names_without_touching_logs(tmp_path, name):
    directory = tmp_path / "work"
    directory.mkdir()
    (tmp_path / "outside.inp").write_text("print('should not run')")
    (directory / "run.log").write_text("existing result")

    status, _ = calculation.on_run_calculation(str(directory), name, sys.executable)

    assert "color:red" in status
    assert not (tmp_path / "outside.log").exists()
    assert (directory / "run.log").read_text() == "existing result"


def _generate(tmp_path, name, *, use_mm=False):
    return calculation.on_generate_input_file(
        str(tmp_path), "water.xyz", "Single-Point", calculation.GROUND_STATE,
        use_mm, "MMFF", 200, False, "cpcm", "water", "DFT", "", "B3LYP",
        "3-21G", 10, 5, 1, False, 0, 1, 1, 1, name)


def test_input_name_with_extension_generates_and_selects_one_inp(tmp_path, monkeypatch, water_mol):
    monkeypatch.setattr(calculation, "mol_from_structure_file", lambda path: water_mol)

    status, files = _generate(tmp_path, "water.inp")

    assert "color:green" in status
    assert files == ["water.inp"]
    _, dropdown = calculation.on_working_directory_file_list_change(files, "water.inp")
    assert dropdown["value"] == "water.inp"


@pytest.mark.parametrize("name", ["", ".inp", "../escape"])
def test_generate_rejects_empty_or_external_output_name(tmp_path, monkeypatch, water_mol, name):
    monkeypatch.setattr(calculation, "mol_from_structure_file", lambda path: water_mol)

    status, files = _generate(tmp_path, name)

    assert "color:red" in status
    assert files == []
    assert not (tmp_path.parent / "escape.inp").exists()


def test_mm_parameter_failure_does_not_generate_misleading_input(tmp_path, monkeypatch, water_mol):
    monkeypatch.setattr(calculation, "mol_from_structure_file", lambda path: water_mol)
    monkeypatch.setattr(calculation.AllChem, "MMFFOptimizeMolecule", lambda *args, **kwargs: -1)

    status, files = _generate(tmp_path, "water", use_mm=True)

    assert "parameters are unavailable" in status
    assert files == []
