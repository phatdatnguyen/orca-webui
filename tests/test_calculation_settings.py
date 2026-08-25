"""Tests for the Calculation tab's settings: ORCA path, cores/memory, and menus.

These cover behaviour that is easy to break silently — the ORCA binary must be invoked
by absolute path or parallel runs refuse to start, ``%maxcore`` is per process so the
memory ceiling has to shrink as cores grow, and the functional/basis menus are rebuilt
per calculation type from module-level lists that used to be duplicated inline.
"""
import os
import sys

import pytest

gradio = pytest.importorskip("gradio")
from orca_webui import calculation


# --- resolve_orca_executable ------------------------------------------------

def test_resolve_orca_executable_returns_absolute_path_for_relative_entry(tmp_path, monkeypatch):
    """ORCA aborts parallel runs unless called by full pathname, so nothing relative
    may reach subprocess."""
    exe = tmp_path / "orca"
    exe.write_text("#!/bin/sh\n")
    monkeypatch.chdir(tmp_path)  # monkeypatch restores the cwd for the next test

    resolved = calculation.resolve_orca_executable("./orca")

    assert os.path.isabs(resolved)
    assert os.path.samefile(resolved, exe)


def test_resolve_orca_executable_passes_absolute_path_through():
    absolute = os.path.abspath(os.path.join(os.sep, "opt", "orca", "orca"))
    assert calculation.resolve_orca_executable(absolute) == absolute


def test_resolve_orca_executable_blank_falls_back_to_path_lookup(monkeypatch):
    monkeypatch.setattr(calculation, "default_orca_path", "/opt/orca/orca")
    assert calculation.resolve_orca_executable("") == "/opt/orca/orca"
    assert calculation.resolve_orca_executable(None) == "/opt/orca/orca"


def test_resolve_orca_executable_empty_when_nothing_found(monkeypatch):
    monkeypatch.setattr(calculation, "default_orca_path", "")
    assert calculation.resolve_orca_executable("") == ""
    assert calculation.resolve_orca_executable(None) == ""


# --- per-core memory ceiling ------------------------------------------------

def test_max_memory_per_core_keeps_total_within_physical_ram():
    """%maxcore is per process: nprocs x maxcore must not exceed physical RAM."""
    for n_cores in range(1, calculation.max_cores + 1):
        per_core = calculation.max_memory_per_core(n_cores)
        assert per_core >= 1
        assert n_cores * per_core <= calculation.max_memory


def test_max_memory_per_core_never_returns_zero():
    """A machine with less RAM than cores must still offer a usable slider."""
    assert calculation.max_memory_per_core(10 ** 6) == 1
    assert calculation.max_memory_per_core(0) >= 1  # guard against division by zero


def test_on_n_cores_change_clamps_value_down_but_not_up():
    lowered = calculation.on_n_cores_change(calculation.max_cores, calculation.max_memory)
    ceiling = calculation.max_memory_per_core(calculation.max_cores)
    assert lowered["maximum"] == ceiling
    assert lowered["value"] == ceiling  # clamped down from the full-RAM value

    kept = calculation.on_n_cores_change(1, 1)
    assert kept["value"] == 1  # already under the new ceiling, so left alone


def test_core_slider_maximum_is_half_the_logical_cpus():
    """ORCA fails to start when asked for more than half the logical CPU count."""
    import multiprocessing
    assert calculation.max_cores == max(1, multiprocessing.cpu_count() // 2)


# --- functional / basis-set menus per calculation type ----------------------

CALCULATION_TYPES = [
    "Single-Point",
    "Geometry Optimization",
    "Frequency",
    "Time-Dependent Density Functional Theory",
    "NMR Spectrum",
]


def _basis_update(calculation_type, current_basis="3-21G"):
    """The basis-set gr.update() is the 4th value returned by the handler."""
    return calculation.on_change_calculation_type(
        calculation_type, current_basis, calculation.GROUND_STATE)[3]


@pytest.mark.parametrize("calculation_type", CALCULATION_TYPES)
def test_pcsseg_offered_only_for_nmr(calculation_type):
    """Jensen's segmented sets are NMR-specific and must not clutter other types."""
    choices = _basis_update(calculation_type)["choices"]
    has_pcsseg = any(b.startswith("pcSseg") for b in choices)
    assert has_pcsseg == (calculation_type == "NMR Spectrum")


@pytest.mark.parametrize("functional", ["wB97X-D3BJ", "wB97X-D4"])
def test_new_functionals_available_everywhere(functional):
    assert functional in calculation.FUNCTIONALS
    assert functional in calculation.TDDFT_FUNCTIONALS


def test_basis_set_selection_survives_switch_when_still_offered():
    assert _basis_update("NMR Spectrum", "cc-pVTZ")["value"] == "cc-pVTZ"


def test_nmr_only_basis_is_dropped_when_leaving_nmr():
    """The dropdown allows custom values, so a stale pcSseg choice would otherwise be
    written into a non-NMR input file."""
    assert _basis_update("NMR Spectrum", "pcSseg-2")["value"] == "pcSseg-2"
    assert _basis_update("Single-Point", "pcSseg-2")["value"] == calculation.DEFAULT_BASIS_SET
