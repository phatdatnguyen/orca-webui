"""Tests for ORCA input-file generation (the write_*_orca_input family and helpers)."""
import io

import pytest

from orca_webui import utils


# --- get_orca_method_keyword -------------------------------------------------

@pytest.mark.parametrize("method_type, method_name, functional, basis, expected", [
    ("HF", "", "B3LYP", "6-31G", "HF 6-31G"),
    ("DFT", "", "PBE0", "def2-SVP", "PBE0 def2-SVP"),
    ("MP2", "", "B3LYP", "cc-pVDZ", "MP2 cc-pVDZ"),
    ("CCSD", "", "B3LYP", "cc-pVDZ", "CCSD cc-pVDZ"),
    ("CCSD(T)", "", "B3LYP", "cc-pVDZ", "CCSD(T) cc-pVDZ"),
    ("Semi-empirical", "AM1", "B3LYP", "", "AM1 "),
])
def test_method_keyword(method_type, method_name, functional, basis, expected):
    assert utils.get_orca_method_keyword(method_type, method_name, functional, basis) == expected


# --- get_orca_header ---------------------------------------------------------

def test_header_writes_pal_and_maxcore():
    buf = io.StringIO()
    utils.get_orca_header(buf, n_proc=4, memory=2, solvation=False,
                          solvation_model=None, solvent=None)
    text = buf.getvalue()
    assert "%pal" in text and "nprocs 4" in text
    # memory is GB -> maxcore in MB
    assert "%maxcore 2048" in text


def test_header_cpcm_solvation():
    buf = io.StringIO()
    utils.get_orca_header(buf, n_proc=1, memory=1, solvation=True,
                          solvation_model="cpcm", solvent="water")
    text = buf.getvalue()
    assert "%cpcm" in text
    assert 'solvent "water"' in text
    assert "smd true" not in text


def test_header_smd_solvation():
    buf = io.StringIO()
    utils.get_orca_header(buf, n_proc=1, memory=1, solvation=True,
                          solvation_model="smd", solvent="dmso")
    text = buf.getvalue()
    assert "smd true" in text
    assert 'smdsolvent "dmso"' in text


def test_header_tddft_block():
    buf = io.StringIO()
    utils.get_orca_header(buf, n_proc=1, memory=1, solvation=False,
                          solvation_model=None, solvent=None,
                          tddft=True, n_states=15, maxdim=7, root=2)
    text = buf.getvalue()
    assert "%tddft" in text
    assert "nroots 15" in text
    assert "maxdim 7" in text
    assert "iroot 2" in text


def test_header_omits_iroot_for_a_single_point():
    """iroot selects which state's gradient to follow, so it is Opt-only.

    Writing it for a single-point TD-DFT would be noise; ORCA needs it only when a
    geometry is being optimized.
    """
    buf = io.StringIO()
    utils.get_orca_header(buf, n_proc=1, memory=1, solvation=False,
                          solvation_model=None, solvent=None,
                          tddft=True, n_states=5, maxdim=5)
    text = buf.getvalue()
    assert "nroots 5" in text
    assert "iroot" not in text


def test_header_omits_maxcore_when_memory_none():
    buf = io.StringIO()
    utils.get_orca_header(buf, n_proc=1, memory=None, solvation=False,
                          solvation_model=None, solvent=None)
    assert "%maxcore" not in buf.getvalue()


# --- get_orca_coordinates ----------------------------------------------------

def test_coordinates_block(water_mol):
    buf = io.StringIO()
    utils.get_orca_coordinates(buf, water_mol, charge=0, multiplicity=1)
    text = buf.getvalue()
    assert text.startswith("* xyz 0 1")
    assert text.rstrip().endswith("*")
    # water = 1 O + 2 H: exactly three atom lines between header and closing "*"
    lines = text.strip().splitlines()
    assert lines[0] == "* xyz 0 1"
    assert lines[-1] == "*"
    atom_lines = lines[1:-1]
    assert len(atom_lines) == 3
    assert atom_lines[0].startswith("O ")
    assert sum(l.startswith("H ") for l in atom_lines) == 2


# --- write_*_orca_input ------------------------------------------------------

def _read(path):
    with open(path) as f:
        return f.read()


def test_write_sp_input(tmp_path, water_mol):
    base = str(tmp_path / "sp")
    utils.write_sp_orca_input(water_mol, base, "DFT", "", functional="B3LYP",
                              basis="3-21G", charge=0, multiplicity=1,
                              n_proc=1, memory=1)
    text = _read(base + ".inp")
    assert text.startswith("! B3LYP 3-21G SP")
    assert "* xyz 0 1" in text


def test_write_opt_input(tmp_path, water_mol):
    base = str(tmp_path / "opt")
    utils.write_opt_orca_input(water_mol, base, "HF", "", basis="6-31G",
                               n_proc=1, memory=1)
    assert _read(base + ".inp").startswith("! HF 6-31G Opt")


def test_write_opt_freq_input(tmp_path, water_mol):
    base = str(tmp_path / "freq")
    utils.write_opt_freq_orca_input(water_mol, base, "DFT", "", functional="PBE",
                                    basis="def2-SVP", n_proc=1, memory=1)
    assert _read(base + ".inp").startswith("! PBE def2-SVP Opt Freq")


def test_write_tddft_input_is_a_single_point(tmp_path, water_mol):
    """TD-DFT must not optimize: the spectrum belongs to the geometry it was given.

    The old writer appended Opt, which (with iroot 0) ORCA turned into a no-op
    optimization that still cost a full extra pass, and which made ORCA reject
    B88-containing functionals such as B3LYP outright.
    """
    base = str(tmp_path / "td")
    utils.write_tddft_orca_input(water_mol, base, "DFT", "", functional="PBE0",
                                 basis="def2-SVP", n_states=8, maxdim=4,
                                 n_proc=1, memory=1)
    text = _read(base + ".inp")
    assert text.startswith("! PBE0 def2-SVP RIJCOSX TightSCF")
    # The method line must carry no Opt (and no Freq).
    assert "Opt" not in text.splitlines()[0]
    assert "nroots 8" in text and "maxdim 4" in text
    assert "iroot" not in text



def test_write_opt_input_ground_state_has_no_tddft_block(tmp_path, water_mol):
    base = str(tmp_path / "opt")
    utils.write_opt_orca_input(water_mol, base, "DFT", "", functional="PBE",
                               basis="def2-SVP", n_proc=1, memory=1)
    text = _read(base + ".inp")
    assert text.startswith("! PBE def2-SVP Opt")
    assert "%tddft" not in text


def test_write_opt_input_excited_state_follows_the_root(tmp_path, water_mol):
    """Excited-state optimization is what gives the geometry emission is measured at."""
    base = str(tmp_path / "opt_es")
    utils.write_opt_orca_input(water_mol, base, "DFT", "", functional="PBE",
                               basis="def2-SVP", n_proc=1, memory=1,
                               excited_state=True, n_states=5, maxdim=5, root=1)
    text = _read(base + ".inp")
    assert text.startswith("! PBE def2-SVP Opt")
    assert "%tddft" in text
    assert "nroots 5" in text
    # iroot 0 would optimize the ground state, making the excited-state option a no-op.
    assert "iroot 1" in text


def test_write_nmr_input_without_coupling(tmp_path, water_mol):
    base = str(tmp_path / "nmr")
    utils.write_nmr_orca_input(water_mol, base, "DFT", functional="B3LYP",
                               basis="6-31G", spin_spin_coupling=False,
                               n_proc=1, memory=1)
    text = _read(base + ".inp")
    assert text.startswith("! B3LYP 6-31G NMR")
    assert "%eprnmr" not in text


def test_write_nmr_input_with_coupling(tmp_path, water_mol):
    base = str(tmp_path / "nmr_ss")
    utils.write_nmr_orca_input(water_mol, base, "DFT", functional="B3LYP",
                               basis="6-31G", spin_spin_coupling=True,
                               n_proc=1, memory=1)
    text = _read(base + ".inp")
    assert "%eprnmr" in text
    assert "ssall" in text
    assert "SpinSpinRThresh" in text


def test_write_input_applies_solvation(tmp_path, water_mol):
    base = str(tmp_path / "sp_solv")
    utils.write_sp_orca_input(water_mol, base, "DFT", "", functional="B3LYP",
                              basis="3-21G", solvation=True,
                              solvation_model="cpcm", solvent="ethanol",
                              n_proc=1, memory=1)
    assert 'solvent "ethanol"' in _read(base + ".inp")
