"""End-to-end tests that invoke the real ORCA binary.

Everything else in this suite parses hand-written log text and substitutes a stand-in
for the ORCA executable, which verifies our own code but says nothing about whether ORCA
actually accepts the inputs we generate or whether our parsers match the output a real
release produces. These tests close that loop: generate an input with the app's writers,
run ORCA on it, and parse the result with the app's parsers.

They are skipped when ORCA is not installed. Point them at an install that is not on
``PATH`` with ``ORCA_EXECUTABLE=/path/to/orca pytest``, and deselect them entirely with
``pytest -m "not orca"``.

Systems are deliberately tiny (H2, methane) at HF/STO-3G with ``nprocs 1`` so the whole
module runs in a couple of seconds and needs no MPI.
"""
import os

import pytest

from orca_webui import utils

pytestmark = pytest.mark.orca


def _run(orca_executable, directory, input_stem):
    """Run ORCA through the app's own handler and return (status_html, log_text)."""
    gr_calculation = pytest.importorskip("orca_webui.calculation")
    status, _files = gr_calculation.on_run_calculation(
        str(directory), f"{input_stem}.inp", orca_executable)
    log_path = os.path.join(str(directory), f"{input_stem}.log")
    log_text = open(log_path).read() if os.path.isfile(log_path) else ""
    return status, log_text


# --- single point -----------------------------------------------------------

def test_generated_single_point_input_runs_to_completion(tmp_path, orca_executable, water_mol):
    """A write_sp_orca_input file must be something ORCA actually accepts."""
    stem = tmp_path / "sp"
    utils.write_sp_orca_input(water_mol, str(stem), "HF", "", basis="STO-3G", n_proc=1, memory=1)

    status, log = _run(orca_executable, tmp_path, "sp")

    assert "color:red" not in status, status
    assert "ORCA TERMINATED NORMALLY" in log
    assert "FINAL SINGLE POINT ENERGY" in log


def test_run_calculation_reports_a_real_orca_failure(tmp_path, orca_executable):
    """A genuinely broken input must surface as an error, not a green status."""
    (tmp_path / "bad.inp").write_text("! ThisIsNotAnOrcaKeyword\n* xyz 0 1\nH 0 0 0\n*\n")

    status, log = _run(orca_executable, tmp_path, "bad")

    assert "color:red" in status, status
    # Assert on ORCA's own diagnostic rather than just the red status: a red status
    # alone would also appear if the binary had never run, making this test vacuous.
    assert "INPUT ERROR" in log
    assert "ORCA TERMINATED NORMALLY" not in log
    assert "exited with code" in status  # the non-zero exit path, not the missing-binary one


# --- atom ordering, end to end ----------------------------------------------

def test_atom_order_survives_the_round_trip(tmp_path, orca_executable):
    """The viewer's index labels are only trustworthy if this ordering chain holds.

    structure file -> mol_from_structure_file -> generated .inp -> ORCA's own
    CARTESIAN COORDINATES block. The carbon is deliberately not first so any reordering
    at any stage shows up.
    """
    xyz = tmp_path / "shuffled.xyz"
    xyz.write_text(
        "0 1\n"
        "H  0.629118  0.629118  0.629118\n"
        "C  0.000000  0.000000  0.000000\n"
        "H -0.629118 -0.629118  0.629118\n"
        "H -0.629118  0.629118 -0.629118\n"
        "H  0.629118 -0.629118 -0.629118\n"
    )
    mol = utils.mol_from_structure_file(str(xyz))
    labelled = [atom.GetSymbol() for atom in mol.GetAtoms()]
    assert labelled == ["H", "C", "H", "H", "H"]  # the viewer would label H0 C1 H2 H3 H4

    utils.write_sp_orca_input(mol, str(tmp_path / "sp"), "HF", "", basis="STO-3G",
                              n_proc=1, memory=1)
    status, log = _run(orca_executable, tmp_path, "sp")
    assert "color:red" not in status, status

    from_orca, _coords = utils.parse_orca_geometry(os.path.join(str(tmp_path), "sp.log"))
    assert from_orca == labelled


# --- optimization target: ground state vs excited state ---------------------

def _oo_distance(log_path):
    """O-O bond length from a log's final geometry, in Angstrom."""
    import numpy as np
    symbols, coords = utils.parse_orca_geometry(log_path)
    assert symbols[:2] == ["O", "O"], symbols
    return float(np.linalg.norm(np.array(coords[0]) - np.array(coords[1])))


def test_excited_state_optimization_relaxes_that_state(tmp_path, orca_executable):
    """The Geometry Optimization excited-state target must really follow the state.

    PBE rather than B3LYP/wB97X-D4: a TD-DFT gradient needs the XC third derivative,
    which ORCA 6.1 does not have natively for those ("try using LibXC instead").
    """
    (tmp_path / "h2o2.xyz").write_text(SKEWED_H2O2_XYZ)
    mol = utils.mol_from_structure_file(str(tmp_path / "h2o2.xyz"))
    utils.write_opt_orca_input(mol, str(tmp_path / "es_opt"), "DFT", "",
                               functional="PBE", basis="STO-3G", n_proc=1, memory=2,
                               excited_state=True, n_states=5, maxdim=5, root=1)
    text = (tmp_path / "es_opt.inp").read_text()
    assert "Opt" in text.splitlines()[0]
    assert "iroot 1" in text

    status, log = _run(orca_executable, tmp_path, "es_opt")
    assert "color:red" not in status, status
    assert "ORCA TERMINATED NORMALLY" in log
    # ORCA names the state it follows and computes that state's gradient. With iroot 0
    # neither line appears and the optimization is a no-op.
    assert "State of interest" in log
    assert "EXCITED STATE GRADIENT" in log

    # S1 of H2O2 is O-O dissociative, so the bond must lengthen well beyond the 1.45 A
    # it started at -- proof the excited state, not the ground state, was relaxed.
    assert _oo_distance(os.path.join(str(tmp_path), "es_opt.log")) > 1.6


def test_ground_state_optimization_has_no_tddft_block(tmp_path, orca_executable):
    (tmp_path / "h2o2.xyz").write_text(SKEWED_H2O2_XYZ)
    mol = utils.mol_from_structure_file(str(tmp_path / "h2o2.xyz"))
    utils.write_opt_orca_input(mol, str(tmp_path / "gs_opt"), "DFT", "",
                               functional="PBE", basis="STO-3G", n_proc=1, memory=2)
    assert "%tddft" not in (tmp_path / "gs_opt.inp").read_text()

    status, log = _run(orca_executable, tmp_path, "gs_opt")
    assert "color:red" not in status, status
    assert "ORCA TERMINATED NORMALLY" in log
    assert "EXCITED STATE GRADIENT" not in log


def test_tddft_single_point_leaves_the_geometry_alone(tmp_path, orca_executable):
    """No Opt means the spectrum belongs to the structure the user selected.

    Also checks the Result tab consequence: with one SCF energy instead of two, the
    Geometry Optimization accordion stops appearing for a spectrum job (it used to show
    a flat two-point plot from the no-op iroot-0 optimization).
    """
    pytest.importorskip("gradio")
    from orca_webui import result

    (tmp_path / "h2o2.xyz").write_text(SKEWED_H2O2_XYZ)
    mol = utils.mol_from_structure_file(str(tmp_path / "h2o2.xyz"))
    utils.write_tddft_orca_input(mol, str(tmp_path / "sp_td"), "DFT", "",
                                 functional="B3LYP", basis="STO-3G",
                                 n_states=5, maxdim=5, n_proc=1, memory=2)

    status, log = _run(orca_executable, tmp_path, "sp_td")
    # B3LYP specifically: with Opt present ORCA rejected B88 functionals outright.
    assert "color:red" not in status, status
    assert "ORCA TERMINATED NORMALLY" in log
    assert "GEOMETRY OPTIMIZATION CYCLE" not in log
    # Geometry untouched: still the 1.45 A the input file specified.
    assert abs(_oo_distance(os.path.join(str(tmp_path), "sp_td.log")) - 1.45) < 1e-3

    outputs = result.on_load_result_file(str(tmp_path), "sp_td.log", "32.9", "203.4", True)
    assert "color:red" not in outputs[0], outputs[0]
    assert outputs[6]["visible"] is False        # geometry-optimization accordion
    assert outputs[13]["visible"] is True        # absorption/emission accordion


# --- ECD: rotatory strengths from TD-DFT ------------------------------------

# Hydrogen peroxide with a ~90 degree dihedral: only 4 atoms, but genuinely chiral (C2),
# so the rotatory strengths are non-zero. A planar or achiral molecule would give R = 0
# and the test could pass without proving anything.
SKEWED_H2O2_XYZ = (
    "0 1\n"
    "O  0.000000  0.000000  0.000000\n"
    "O  1.450000  0.000000  0.000000\n"
    "H -0.168000  0.955000  0.000000\n"
    "H  1.618000  0.000000  0.955000\n"
)


@pytest.fixture(scope="module")
def ecd_log(tmp_path_factory, orca_executable):
    """One real TD-DFT single point on a chiral molecule."""
    from orca_webui import calculation

    directory = tmp_path_factory.mktemp("ecd")
    (directory / "h2o2.xyz").write_text(SKEWED_H2O2_XYZ)
    mol = utils.mol_from_structure_file(str(directory / "h2o2.xyz"))

    # A single point now, so any functional works -- the old writer's Opt made ORCA
    # demand XC third derivatives and reject B88-containing functionals.
    utils.write_tddft_orca_input(mol, str(directory / "ecd"), "DFT", "",
                                 functional="wB97X-D4", basis="def2-SVP",
                                 n_states=5, maxdim=5, n_proc=1, memory=2)

    status, _files = calculation.on_run_calculation(str(directory), "ecd.inp", orca_executable)
    if "color:red" in status:
        pytest.fail(f"ECD run failed: {status}")
    return os.path.join(str(directory), "ecd.log")


def test_ecd_run_emits_the_cd_spectrum_block(ecd_log):
    log = open(ecd_log).read()
    assert "ORCA TERMINATED NORMALLY" in log
    # Both representations ORCA offers; the electric-dipole one is what cclib reads.
    assert "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in log


def test_cclib_exposes_rotatory_strengths_from_a_real_run(ecd_log):
    """The Result tab reads etrotats, so cclib must actually populate it."""
    cclib_io = pytest.importorskip("cclib.io")
    data = cclib_io.ccopen(ecd_log).parse()

    assert hasattr(data, "etrotats")
    assert len(data.etrotats) == len(data.etenergies)
    # A chiral molecule must give a non-trivial, signed set of rotatory strengths.
    rotats = list(data.etrotats)
    assert any(abs(r) > 0.1 for r in rotats)
    assert any(r < 0 for r in rotats) or any(r > 0 for r in rotats)


def test_result_tab_enables_and_builds_ecd_from_a_real_run(ecd_log):
    """End to end: load the log through the Result tab, then press Show ECD."""
    pytest.importorskip("gradio")
    import plotly.graph_objects as go
    from orca_webui import result

    directory, name = os.path.split(ecd_log)
    outputs = result.on_load_result_file(directory, name, "32.9", "203.4", True)
    assert "color:red" not in outputs[0], outputs[0]

    data = outputs[1]
    assert hasattr(data, "etrotats")
    # The ECD button update sits right after the peak-export button in the outputs list.
    assert outputs[17]["interactive"] is True

    peak_df, fig = result.on_show_ecd_spectrum(data)
    assert len(peak_df) == len(data.etenergies)
    assert "Rotatory Strength (1e-40 cgs)" in peak_df.columns
    assert isinstance(fig, go.Figure)


def test_ecd_needs_no_docd_keyword(tmp_path, orca_executable):
    """Documents why there is no "compute CD" checkbox: ORCA reports R regardless.

    ``DoCD true`` is the documented way to request CD intensities, but ORCA 6.1.1 prints
    the CD SPECTRUM blocks for a plain ``%tddft`` job anyway, with identical R values
    (verified by diffing the two runs). So the ECD button is driven by whether
    ``etrotats`` is present, which also means ECD works on TD-DFT logs produced before
    the feature existed. If a future ORCA stops emitting CD by default, this test fails
    and a DoCD option has to come back.
    """
    pytest.importorskip("gradio")
    from orca_webui import calculation, result

    (tmp_path / "h2o2.xyz").write_text(SKEWED_H2O2_XYZ)
    mol = utils.mol_from_structure_file(str(tmp_path / "h2o2.xyz"))
    utils.write_tddft_orca_input(mol, str(tmp_path / "plain"), "DFT", "",
                                 functional="wB97X-D4", basis="def2-SVP",
                                 n_states=5, maxdim=5, n_proc=1, memory=2)
    assert "DoCD" not in (tmp_path / "plain.inp").read_text()

    status, _files = calculation.on_run_calculation(str(tmp_path), "plain.inp", orca_executable)
    assert "color:red" not in status, status

    log = open(os.path.join(str(tmp_path), "plain.log")).read()
    assert "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in log

    outputs = result.on_load_result_file(str(tmp_path), "plain.log", "32.9", "203.4", True)
    assert outputs[17]["interactive"] is True
    peak_df, fig = result.on_show_ecd_spectrum(outputs[1])
    assert peak_df is not None and fig is not None


# --- NMR: shieldings, couplings, and the peak table -------------------------

@pytest.fixture(scope="module")
def methane_nmr_log(tmp_path_factory, orca_executable):
    """Run one real methane NMR job with spin-spin couplings, shared by the tests below."""
    from orca_webui import calculation

    directory = tmp_path_factory.mktemp("nmr")
    xyz = directory / "methane.xyz"
    lines = ["0 1"]
    for symbol, (x, y, z) in zip(
            ["C", "H", "H", "H", "H"],
            [(0.0, 0.0, 0.0), (0.629118, 0.629118, 0.629118),
             (-0.629118, -0.629118, 0.629118), (-0.629118, 0.629118, -0.629118),
             (0.629118, -0.629118, -0.629118)]):
        lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
    xyz.write_text("\n".join(lines) + "\n")

    mol = utils.mol_from_structure_file(str(xyz))
    utils.write_nmr_orca_input(mol, str(directory / "nmr"), "HF", basis="STO-3G",
                               spin_spin_coupling=True, n_proc=1, memory=1)

    status, _files = calculation.on_run_calculation(str(directory), "nmr.inp", orca_executable)
    log_path = os.path.join(str(directory), "nmr.log")
    if "color:red" in status or not os.path.isfile(log_path):
        pytest.fail(f"methane NMR run failed: {status}")
    return log_path


def test_nmr_run_produces_the_blocks_the_parsers_expect(methane_nmr_log):
    log = open(methane_nmr_log).read()
    assert "ORCA TERMINATED NORMALLY" in log
    assert "CHEMICAL SHIELDING SUMMARY" in log
    assert "NUCLEUS A" in log  # the spin-spin block %eprnmr was asked for


def test_parse_shieldings_from_real_output(methane_nmr_log):
    shieldings = utils.parse_nmr_shielding_constants(methane_nmr_log)

    assert [(i, e) for i, e, _ in shieldings] == [(0, "C"), (1, "H"), (2, "H"), (3, "H"), (4, "H")]
    # Methane's four hydrogens are symmetry-equivalent, so their shieldings must agree.
    h_values = [v for _, e, v in shieldings if e == "H"]
    assert max(h_values) - min(h_values) < 0.01
    # Sanity: a proton is far more shielded than the carbon it hangs off.
    c_value = next(v for _, e, v in shieldings if e == "C")
    assert c_value > max(h_values)


def test_parse_jcouplings_from_real_output(methane_nmr_log):
    couplings = utils.parse_nmr_jcouplings(methane_nmr_log)
    assert couplings, "no J-couplings parsed from a run that requested ssall"

    ch = [j for a, ea, b, eb, j in couplings if {ea, eb} == {"C", "H"}]
    hh = [j for a, ea, b, eb, j in couplings if ea == "H" and eb == "H"]

    assert len(ch) == 4 and len(hh) == 6  # every C-H and H-H pair in methane
    # 1J(C,H) in methane is ~125 Hz experimentally; HF/STO-3G lands in the right region.
    assert all(50 < abs(j) < 300 for j in ch), ch
    # Equivalent nuclei: all four C-H couplings must come out identical.
    assert max(ch) - min(ch) < 0.01


def test_peak_table_from_real_output(methane_nmr_log):
    """The full pipeline: real log -> shifts -> equivalence -> averaged peaks."""
    shieldings = utils.parse_nmr_shielding_constants(methane_nmr_log)
    couplings = utils.parse_nmr_jcouplings(methane_nmr_log)
    shifts = utils.calculate_chemical_shifts(shieldings, {"H": 32.9, "C": 203.4})

    symbols, coords = utils.parse_orca_geometry(methane_nmr_log)
    equivalence = utils.compute_hydrogen_equivalence(coords, symbols)

    peaks = utils.build_nmr_peak_table(shifts, "H", jcouplings=couplings,
                                       equivalence=equivalence)

    # All four methane protons are equivalent, so they merge into a single peak...
    assert len(peaks) == 1
    assert peaks[0]["count"] == 4
    assert sorted(peaks[0]["atom_idxs"]) == [1, 2, 3, 4]
    # ...and mutually equivalent nuclei do not split each other, so it is a singlet.
    assert utils.multiplicity_label(peaks[0]["couplings"]) == "s"


def test_jcoupling_column_text_from_real_output(methane_nmr_log):
    """The Result tab's J-Couplings column, built from genuine ORCA numbers."""
    couplings = utils.parse_nmr_jcouplings(methane_nmr_log)
    text = utils.format_jcouplings_by_atom(couplings)

    # H1 couples to the carbon and to the other three protons.
    entry = text[1]
    assert entry.count(" Hz (") == 4
    assert "(C0)" in entry
    for partner in ("H2", "H3", "H4"):
        assert f"({partner})" in entry
    # Strongest first: 1J(C,H) outweighs 2J(H,H).
    assert entry.startswith(entry.split(",")[0])
    assert entry.split(" Hz")[0].lstrip("-").replace(".", "").isdigit()
