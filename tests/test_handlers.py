"""Regression tests for the Gradio event handlers.

These modules import Gradio (and, for the working directory, nglview) at import time,
so the whole module is skipped when those optional UI deps are absent. Run them in the
full application environment (the venv) to exercise the fixed bugs end-to-end.
"""
import os
import sys

import pytest


class FakeProgress:
    """Stand-in for gr.Progress().

    ``.tqdm()`` passes the iterable straight through, and the instance is callable
    because the conformer search reports fractional progress with ``progress(x, desc=)``
    as well as wrapping loops in ``.tqdm``.
    """
    def tqdm(self, iterable, total=None, desc=None):
        return iterable

    def __call__(self, fraction, desc=None):
        return None


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

    status, _files, _table = conformer_generation.on_generate_conformers(
        str(tmp_path), "not_a_valid_smiles", 0, 1, 3, 0.1, 0.5, "conformer", "xyz",
        progress=FakeProgress())
    assert "invalid SMILES" in status
    assert "color:red" in status
    assert not list(tmp_path.glob("conformer*"))


def test_generate_conformers_writes_expected_files(tmp_path):
    pytest.importorskip("gradio")
    pytest.importorskip("gradio_molecule2d")
    from orca_webui import conformer_generation

    status, _files, table = conformer_generation.on_generate_conformers(
        str(tmp_path), "CCO", 0, 1, 3, 0.1, 0.5, "conformer", "xyz",
        progress=FakeProgress())
    # Green when all 3 were kept, orange when deduplication left fewer; either way the
    # run succeeded, and only the error path is coloured red.
    assert "color:red" not in status
    produced = sorted(p.name for p in tmp_path.glob("conformer_*.xyz"))
    # RDKit may embed fewer than requested; every produced file is named
    # sequentially and no GetConformer(bad_id) crash occurs (the fixed bug).
    assert len(produced) >= 1
    assert produced[0] == "conformer_1.xyz"
    for name in produced:
        assert (tmp_path / name).read_text().strip() != ""
    # One table row per file written.
    assert len(table) == len(produced)


# --- no working directory open: shared states are None ----------------------
# os.listdir(None) silently lists the *current* directory rather than raising, so a
# None path used to yield a plausible-looking file list that blew up downstream.

def test_file_list_handlers_tolerate_no_working_directory():
    pytest.importorskip("gradio")
    pytest.importorskip("nglview")
    from orca_webui import calculation, result, utils, working_directory

    assert utils.get_files_in_working_directory(None) == []
    assert utils.get_files_in_working_directory("") == []

    empty = working_directory.on_file_list_change(None)
    assert len(empty) == 0
    assert list(empty.columns) == ["File", "Type", "Modified"]

    structures, inputs = calculation.on_working_directory_file_list_change(None, None)
    assert structures["choices"] == [] and inputs["choices"] == []

    assert result.on_working_directory_file_list_change(None)["choices"] == []


# --- ECD: the TD-DFT checkbox and the Result-tab button ---------------------

class FakeTDDFTData:
    """Stand-in for the cclib object stashed in the Result tab's data state."""
    def __init__(self, with_rotats):
        self.etenergies = [53985.7, 54496.6, 61380.2]   # cm^-1
        self.etoscs = [0.0007, 0.0031, 0.0001]
        self.etsyms = ["Singlet", "Singlet", "Singlet"]
        self.etsecs = [[((4,), (5,), 0.7)], [((4,), (6,), 0.6)], [((3,), (5,), 0.5)]]
        if with_rotats:
            self.etrotats = [-5.49159, 1.88894, 3.39576]


ALL_CALCULATION_TYPES = ("Single-Point", "Geometry Optimization", "Frequency",
                        "NMR Spectrum", "Time-Dependent Density Functional Theory")

# Indices into on_change_calculation_type's return tuple, matching the outputs list
# wired in calculation_tab_content.
OPT_TARGET, N_STATES, MAXDIM, ROOT = 4, 5, 6, 7


def test_optimize_target_offered_only_for_geometry_optimization():
    pytest.importorskip("gradio")
    from orca_webui import calculation

    for calculation_type in ALL_CALCULATION_TYPES:
        returned = calculation.on_change_calculation_type(
            calculation_type, "3-21G", calculation.GROUND_STATE)
        visible = returned[OPT_TARGET]["visible"]
        assert visible == (calculation_type == "Geometry Optimization"), \
            f"{calculation_type} -> {visible}"


def test_state_controls_shown_for_tddft_and_excited_state_optimization():
    """nroots/maxdim serve both; iroot is only meaningful when optimizing a state."""
    pytest.importorskip("gradio")
    from orca_webui import calculation

    tddft = calculation.on_change_calculation_type(
        "Time-Dependent Density Functional Theory", "3-21G", calculation.GROUND_STATE)
    assert tddft[N_STATES]["visible"] and tddft[MAXDIM]["visible"]
    assert not tddft[ROOT]["visible"]  # a single point has no gradient to follow

    ground_opt = calculation.on_change_calculation_type(
        "Geometry Optimization", "3-21G", calculation.GROUND_STATE)
    assert not ground_opt[N_STATES]["visible"] and not ground_opt[ROOT]["visible"]

    excited_opt = calculation.on_change_calculation_type(
        "Geometry Optimization", "3-21G", calculation.EXCITED_STATE)
    assert excited_opt[N_STATES]["visible"] and excited_opt[ROOT]["visible"]


def test_optimize_target_change_toggles_the_state_controls():
    pytest.importorskip("gradio")
    from orca_webui import calculation

    n_states, maxdim, root = calculation.on_optimize_target_change(
        "Geometry Optimization", calculation.EXCITED_STATE)
    assert n_states["visible"] and maxdim["visible"] and root["visible"]

    n_states, maxdim, root = calculation.on_optimize_target_change(
        "Geometry Optimization", calculation.GROUND_STATE)
    assert not n_states["visible"] and not root["visible"]

    # Switching the radio must not reveal anything for a non-optimization type.
    n_states, _maxdim, root = calculation.on_optimize_target_change(
        "Single-Point", calculation.EXCITED_STATE)
    assert not n_states["visible"] and not root["visible"]


def test_root_ceiling_tracks_the_number_of_excited_states():
    """iroot N needs at least N roots solved, so lowering n_states must lower iroot."""
    pytest.importorskip("gradio")
    from orca_webui import calculation

    lowered = calculation.on_n_states_change(5, 9)
    assert lowered["maximum"] == 5
    assert lowered["value"] == 5          # clamped down from an unreachable state

    kept = calculation.on_n_states_change(20, 3)
    assert kept["maximum"] == 20
    assert kept["value"] == 3             # already reachable, so left alone


def test_show_ecd_spectrum_needs_rotatory_strengths():
    """A result with no rotatory strengths must produce no ECD output."""
    pytest.importorskip("gradio")
    from orca_webui import result

    assert result.on_show_ecd_spectrum(FakeTDDFTData(with_rotats=False)) == (None, None)
    assert result.on_show_ecd_spectrum(object()) == (None, None)


def test_show_ecd_spectrum_builds_signed_table_and_figure():
    pytest.importorskip("gradio")
    import plotly.graph_objects as go
    from orca_webui import result

    peak_df, fig = result.on_show_ecd_spectrum(FakeTDDFTData(with_rotats=True))

    assert list(peak_df.columns) == ["Wavelength (nm)", "Rotatory Strength (1e-40 cgs)",
                                     "Oscillator Strength", "Transitions", "Symmetry"]
    # Signs must survive into the table; a negative Cotton effect is real information.
    assert (peak_df["Rotatory Strength (1e-40 cgs)"] < 0).any()
    assert (peak_df["Rotatory Strength (1e-40 cgs)"] > 0).any()
    # Sorted by wavelength, like the absorption table.
    assert list(peak_df["Wavelength (nm)"]) == sorted(peak_df["Wavelength (nm)"])
    assert isinstance(fig, go.Figure)


# --- file dropdowns are sorted alphabetically -------------------------------
# os.listdir returns an arbitrary order, so the lists are sorted before display.

def test_file_dropdowns_sorted_alphabetically():
    pytest.importorskip("gradio")
    from orca_webui import calculation, result

    unsorted_files = [
        "zinc.xyz", "Benzene.pdb", "alpha.mol", "opt.log",
        "run_b.inp", "Run_a.inp", "acetone.log",
    ]

    structures, inputs = calculation.on_working_directory_file_list_change(
        unsorted_files, "")
    assert structures["choices"] == [
        "acetone.log", "alpha.mol", "Benzene.pdb", "opt.log", "zinc.xyz"]
    # Case-insensitive: "Run_a.inp" must not sort ahead of every lowercase name.
    assert inputs["choices"] == ["Run_a.inp", "run_b.inp"]
    # The pre-selected value is the first entry of the sorted list.
    assert structures["value"] == "acetone.log"
    assert inputs["value"] == "Run_a.inp"

    results = result.on_working_directory_file_list_change(unsorted_files)
    assert results["choices"] == ["acetone.log", "opt.log"]
    assert results["value"] == "acetone.log"


def test_input_file_selection_survives_sorting():
    """A name typed in the file-name box still wins over the alphabetical first."""
    pytest.importorskip("gradio")
    from orca_webui import calculation

    _structures, inputs = calculation.on_working_directory_file_list_change(
        ["zzz.inp", "aaa.inp"], "zzz")
    assert inputs["choices"] == ["aaa.inp", "zzz.inp"]
    assert inputs["value"] == "zzz.inp"


# --- calculation.on_run_calculation: the .log is not clobbered --------------
# sys.executable stands in for the ORCA binary: the handler runs "<exe> <input file>",
# so an .inp holding valid Python runs and exits 0 on any platform.

def test_run_calculation_backs_up_previous_log(tmp_path):
    pytest.importorskip("gradio")
    from orca_webui import calculation

    (tmp_path / "run.inp").write_text('print("SECOND RUN")\n')
    (tmp_path / "run.log").write_text("FIRST RUN RESULTS")

    status, _files = calculation.on_run_calculation(str(tmp_path), "run.inp", sys.executable)

    assert "color:red" not in status
    assert (tmp_path / "run.log.bak").read_text().strip() == "FIRST RUN RESULTS"
    assert "SECOND RUN" in (tmp_path / "run.log").read_text()
    assert "run.log.bak" in status  # the rescue is reported to the user


def test_run_calculation_makes_no_backup_on_first_run(tmp_path):
    pytest.importorskip("gradio")
    from orca_webui import calculation

    (tmp_path / "run.inp").write_text('print("FIRST RUN")\n')

    status, _files = calculation.on_run_calculation(str(tmp_path), "run.inp", sys.executable)

    assert "color:red" not in status
    assert not (tmp_path / "run.log.bak").exists()
    assert ".bak" not in status


def test_run_calculation_reports_missing_executable(tmp_path):
    pytest.importorskip("gradio")
    from orca_webui import calculation

    (tmp_path / "run.inp").write_text("dummy\n")

    status, _files = calculation.on_run_calculation(
        str(tmp_path), "run.inp", os.path.join(os.sep, "nonexistent", "orca"))

    assert "color:red" in status
    # Bailing out before the run means the previous log is left untouched.
    assert not (tmp_path / "run.log").exists()


# --- working_directory.on_view_structure_file: atom labels ------------------

ETHANOL_XYZ = """0 1
C 0.9026489597488216 -0.05000028588257832 0.021394609774926115
C -0.5742031126818012 -0.385505872155332 0.04058986810948101
O -1.3316454554006916 0.7756870576731353 0.3485124118176326
H 1.2173201961756799 0.35934440949711477 0.9869668094581394
H 1.116183702004607 0.7142026160257686 -0.7331640213952032
H 1.5014612719973133 -0.9382867335837253 -0.19824701389828622
H -0.7850341567874283 -1.1421184407750054 0.8019259675018772
H -0.9017076062007929 -0.7714452849250665 -0.9295025768488222
H -1.1450237988557044 1.4381225341256896 -0.33847605451974583
"""


def test_structure_viewer_labels_atoms_by_file_position(tmp_path, monkeypatch):
    """Labels must number atoms by their position in the file, not per element.

    NGL's "atomname" labels come from the PDB block nglview generates, and RDKit
    numbers those per element -- which rendered ethanol's third atom as "O1" (the first
    oxygen) instead of "O2" (the third coordinate line of the generated .inp).
    Numbering is 0-based to match ORCA's own nucleus indices.
    """
    pytest.importorskip("gradio")
    pytest.importorskip("nglview")
    from orca_webui import working_directory

    (tmp_path / "conformer_1.xyz").write_text(ETHANOL_XYZ)
    (tmp_path / "static").mkdir()
    monkeypatch.chdir(tmp_path)

    html = working_directory.on_view_structure_file(str(tmp_path), "conformer_1.xyz")
    assert html is not None

    doc = (tmp_path / "static" / "structure.html").read_text()
    for expected in ["C0", "C1", "O2", "H3", "H4", "H5", "H6", "H7", "H8"]:
        assert f'"{expected}"' in doc, f"missing label {expected}"
    # The per-element PDB naming that caused the bug must not be in use.
    assert '"atomname"' not in doc


def test_viewer_atom_order_matches_generated_input_file(tmp_path):
    """The labels are only meaningful if they match the .inp coordinate block."""
    pytest.importorskip("gradio")
    from orca_webui import utils

    structure = tmp_path / "conformer_1.xyz"
    structure.write_text(ETHANOL_XYZ)
    mol = utils.mol_from_structure_file(str(structure))

    inp = tmp_path / "job"
    utils.write_sp_orca_input(mol, str(inp), "DFT", "", functional="B3LYP", basis="3-21G")

    block = (tmp_path / "job.inp").read_text().split("* xyz")[1].split("*")[0]
    written = [line.split()[0] for line in block.strip().splitlines()[1:]]
    labelled = [atom.GetSymbol() for atom in mol.GetAtoms()]

    assert written == labelled == ["C", "C", "O", "H", "H", "H", "H", "H", "H"]


def test_run_calculation_reports_nonzero_exit(tmp_path):
    pytest.importorskip("gradio")
    from orca_webui import calculation

    (tmp_path / "run.inp").write_text("raise SystemExit(1)\n")

    status, _files = calculation.on_run_calculation(str(tmp_path), "run.inp", sys.executable)

    assert "color:red" in status
    assert "exited with code 1" in status
