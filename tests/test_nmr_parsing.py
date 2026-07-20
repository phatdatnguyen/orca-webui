"""Tests for the ORCA log text parsers used by the NMR pipeline."""
from orca_webui import utils


def test_parse_shielding_constants(tmp_path, nmr_log_text):
    log = tmp_path / "nmr.log"
    log.write_text(nmr_log_text)
    data = utils.parse_nmr_shielding_constants(str(log))
    assert data == [
        (0, "C", 190.1234),
        (1, "H", 30.5000),
        (2, "H", 31.5000),
    ]


def test_parse_shielding_returns_empty_when_absent(tmp_path):
    log = tmp_path / "empty.log"
    log.write_text("no shielding data here\n")
    assert utils.parse_nmr_shielding_constants(str(log)) == []


def test_parse_jcouplings(tmp_path, nmr_log_text):
    log = tmp_path / "nmr.log"
    log.write_text(nmr_log_text)
    jc = utils.parse_nmr_jcouplings(str(log))
    # Stored as (idx_a, elem_a, idx_b, elem_b, J_iso)
    assert (0, "C", 1, "H", 44.142) in jc
    assert (1, "H", 2, "H", 7.000) in jc


def test_parse_jcouplings_ignores_reduced_K_block(tmp_path, nmr_log_text):
    """Only J[..](Total) lines count; the reduced-coupling K[..] line (iso=999) must not."""
    log = tmp_path / "nmr.log"
    log.write_text(nmr_log_text)
    jc = utils.parse_nmr_jcouplings(str(log))
    assert all(val != 999.0 for *_rest, val in jc)


def test_parse_orca_geometry(tmp_path, nmr_log_text):
    log = tmp_path / "nmr.log"
    log.write_text(nmr_log_text)
    symbols, coords = utils.parse_orca_geometry(str(log))
    assert symbols == ["C", "H", "H"]
    assert coords.shape == (3, 3)
    assert abs(coords[1][0] - 0.629118) < 1e-6


def test_parse_orca_geometry_keeps_last_block(tmp_path):
    text = (
        "CARTESIAN COORDINATES (ANGSTROEM)\n"
        "---\n"
        "C 0.0 0.0 0.0\n"
        "\n"
        "CARTESIAN COORDINATES (ANGSTROEM)\n"
        "---\n"
        "O 1.0 1.0 1.0\n"
        "H 2.0 2.0 2.0\n"
    )
    log = tmp_path / "geo.log"
    log.write_text(text)
    symbols, coords = utils.parse_orca_geometry(str(log))
    assert symbols == ["O", "H"]  # last block wins


def test_calculate_chemical_shifts():
    shielding = [(0, "C", 190.0), (1, "H", 30.0), (2, "H", 31.0)]
    refs = {"C": 200.0, "H": 32.0}
    shifts = utils.calculate_chemical_shifts(shielding, refs)
    assert shifts == [(0, "C", 10.0), (1, "H", 2.0), (2, "H", 1.0)]


def test_calculate_chemical_shifts_skips_unknown_element():
    shielding = [(0, "F", 100.0), (1, "H", 30.0)]
    refs = {"H": 32.0}
    shifts = utils.calculate_chemical_shifts(shielding, refs)
    assert shifts == [(1, "H", 2.0)]  # F skipped, no reference
