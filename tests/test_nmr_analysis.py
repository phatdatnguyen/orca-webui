"""Tests for the bond-perception, symmetry-equivalence, and multiplet logic."""
from orca_webui import utils


# --- perceive_bonds / morgan_ranks ------------------------------------------

def test_perceive_bonds_methane(methane_geometry):
    symbols, coords = methane_geometry
    adjacency = utils.perceive_bonds(coords, symbols)
    # C (index 0) bonded to all four H; each H bonded only to C.
    assert adjacency[0] == {1, 2, 3, 4}
    for h in range(1, 5):
        assert adjacency[h] == {0}


def test_morgan_ranks_group_equivalent_hydrogens(methane_geometry):
    symbols, coords = methane_geometry
    _, ranks = utils.morgan_ranks(coords, symbols)
    h_ranks = {ranks[i] for i in range(1, 5)}
    assert len(h_ranks) == 1               # all four H equivalent
    assert ranks[0] not in h_ranks          # C distinct from H


# --- equivalence maps --------------------------------------------------------

def test_hydrogen_equivalence_methane(methane_geometry):
    symbols, coords = methane_geometry
    equiv = utils.compute_hydrogen_equivalence(coords, symbols)
    assert set(equiv.keys()) == {1, 2, 3, 4}   # only hydrogens keyed
    assert len(set(equiv.values())) == 1        # all four collapse to one group


def test_carbon_equivalence_only_keys_carbons(methane_geometry):
    symbols, coords = methane_geometry
    equiv = utils.compute_carbon_equivalence(coords, symbols)
    assert set(equiv.keys()) == {0}


# --- build_nmr_peak_table + multiplicity_label -------------------------------

def _ethyl_inputs():
    """A synthetic CH3-CH2 ethyl fragment: group A = 3 H, group B = 2 H, J = 7 Hz.

    Uses explicit shifts_data + equivalence + jcouplings so the multiplet logic is
    tested in isolation (no geometry / DFT numbers involved).
    """
    shifts_data = [
        (0, "H", 1.20), (1, "H", 1.20), (2, "H", 1.20),   # CH3
        (3, "H", 3.60), (4, "H", 3.60),                    # CH2
    ]
    equivalence = {0: ("grp", "A"), 1: ("grp", "A"), 2: ("grp", "A"),
                   3: ("grp", "B"), 4: ("grp", "B")}
    jcouplings = []
    for a in (0, 1, 2):
        for b in (3, 4):
            jcouplings.append((a, "H", b, "H", 7.0))
    return shifts_data, equivalence, jcouplings


def test_peak_table_averages_equivalent_and_assigns_multiplicity():
    shifts_data, equivalence, jcouplings = _ethyl_inputs()
    peaks = utils.build_nmr_peak_table(shifts_data, "H",
                                       jcouplings=jcouplings, equivalence=equivalence)
    assert len(peaks) == 2
    peaks.sort(key=lambda p: p["shift"])
    ch3, ch2 = peaks

    # CH3: 3 protons, split by the 2-proton partner -> triplet
    assert ch3["count"] == 3
    assert abs(ch3["shift"] - 1.20) < 1e-9
    assert len(ch3["couplings"]) == 1
    assert ch3["couplings"][0]["n"] == 2
    assert abs(ch3["couplings"][0]["J"] - 7.0) < 1e-9
    assert utils.multiplicity_label(ch3["couplings"]) == "t"

    # CH2: 2 protons, split by the 3-proton partner -> quartet
    assert ch2["count"] == 2
    assert ch2["couplings"][0]["n"] == 3
    assert utils.multiplicity_label(ch2["couplings"]) == "q"


def test_peak_table_drops_negligible_coupling():
    shifts_data, equivalence, _ = _ethyl_inputs()
    tiny = [(a, "H", b, "H", 0.1) for a in (0, 1, 2) for b in (3, 4)]
    peaks = utils.build_nmr_peak_table(shifts_data, "H",
                                       jcouplings=tiny, equivalence=equivalence,
                                       j_threshold=0.5)
    for p in peaks:
        assert p["couplings"] == []       # 0.1 Hz < threshold -> unresolved singlet
        assert utils.multiplicity_label(p["couplings"]) == "s"


def test_peak_table_empty_for_missing_element():
    shifts_data = [(0, "H", 1.0)]
    assert utils.build_nmr_peak_table(shifts_data, "C") == []


def test_multiplicity_label_singlet_when_no_coupling():
    assert utils.multiplicity_label([]) == "s"


def test_multiplicity_label_multiple_groups_concatenate():
    # n=2 -> 't', n=1 -> 'd', concatenated in list order.
    couplings = [{"n": 2, "J": 7.0}, {"n": 1, "J": 3.0}]
    assert utils.multiplicity_label(couplings) == "td"


def test_multiplicity_label_large_group_is_m():
    # Group sizes outside the lookup table fall back to 'm'.
    assert utils.multiplicity_label([{"n": 99, "J": 7.0}]) == "m"
