"""Tests for molecule I/O helpers: xyz round-trip, bond perception, file listing."""
import os

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from orca_webui import utils


# --- get_files_in_working_directory -----------------------------------------

def test_get_files_excludes_zone_identifier(tmp_path):
    (tmp_path / "a.xyz").write_text("x")
    (tmp_path / "b.inp").write_text("y")
    (tmp_path / "a.xyz.Zone.Identifier").write_text("junk")
    files = utils.get_files_in_working_directory(str(tmp_path))
    assert "a.xyz" in files
    assert "b.inp" in files
    assert all(not f.endswith("Zone.Identifier") for f in files)


# --- conformer_to_xyz_file + mol_from_xyz_file round trip --------------------

def test_xyz_roundtrip_preserves_atoms_and_charge(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("CO"))  # methanol: C, O, 4 H
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    path = str(tmp_path / "mol.xyz")
    utils.conformer_to_xyz_file(mol, 0, path, charge=-1, multiplicity=2)

    parsed, charge, mult = utils.mol_from_xyz_file(path, return_charge_and_multiplicity=True)
    assert charge == -1
    assert mult == 2
    assert parsed.GetNumAtoms() == mol.GetNumAtoms()
    # element multiset preserved
    orig = sorted(a.GetSymbol() for a in mol.GetAtoms())
    got = sorted(a.GetSymbol() for a in parsed.GetAtoms())
    assert orig == got


def test_xyz_coordinates_preserved(tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles("O"))
    assert AllChem.EmbedMolecule(mol, randomSeed=2) == 0
    path = str(tmp_path / "w.xyz")
    utils.conformer_to_xyz_file(mol, 0, path)
    parsed = utils.mol_from_xyz_file(path)

    conf_a = mol.GetConformer()
    conf_b = parsed.GetConformer()
    for i in range(mol.GetNumAtoms()):
        pa, pb = conf_a.GetAtomPosition(i), conf_b.GetAtomPosition(i)
        assert abs(pa.x - pb.x) < 1e-4
        assert abs(pa.y - pb.y) < 1e-4
        assert abs(pa.z - pb.z) < 1e-4


def test_xyz_defaults_charge_multiplicity_when_absent(tmp_path):
    # An xyz whose first line is a bare count (not "charge multiplicity")
    path = tmp_path / "plain.xyz"
    path.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n")
    mol, charge, mult = utils.mol_from_xyz_file(str(path), return_charge_and_multiplicity=True)
    assert (charge, mult) == (0, 1)
    assert mol.GetNumAtoms() == 3


# --- add_bonds ---------------------------------------------------------------

def test_add_bonds_infers_connectivity():
    # methanol geometry: RDKit builds bonds; strip them, re-perceive from distance.
    mol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    assert AllChem.EmbedMolecule(mol, randomSeed=3) == 0

    # Rebuild a bond-free mol carrying only atoms + conformer.
    rw = Chem.RWMol()
    for atom in mol.GetAtoms():
        rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
    rw.AddConformer(mol.GetConformer(), assignId=True)
    bare = rw.GetMol()
    assert bare.GetNumBonds() == 0

    bonded = utils.add_bonds(bare)
    # methanol has 5 bonds (C-O, O-H, 3x C-H)
    assert bonded.GetNumBonds() == 5


def test_add_bonds_charges_tetravalent_boron(tmp_path):
    """A 4-coordinate boron (borate) must receive a -1 formal charge so the molecule
    passes RDKit sanitization instead of raising "Explicit valence for atom B, 4, is
    greater than permitted". This is the BF2 bridge case in BODIPY dyes; modeled here
    with a tetrahedral BF4- so the test needs no external file.
    """
    path = tmp_path / "bf4.xyz"
    path.write_text(
        "0 1\nBF4\n"
        "B  0.000  0.000  0.000\n"
        "F  0.810  0.810  0.810\n"
        "F -0.810 -0.810  0.810\n"
        "F -0.810  0.810 -0.810\n"
        "F  0.810 -0.810 -0.810\n"
    )
    mol = utils.add_bonds(utils.mol_from_xyz_file(str(path)))
    boron = next(a for a in mol.GetAtoms() if a.GetSymbol() == "B")
    assert boron.GetDegree() == 4
    assert boron.GetFormalCharge() == -1
    Chem.SanitizeMol(mol)  # must not raise


def test_add_bonds_leaves_normal_valences_uncharged():
    """Ordinary (non-over-coordinated) atoms keep a zero formal charge."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CO"))  # methanol: normal C, O valences
    assert AllChem.EmbedMolecule(mol, randomSeed=5) == 0
    rw = Chem.RWMol()
    for atom in mol.GetAtoms():
        rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
    rw.AddConformer(mol.GetConformer(), assignId=True)
    bonded = utils.add_bonds(rw.GetMol())
    assert all(a.GetFormalCharge() == 0 for a in bonded.GetAtoms())


def test_add_bonds_preserves_conformer_geometry():
    """add_bonds must keep the input geometry (exactly one conformer, same coords).

    Regression guard: the structure viewer only re-embeds when GetNumConformers()==0,
    so a parsed file keeps its REAL coordinates. If add_bonds dropped the conformer,
    the viewer would fall back to a fresh (all-single-bond) embedding that puckers
    planar/aromatic systems into a bogus sp3-looking geometry.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))  # benzene: planar
    assert AllChem.EmbedMolecule(mol, randomSeed=11) == 0

    rw = Chem.RWMol()
    for atom in mol.GetAtoms():
        rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
    rw.AddConformer(mol.GetConformer(), assignId=True)
    bonded = utils.add_bonds(rw.GetMol())

    assert bonded.GetNumConformers() == 1  # geometry retained, viewer won't re-embed
    src, dst = mol.GetConformer(), bonded.GetConformer()
    for i in range(mol.GetNumAtoms()):
        a, b = src.GetAtomPosition(i), dst.GetAtomPosition(i)
        assert abs(a.x - b.x) < 1e-6 and abs(a.y - b.y) < 1e-6 and abs(a.z - b.z) < 1e-6


# --- mol_from_orca_file (requires cclib) ------------------------------------

def test_mol_from_orca_file_requires_cclib_only_when_called():
    """utils must import even without cclib; the dependency is used lazily."""
    # Importing utils already succeeded at module load, which proves the point.
    # If cclib is installed, a bogus path should raise from cclib, not ImportError.
    pytest.importorskip("cclib")
    with pytest.raises(Exception) as excinfo:
        utils.mol_from_orca_file(os.path.join("does", "not", "exist.log"))
    assert not isinstance(excinfo.value, ImportError)


# --- mol_from_orca_geometry (ORCA nucleus order, no cclib) -------------------

def test_mol_from_orca_geometry_matches_orca_nucleus_numbering(tmp_path, nmr_log_text):
    """Atom indices must line up with the nucleus indices ORCA reports.

    The structure viewer labels atoms "<element> (index=<n>)" so a user can pick
    J-coupling partners by eye; those numbers are only meaningful if index n really is
    ORCA nucleus n. cclib reorders atoms, hence this cclib-free path.
    """
    log = tmp_path / "nmr.log"
    log.write_text(nmr_log_text)

    mol = utils.mol_from_orca_geometry(str(log))
    viewer = [(a.GetIdx(), a.GetSymbol()) for a in mol.GetAtoms()]
    orca = [(idx, elem) for idx, elem, _ in utils.parse_nmr_shielding_constants(str(log))]

    assert viewer == orca == [(0, "C"), (1, "H"), (2, "H")]


def test_mol_from_orca_geometry_keeps_coordinates(tmp_path, nmr_log_text):
    """The conformer must carry the log geometry, not a re-embedded one."""
    log = tmp_path / "nmr.log"
    log.write_text(nmr_log_text)

    mol = utils.mol_from_orca_geometry(str(log))
    _, coords = utils.parse_orca_geometry(str(log))
    conf = mol.GetConformer()

    assert mol.GetNumConformers() == 1
    for i, (x, y, z) in enumerate(coords):
        pos = conf.GetAtomPosition(i)
        assert abs(pos.x - x) < 1e-6 and abs(pos.y - y) < 1e-6 and abs(pos.z - z) < 1e-6
