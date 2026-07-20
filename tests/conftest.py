"""Shared pytest fixtures and path setup for the orca-webui test suite.

The application modules live in the ``src/orca_webui`` package, so ``src`` is added to
sys.path here (mirroring ``pythonpath = src`` in pytest.ini) to keep
``from orca_webui import ...`` working no matter where pytest is invoked from.

The synthetic ORCA-log fixtures below deliberately reproduce only the small slices of
real ORCA output that the parsers in utils.py key on (the exact block headers and line
formats documented in those functions). They are hand-written on purpose: the real
outputs under data/ are gitignored and must not be relied on as test fixtures.
"""
import os
import sys

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))


@pytest.fixture
def water_mol():
    """A 3D water molecule (has one conformer), for ORCA coordinate/input tests."""
    mol = Chem.AddHs(Chem.MolFromSmiles("O"))
    assert AllChem.EmbedMolecule(mol, randomSeed=0xF00D) == 0
    return mol


# A tetrahedral methane geometry (Angstrom) in ORCA nucleus order: C then 4 H.
# C-H ~1.09 A (bonded); H-H ~1.78 A (not bonded). All four H are equivalent.
METHANE_SYMBOLS = ["C", "H", "H", "H", "H"]
METHANE_COORDS = [
    [0.000000, 0.000000, 0.000000],
    [0.629118, 0.629118, 0.629118],
    [-0.629118, -0.629118, 0.629118],
    [-0.629118, 0.629118, -0.629118],
    [0.629118, -0.629118, -0.629118],
]


@pytest.fixture
def methane_geometry():
    return list(METHANE_SYMBOLS), [list(c) for c in METHANE_COORDS]


@pytest.fixture
def nmr_log_text():
    """Minimal ORCA NMR output covering every block the parsers scan for.

    - CARTESIAN COORDINATES (ANGSTROEM): consumed by parse_orca_geometry
    - CHEMICAL SHIELDINGS / CHEMICAL SHIELDING SUMMARY: parse_nmr_shielding_constants
    - NUCLEUS A/B + J[..](Total) iso=: parse_nmr_jcouplings
    """
    return """\
------------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
------------------------------------
  C    0.000000    0.000000    0.000000
  H    0.629118    0.629118    0.629118
  H   -0.629118   -0.629118    0.629118

CHEMICAL SHIELDINGS (ppm)
=========================

CHEMICAL SHIELDING SUMMARY (ppm)
--------------------------------------
  Nucleus  Element    Isotropic     Anisotropy
  0   C     190.1234    10.0
  1   H      30.5000     5.0
  2   H      31.5000     5.0

 NUCLEUS A = C    0 NUCLEUS B = H    1
 something something
 J[0,1](Total)   37.220   46.839   48.366  iso=   44.142
 NUCLEUS A = H    1 NUCLEUS B = H    2
 K[1,2](Total)   1.0      1.0      1.0     iso=   999.0
 J[1,2](Total)   6.900    7.100    7.000   iso=    7.000
"""
