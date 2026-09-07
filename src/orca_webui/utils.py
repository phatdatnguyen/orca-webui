"""Core, UI-free engine for orca-webui.

Groups two responsibilities used by every Gradio tab:

1. Molecule I/O and ORCA input-file generation (``mol_from_*``, ``add_bonds``,
   ``write_*_orca_input`` and helpers).
2. Spectrum simulation and ORCA-output text parsing, including the NMR
   shielding/coupling -> chemical-shift -> multiplet pipeline.

Everything here is pure Python with no Gradio dependency so it can be tested and
reused directly.
"""
from __future__ import annotations

import os
import re
from math import comb
from typing import TextIO

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import plotly.graph_objects as go


def resolve_working_file(working_directory_path: str, file_name: str) -> str:
    """Resolve a single file name inside an existing working directory.

    UI selections and editable output names must not escape the chosen directory,
    including through a symlink. Return an absolute path for subprocess callers.
    """
    if not working_directory_path or not os.path.isdir(working_directory_path):
        raise ValueError("Please open an existing working directory.")
    if (not isinstance(file_name, str) or not file_name.strip()
            or file_name in {'.', '..'}
            or any(c in file_name for c in ('/', '\\', '\x00'))
            or os.path.isabs(file_name) or os.path.splitdrive(file_name)[0]):
        raise ValueError("Please use a file name without directory separators.")
    directory = os.path.realpath(working_directory_path)
    file_path = os.path.join(directory, file_name)
    if os.path.commonpath([directory, os.path.realpath(file_path)]) != directory:
        raise ValueError("The selected file points outside the working directory.")
    return file_path


def get_files_in_working_directory(working_directory_path: str) -> list[str]:
    """Return the file names in ``working_directory_path``, or ``[]`` if it is unset.

    Windows ``Zone.Identifier`` alternate-data-stream files (created when files are
    downloaded) are filtered out so they never appear in the UI file list.

    The ``None`` guard matters: the shared path ``gr.State`` starts empty and is reset
    to ``None`` when a directory fails to open, and ``os.listdir(None)`` quietly lists
    the *current* directory instead of raising — which surfaced as callers building
    paths out of a ``None`` directory.
    """
    if not working_directory_path:
        return []
    try:
        with os.scandir(working_directory_path) as entries:
            files = []
            for entry in entries:
                if entry.name.endswith('Zone.Identifier') or not entry.is_file():
                    continue
                try:
                    resolve_working_file(working_directory_path, entry.name)
                except ValueError:
                    continue
                files.append(entry.name)
        return sorted(files, key=str.lower)
    except OSError:
        # Error handlers also refresh this state; a directory deleted outside the
        # app must not raise again while they are reporting the original error.
        return []

def conformer_to_xyz_file(mol: Chem.Mol, conf_id: int, file_path: str,
                          charge: int = 0, multiplicity: int = 1) -> None:
    """Write one conformer of ``mol`` to an XYZ file.

    The first line is ``"<charge> <multiplicity>"`` (the convention this app's XYZ
    reader expects) rather than the standard atom count, followed by one
    ``"<symbol> <x> <y> <z>"`` line per atom taken from conformer ``conf_id``.
    """
    # Get atom information
    atoms = mol.GetAtoms()
    xyz_lines = []
    for atom in atoms:
        pos = mol.GetConformer(conf_id).GetAtomPosition(atom.GetIdx())
        xyz_lines.append(f"{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}")

    # Construct the XYZ string
    xyz_string = f"{charge} {multiplicity}\n" + "\n".join(xyz_lines)
    
    with open(file_path, 'w') as file:
        file.write(xyz_string)

def add_bonds(mol: Chem.Mol, bond_factor: float = 1.25) -> Chem.Mol:
    """Return a copy of ``mol`` with bonds inferred from interatomic distances.

    Structure files such as XYZ carry only atoms and coordinates, so RDKit sees no
    connectivity. Two atoms are joined by a single bond when their distance is below
    ``(r_cov_i + r_cov_j) * bond_factor`` using RDKit's covalent radii. Bond orders
    are not perceived (all bonds are single); this is enough for visualization and
    for feeding coordinates into ORCA input generation.

    Formal charges are assigned to hypervalent main-group atoms so the result passes
    ``Chem.SanitizeMol`` — currently 4-coordinate boron (see ``_assign_formal_charges``),
    e.g. the ``BF2`` bridge in BODIPY dyes. These are per-atom bookkeeping charges only;
    the overall charge used for an ORCA calculation is set separately by the user.
    """
    # Create a new empty molecule
    mol_new = Chem.RWMol()

    # Add atoms
    for atom in mol.GetAtoms():
        mol_new.AddAtom(atom)

    # Add conformer
    conf = mol.GetConformer()
    mol_new.AddConformer(conf)

    # Add bonds based on covalent radii
    for i in range(mol.GetNumAtoms()):
        for j in range(i + 1, mol.GetNumAtoms()):
            dist = np.linalg.norm(np.array(conf.GetAtomPosition(i)) - np.array(conf.GetAtomPosition(j)))
            pt = Chem.GetPeriodicTable()
            r_cov_i = pt.GetRcovalent(mol.GetAtomWithIdx(i).GetSymbol())
            r_cov_j = pt.GetRcovalent(mol.GetAtomWithIdx(j).GetSymbol())
            if dist < (r_cov_i + r_cov_j) * bond_factor:
                mol_new.AddBond(i, j, Chem.BondType.SINGLE)

    # Reconcile valences that would otherwise fail sanitization (e.g. borate boron).
    _assign_formal_charges(mol_new)

    # Convert to a regular Mol object
    mol_new = mol_new.GetMol()

    return mol_new

# Neutral bonding valence for main-group elements whose over-coordination (with the
# all-single-bond perception used here) implies a formal charge. degree - valence gives
# the charge: 4-coordinate boron -> -1 (borate), a 4th bond on N -> +1 (ammonium), etc.
_NEUTRAL_VALENCE = {'B': 3, 'N': 3, 'O': 2}

def _assign_formal_charges(mol: Chem.RWMol) -> None:
    """Assign formal charges so over-coordinated main-group atoms pass sanitization.

    Only atoms bonded to MORE neighbors than their neutral valence are touched, and the
    charge is the difference (boron gains electrons -> negative; N/O lose -> positive).
    This fixes structures like BODIPY's 4-coordinate borate boron, which RDKit otherwise
    rejects with "Explicit valence for atom B, 4, is greater than permitted."
    """
    for atom in mol.GetAtoms():
        neutral = _NEUTRAL_VALENCE.get(atom.GetSymbol())
        if neutral is None:
            continue
        excess = atom.GetDegree() - neutral
        if excess > 0:
            # Group 13 (B) is electron-deficient: an extra bond makes it anionic;
            # groups 15/16 (N, O) become cationic.
            sign = -1 if atom.GetSymbol() == 'B' else 1
            atom.SetFormalCharge(sign * excess)

def mol_from_xyz_file(file_path: str, return_charge_and_multiplicity: bool = False):
    """Read an XYZ file into a bond-free RDKit ``Mol`` with one conformer.

    Supports standard XYZ (atom count, comment, coordinates) and the legacy
    ``"<charge> <multiplicity>"`` header written by :func:`conformer_to_xyz_file`.
    Standard XYZ comments are never interpreted as atoms. Malformed or truncated
    coordinates raise ``ValueError`` rather than silently dropping atoms.

    Returns the ``Mol`` alone, or ``(mol, charge, multiplicity)`` when
    ``return_charge_and_multiplicity`` is True.
    """
    with open(file_path, 'r') as file:
        xyz_string = file.read()

    lines = xyz_string.splitlines()
    element_symbols = {Chem.GetPeriodicTable().GetElementSymbol(i) for i in range(1, 119)}
    if not lines:
        raise ValueError("The XYZ file is empty.")
    header = lines[0].split()
    charge, multiplicity = 0, 1
    atom_count = None
    if len(header) == 1:
        try:
            atom_count = int(header[0])
        except ValueError as exc:
            raise ValueError("The XYZ header must contain an atom count or charge and multiplicity.") from exc
        if atom_count < 1:
            raise ValueError("The XYZ atom count must be positive.")
        coordinate_lines = lines[2:]
    elif len(header) == 2:
        try:
            charge, multiplicity = map(int, header)
        except ValueError as exc:
            raise ValueError("The XYZ charge and multiplicity must be integers.") from exc
        if multiplicity < 1:
            raise ValueError("The XYZ multiplicity must be positive.")
        coordinate_lines = lines[1:]
        # Older files may include a comment after the charge/multiplicity header.
        if coordinate_lines:
            first = coordinate_lines[0].split()
            is_coordinate = bool(first) and first[0].capitalize() in element_symbols
            if len(first) == 4:
                try:
                    tuple(map(float, first[1:]))
                    is_coordinate = True  # Unknown element symbols must be rejected below.
                except ValueError:
                    pass
            if not is_coordinate:
                coordinate_lines = coordinate_lines[1:]
    else:
        raise ValueError("The XYZ header must contain an atom count or charge and multiplicity.")

    coordinate_lines = [line for line in coordinate_lines if line.strip()]
    if atom_count is not None and len(coordinate_lines) != atom_count:
        raise ValueError(f"XYZ atom count is {atom_count}, but found {len(coordinate_lines)} coordinate lines.")
    if not coordinate_lines:
        raise ValueError("The XYZ file contains no atomic coordinates.")

    mol = Chem.RWMol()
    conf = Chem.Conformer(len(coordinate_lines))
    for atom_idx, line in enumerate(coordinate_lines):
        parts = line.split()
        if len(parts) != 4:
            raise ValueError(f"Invalid XYZ coordinate line: {line!r}")
        try:
            atom = Chem.Atom(parts[0].capitalize())
            coords = tuple(map(float, parts[1:]))
        except (ValueError, RuntimeError) as exc:
            raise ValueError(f"Invalid XYZ coordinate line: {line!r}") from exc
        if atom.GetAtomicNum() == 0 or not np.all(np.isfinite(coords)):
            raise ValueError(f"Invalid XYZ coordinate line: {line!r}")
        mol.AddAtom(atom)
        conf.SetAtomPosition(atom_idx, coords)
    mol.AddConformer(conf)
                
    # Convert to a regular Mol object
    mol = mol.GetMol()
    
    if return_charge_and_multiplicity:
        return mol, charge, multiplicity
    else:
        return mol

def mol_from_orca_file(filename: str) -> Chem.Mol:
    """Read the final geometry from an ORCA output file into a bond-free ``Mol``.

    Uses cclib to parse the log and takes the last set of coordinates
    (``atomcoords[-1]``), i.e. the optimized geometry for an Opt job. The returned
    molecule has atoms and one conformer but no bonds; pass it through
    :func:`add_bonds` if connectivity is needed.
    """
    # cclib is imported lazily: it is only needed to parse ORCA/Gaussian output
    # files, so the rest of utils (input generation, spectra, NMR text parsing)
    # stays importable without this heavy optional dependency.
    import cclib.io

    # Parse the output file using cclib
    parser = cclib.io.ccopen(filename)
    if parser is None:
        raise ValueError("The output file format could not be recognized.")
    data = parser.parse()

    # Extract atomic numbers and coordinates
    atomic_numbers = data.atomnos  # List of atomic numbers
    coordinates = data.atomcoords[-1]  # Final geometry (last set of coordinates)

    # Create an empty editable RDKit molecule
    mol = Chem.RWMol()

    # Add atoms to the molecule
    for atomic_num in atomic_numbers:
        atom = Chem.Atom(int(atomic_num))
        mol.AddAtom(atom)

    # Create a conformer to store atomic positions
    conf = Chem.Conformer(len(atomic_numbers))
    for i, (x, y, z) in enumerate(coordinates):
        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf)

    # Convert RWMol to Mol before adding bonds
    mol = mol.GetMol()

    return mol

def mol_from_orca_geometry(logfile: str) -> Chem.Mol:
    """Read the final ORCA geometry into a bond-free ``Mol`` in ORCA nucleus order.

    :func:`mol_from_orca_file` goes through cclib, which *reorders* atoms; this builds
    the molecule straight from the log's last ``CARTESIAN COORDINATES`` block instead
    (via :func:`parse_orca_geometry`). Atom indices of the returned ``Mol`` are then
    exactly the indices ORCA uses to number NMR nuclei, so they are safe to show to a
    user who has to name nuclei (e.g. for spin-spin couplings). Like
    :func:`mol_from_orca_file` the result has no bonds; pass it through
    :func:`add_bonds` if connectivity is needed.
    """
    symbols, coords = parse_orca_geometry(logfile)
    if not symbols:
        raise ValueError("No ORCA Cartesian coordinates were found in the output file.")

    mol = Chem.RWMol()
    for symbol in symbols:
        mol.AddAtom(Chem.Atom(symbol.capitalize()))

    conf = Chem.Conformer(len(symbols))
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf)

    return mol.GetMol()

def mol_from_structure_file(file_path: str) -> Chem.Mol:
    """Load any structure file the UI accepts into a sanitized ``Mol`` with a conformer.

    Shared by the structure viewer and ORCA input generation so that both see the same
    atoms in the same order. That ordering is load-bearing: the viewer labels each atom
    with its index, and :func:`get_orca_coordinates` writes the ``* xyz`` block in that
    same index order, so two different loaders here would make the on-screen labels
    disagree with the generated input file.

    ``.log`` files are read straight from the log's coordinate block via
    :func:`mol_from_orca_geometry`; :func:`mol_from_orca_file` (cclib) *reorders* atoms
    and is only used as a fallback when no coordinate block can be parsed.
    """
    file_path = os.fspath(file_path)
    suffix = os.path.splitext(file_path)[1].lower()
    if suffix == '.pdb':
        mol = Chem.MolFromPDBFile(file_path, sanitize=False, removeHs=False)
    elif suffix == '.mol':
        mol = Chem.MolFromMolFile(file_path, sanitize=False, removeHs=False)
    elif suffix == '.mol2':
        mol = Chem.MolFromMol2File(file_path, sanitize=False, removeHs=False)
    elif suffix == '.xyz':
        mol = add_bonds(mol_from_xyz_file(file_path))
    elif suffix == '.log':
        try:
            mol = add_bonds(mol_from_orca_geometry(file_path))
        except Exception:
            mol = add_bonds(mol_from_orca_file(file_path))
    else:
        raise ValueError(f"Unsupported structure file extension: {suffix}")

    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("The structure file contains no readable atoms.")
    Chem.SanitizeMol(mol)
    # Display/write the file's REAL coordinates. Only embed a fresh conformer if the
    # file somehow carried none — re-embedding otherwise would discard the parsed (e.g.
    # ORCA-optimized) geometry and, because add_bonds perceives every bond as single,
    # pucker planar/aromatic systems into a bogus sp3-looking shape.
    if mol.GetNumConformers() == 0:
        if AllChem.EmbedMolecule(mol) < 0:
            raise ValueError("Could not generate coordinates for the structure.")

    return mol

def get_orca_method_keyword(method_type: str, method_name: str,
                            functional: str, basis: str) -> str:
    """Build the ORCA method/basis keyword string for the ``!`` command line.

    Dispatches on ``method_type``: ``HF``/``MP2``/``CCSD``/``CCSD(T)`` use the type
    name, ``DFT`` uses ``functional``, and anything else (semi-empirical) uses
    ``method_name``. The basis set is appended in every case.
    """
    if method_type == "HF":
        return f'HF {basis}'
    elif method_type == "DFT":
        return f'{functional} {basis}'
    elif method_type in ["MP2", "CCSD", "CCSD(T)"]:
        return f'{method_type} {basis}'
    else:  # method_type == "Semi-empirical":
        return f'{method_name} {basis}'

def get_orca_header(f: TextIO, n_proc: int, memory: float | None,
                    solvation: bool, solvation_model: str | None, solvent: str | None,
                    tddft: bool = False, n_states: int = 10, maxdim: int = 5,
                    root: int | None = None) -> None:
    """Write the ORCA ``%``-block header (parallel, memory, TDDFT, solvation) to ``f``.

    Always writes the ``%pal`` block. ``%maxcore`` (per-core memory in MB) is written
    only when ``memory`` is given (``memory`` is in GB). When ``tddft`` is set, a
    ``%tddft`` block with ``nroots``/``maxdim`` is emitted. ``root`` adds ``iroot``,
    which selects the state whose *gradient* is followed, so it only means anything
    alongside ``Opt``: ``iroot 0`` optimizes the ground state, ``iroot N`` (N >= 1)
    optimizes excited state N. Leave it ``None`` for a single-point excited-state
    calculation, where the keyword has no effect. When
    ``solvation`` is on and both a model and solvent are given, a ``%cpcm`` block is
    written (SMD via ``smd true`` + ``smdsolvent``, otherwise plain CPCM ``solvent``).
    """
    f.write(f'%pal\n  nprocs {n_proc}\nend\n')
    if memory is not None:
        maxcore = int(memory * 1024)
        f.write(f'%maxcore {maxcore}\n')
    if tddft:
        f.write('%tddft\n')
        f.write(f'  nroots {n_states}\n')
        f.write(f'  maxdim {maxdim}\n')
        if root is not None:
            f.write(f'  iroot {root}\n')
        f.write('end\n\n')
    if solvation and solvation_model and solvent:
        f.write('%cpcm\n')
        if solvation_model.lower() == 'smd':
            f.write('  smd true\n')
            f.write(f'  smdsolvent "{solvent}"\n')
        else:
            f.write(f'  solvent "{solvent}"\n')
        f.write('end\n\n')

def get_orca_coordinates(f: TextIO, mol: Chem.Mol, charge: int, multiplicity: int) -> None:
    """Write the ORCA ``* xyz ... *`` coordinate block for ``mol`` to ``f``.

    Coordinates are taken from the molecule's active conformer.
    """
    f.write(f'* xyz {charge} {multiplicity}\n')
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        pos = conf.GetAtomPosition(idx)
        f.write(f'{symbol} {pos.x:>12.6f} {pos.y:>12.6f} {pos.z:>12.6f}\n')
    f.write('*\n')

def write_sp_orca_input(mol: Chem.Mol, file_name: str, method_type: str,
                        method_name: str, functional: str = 'B3LYP',
                        basis: str = '6-31G(d)', charge: int = 0, multiplicity: int = 1,
                        solvation: bool = False, solvation_model: str | None = None,
                        solvent: str | None = None, n_proc: int = 4,
                        memory: float = 2) -> None:
    """Write a single-point (``SP``) ORCA input file.

    ``file_name`` is the path *without* extension; ``.inp`` is appended. See
    :func:`get_orca_method_keyword` and :func:`get_orca_header` for how the method
    line and ``%`` blocks are assembled.
    """
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, method_name, functional, basis)
        f.write(f'! {method} SP\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent)
        get_orca_coordinates(f, mol, charge, multiplicity)

def write_opt_orca_input(mol: Chem.Mol, file_name: str, method_type: str,
                         method_name: str, functional: str = 'B3LYP',
                         basis: str = '6-31G(d,p)', charge: int = 0, multiplicity: int = 1,
                         solvation: bool = False, solvation_model: str | None = None,
                         solvent: str | None = None, n_proc: int = 4,
                         memory: float = 2, excited_state: bool = False,
                         n_states: int = 10, maxdim: int = 5, root: int = 1) -> None:
    """Write a geometry-optimization (``Opt``) ORCA input file.

    Optimizes the ground state by default. With ``excited_state`` set, a ``%tddft``
    block makes ORCA follow the gradient of excited state ``root``, giving that state's
    relaxed geometry -- which is the geometry an emission (fluorescence) spectrum has to
    be computed at, since emission happens after the excited state relaxes.

    An excited-state optimization needs the third derivative of the exchange-correlation
    functional, which ORCA 6.1 does not have natively for every functional: B3LYP and
    wB97X-D4 abort with "try using LibXC instead", while e.g. PBE works.

    ``file_name`` is the path without extension; ``.inp`` is appended.
    """
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, method_name, functional, basis)
        f.write(f'! {method} Opt\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent,
                        tddft=excited_state, n_states=n_states, maxdim=maxdim,
                        root=root if excited_state else None)
        get_orca_coordinates(f, mol, charge, multiplicity)

def write_opt_freq_orca_input(mol: Chem.Mol, file_name: str, method_type: str,
                              method_name: str, functional: str = 'B3LYP',
                              basis: str = '6-31G(d,p)', charge: int = 0,
                              multiplicity: int = 1, solvation: bool = False,
                              solvation_model: str | None = None, solvent: str | None = None,
                              n_proc: int = 4, memory: float = 2) -> None:
    """Write an optimization + frequency (``Opt Freq``) ORCA input file.

    The frequency step yields the vibrational/IR and thermochemistry data shown in
    the Result tab. ``file_name`` is the path without extension; ``.inp`` is appended.
    """
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, method_name, functional, basis)
        f.write(f'! {method} Opt Freq\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent)
        get_orca_coordinates(f, mol, charge, multiplicity)

def write_tddft_orca_input(mol: Chem.Mol, file_name: str, method_type: str,
                           method_name: str, functional: str = 'B3LYP',
                           basis: str = '6-31G(d,p)', n_states: int = 10,
                           maxdim: int = 5, charge: int = 0,
                           multiplicity: int = 1, solvation: bool = False,
                           solvation_model: str | None = None, solvent: str | None = None,
                           n_proc: int = 4, memory: float = 2) -> None:
    """Write a single-point TD-DFT input file (absorption/emission/ECD spectra).

    Adds ``RIJCOSX TightSCF`` to the method line and a ``%tddft`` block (``n_states``
    roots, expansion-space ``maxdim``). Deliberately **no** ``Opt``: this is a single
    point on the geometry it is given, so the spectrum belongs to the structure the user
    selected. Optimizing here would silently discard that geometry, spend a full
    optimization on every spectrum, and make ORCA reject B88-containing functionals such
    as B3LYP -- its input check demands XC third derivatives as soon as ``%tddft`` and
    ``Opt`` appear together. Optimize first with :func:`write_opt_orca_input` (including
    for an excited state, when an emission spectrum is wanted) and run this on the
    resulting geometry.

    ``file_name`` is the path without extension; ``.inp`` is appended.
    """
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, method_name, functional, basis)
        f.write(f'! {method} RIJCOSX TightSCF\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent, tddft=True, n_states=n_states, maxdim=maxdim)
        get_orca_coordinates(f, mol, charge, multiplicity)

def write_nmr_orca_input(mol: Chem.Mol, file_name: str, method_type: str,
                         functional: str = 'B3LYP', basis: str = '6-31G(d,p)',
                         spin_spin_coupling: bool = False, charge: int = 0,
                         multiplicity: int = 1, solvation: bool = False,
                         solvation_model: str | None = None, solvent: str | None = None,
                         n_proc: int = 4, memory: float = 2) -> None:
    """Write an NMR (GIAO shielding) ORCA input file.

    When ``spin_spin_coupling`` is set, an ``%eprnmr`` block requests all H and C
    spin-spin couplings (``ssall``) within ``SpinSpinRThresh`` (3 A, enough for the
    geminal/vicinal H-H pairs that produce 1H multiplets). ``file_name`` is the path
    without extension; ``.inp`` is appended.
    """
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, '', functional, basis)
        f.write(f'! {method} NMR\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent)
        get_orca_coordinates(f, mol, charge, multiplicity)
        if spin_spin_coupling:
            # Request all spin-spin coupling contributions for H and C nuclei.
            # SpinSpinRThresh limits pairs by internuclear distance (Angstrom);
            # ~3 A captures geminal and vicinal H-H pairs needed for 1H multiplets.
            f.write('\n%eprnmr\n')
            f.write('  Nuclei = all H { ssall };\n')
            f.write('  Nuclei = all C { ssall };\n')
            f.write('  SpinSpinRThresh 3.0   # Angstrom\n')
            f.write('  PrintLevel 3\n')
            f.write('end\n\n')
        

def lorentzian_ir(wavenumber, position, intensity, width: float = 10):
    """Lorentzian band shape (peak height = ``intensity`` at ``position``).

    ``wavenumber`` may be a scalar or a numpy array (vectorized over the x-axis).
    ``width`` is the half-width at half-maximum in cm^-1.
    """
    return intensity * (width ** 2) / ((wavenumber - position) ** 2 + width ** 2)

def generate_ir_spectrum_interactive(
    frequencies,
    intensities,
    width: float = 10.0,
    points: int = 4000,
    plot_range: tuple[float, float] | None = None,
    normalize: bool = True,
    transmittance: bool = True,
) -> go.Figure | None:
    """Build an interactive IR spectrum from vibrational frequencies + IR intensities.

    Each mode is broadened with a Lorentzian (:func:`lorentzian_ir`) and summed onto
    a wavenumber grid. Optionally normalized and converted to transmittance
    (``1 - absorbance``, with the x-axis reversed per IR convention). Returns a Plotly
    figure, or ``None`` when no frequencies are supplied.
    """
    # Safety check
    if frequencies is None or len(frequencies) == 0:
        print("No IR frequencies detected.")
        return None

    # Define plotting range
    if plot_range:
        min_wn, max_wn = plot_range
    else:
        min_wn = np.min(frequencies) - 100
        max_wn = np.max(frequencies) + 100

    x = np.linspace(min_wn, max_wn, points)
    spectrum = np.zeros_like(x)

    # Generate Lorentzian broadened spectrum
    annotations = []
    for freq, inten in zip(frequencies, intensities):
        spectrum += lorentzian_ir(x, freq, inten, width)

        annotations.append({
            "freq": freq,
            "intensity": inten,
            "text": f"{freq:.1f}"
        })

    # Normalize if requested
    if normalize and np.max(spectrum) != 0:
        spectrum /= np.max(spectrum)

    # Convert to transmittance if requested
    if transmittance:
        spectrum = 1 - spectrum  # simple inversion

    # Create figure
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=x,
        y=spectrum,
        mode='lines',
        line=dict(color='black'),
        name="IR Spectrum"
    ))

    # Add peak annotations
    for ann in annotations:
        fig.add_annotation(
            x=ann["freq"],
            y=max(spectrum),
            text=ann["text"],
            showarrow=True,
            arrowhead=1,
            ax=0,
            ay=-30,
            textangle=-90,
            font=dict(color='red', size=10)
        )

    # Layout
    fig.update_layout(
        title="IR Spectrum",
        xaxis=dict(
            title="Wavenumber (cm⁻¹)",
            autorange="reversed"  # Important for IR
        ),
        yaxis=dict(
            title="Transmittance" if transmittance else "Intensity (a.u.)"
        ),
        showlegend=False
    )

    # Apply manual range if provided
    if plot_range:
        fig.update_xaxes(range=plot_range[::-1])  # reverse manually

    return fig

def gaussian(x, x0, fwhm: float):
    """Unit-height Gaussian centered at ``x0`` with the given full width at half max.

    ``x`` may be a scalar or numpy array (vectorized over the x-axis).
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-(x - x0)**2 / (2 * sigma**2))

def generate_absorption_emission_spectrum_interactive(
    wavelengths, oscs, points: int = 10000,
    plot_range: tuple[float, float] | None = None,
) -> go.Figure | None:
    """Build an interactive UV-Vis spectrum from transition wavelengths + oscillator strengths.

    Each transition is broadened with a Gaussian whose width grows with its oscillator
    strength, weighted by that strength, and summed onto a wavelength grid. Returns a
    normalized Plotly figure, or ``None`` when no transitions are supplied.
    """
    # Check if wavelengths and oscs are valid
    if wavelengths is None or len(wavelengths) == 0:
        print(f"No absorption peak detected.")
        return
    
    # Create the chemical shift axis
    if plot_range:
        min_wavelength, max_wavelength = plot_range
    else:
        min_wavelength = np.min(wavelengths) - 10
        max_wavelength = np.max(wavelengths) + 10
    
    x = np.linspace(min_wavelength, max_wavelength, points)
    
    # Initialize the spectrum
    spectrum = np.zeros_like(x)
    
    # Generate annotations for each peak with Gaussian broadening
    annotations = []
    base_fwhm = 10.0
    for wavelength, osc in zip(wavelengths, oscs):
        fwhm = base_fwhm * (1 + osc)   # stronger osc → broader peak
        spectrum += osc * gaussian(x, wavelength, fwhm)
        annotations.append({
            'wavelength': wavelength,
            'oscilation strength': oscs,
            'text': f"{wavelength:.2f}",
        })

    # Normalize the spectrum (guard against all-zero oscillator strengths)
    max_intensity = np.max(spectrum)
    if max_intensity != 0:
        spectrum /= max_intensity

    # Create the Plotly figure
    fig = go.Figure()

    # Add the spectrum trace
    fig.add_trace(go.Scatter(
        x=x,
        y=spectrum,
        mode='lines',
        line=dict(color='black'),
        name=f'Spectrum'
    ))
    
    # Add peak annotations
    for ann in annotations:
        fig.add_annotation(
            x=ann['wavelength'],
            y=1.05,
            text=ann['text'],
            showarrow=True,
            arrowhead=1,
            ax=0,
            ay=-40,
            yshift=10,
            textangle=-90,
            font=dict(color='red'),
        )

    fig.update_layout(
        xaxis=dict(
            title='Wavelength (nm)',
        ),
        yaxis=dict(
            title='Relative Intensity (a.u.)',
        ),
        title=f'Spectrum',
        showlegend=False,
    )

    # Adjust plot range if specified
    if plot_range:
        fig.update_xaxes(range=plot_range)

    return fig

def generate_ecd_spectrum_interactive(
    wavelengths, rotatory_strengths, points: int = 10000,
    plot_range: tuple[float, float] | None = None,
    fwhm_wavenumber: float = 3000.0,
) -> go.Figure | None:
    """Build an interactive ECD spectrum from wavelengths + rotatory strengths.

    Differs from :func:`generate_absorption_emission_spectrum_interactive` in two ways
    that matter physically:

    * Rotatory strengths are **signed**. They are summed with their sign so positive and
      negative Cotton effects appear above and below zero, and neighbouring bands of
      opposite sign partially cancel. Taking magnitudes here would destroy the
      information that makes ECD useful (two enantiomers differ only by the sign of
      every R, so their spectra are mirror images).
    * Broadening is done in **wavenumber space**, where a Gaussian band shape is the
      physically appropriate model, and only the display axis is wavelength. Broadening
      with a constant width directly in nm would make bands artificially narrow at short
      wavelengths and wide at long ones.

    ``fwhm_wavenumber`` is the band width in cm^-1 (~3000 cm^-1 is a conventional
    starting point for ECD). The curve is scaled by its largest absolute value, so the
    y-axis is relative and signed, in [-1, 1]. Returns ``None`` when no transitions are
    supplied.
    """
    if wavelengths is None or len(wavelengths) == 0:
        print("No ECD transitions to plot.")
        return None

    wavelengths = np.asarray(wavelengths, dtype=float)
    rotatory_strengths = np.asarray(rotatory_strengths, dtype=float)

    if plot_range:
        min_wavelength, max_wavelength = plot_range
    else:
        min_wavelength = np.min(wavelengths) - 50
        max_wavelength = np.max(wavelengths) + 50
    # The nm -> cm^-1 conversion below divides by the axis, so it must stay positive.
    min_wavelength = max(float(min_wavelength), 1.0)
    max_wavelength = max(float(max_wavelength), min_wavelength + 1.0)

    x = np.linspace(min_wavelength, max_wavelength, points)
    grid_wavenumbers = 1e7 / x

    spectrum = np.zeros_like(x)
    annotations = []
    for wavelength, rotatory_strength in zip(wavelengths, rotatory_strengths):
        if wavelength <= 0:
            continue
        spectrum += rotatory_strength * gaussian(grid_wavenumbers, 1e7 / wavelength, fwhm_wavenumber)
        annotations.append({
            'wavelength': wavelength,
            'rotatory_strength': rotatory_strength,
            'text': f"{wavelength:.2f}",
        })

    # Scale by the largest absolute value so the sign is preserved.
    max_intensity = np.max(np.abs(spectrum))
    if max_intensity != 0:
        spectrum /= max_intensity

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=x,
        y=spectrum,
        mode='lines',
        line=dict(color='black'),
        name='ECD',
    ))
    # ECD crosses zero, so the baseline is a meaningful reference and is drawn in.
    fig.add_hline(y=0, line_width=1, line_color='grey')

    # Annotate each transition on the side its own sign points to.
    for ann in annotations:
        above = ann['rotatory_strength'] >= 0
        fig.add_annotation(
            x=ann['wavelength'],
            y=1.05 if above else -1.05,
            text=ann['text'],
            showarrow=True,
            arrowhead=1,
            ax=0,
            ay=-40 if above else 40,
            yshift=10 if above else -10,
            textangle=-90,
            font=dict(color='red' if above else 'blue'),
        )

    fig.update_layout(
        xaxis=dict(title='Wavelength (nm)'),
        # Symmetric range so a mirror-image (opposite enantiomer) spectrum is directly
        # comparable, and so the zero crossing sits in the middle.
        yaxis=dict(title='Relative Δε (a.u.)', range=[-1.3, 1.3]),
        title='ECD Spectrum',
        showlegend=False,
    )

    if plot_range:
        fig.update_xaxes(range=[min_wavelength, max_wavelength])

    return fig

def parse_nmr_shielding_constants(logfile: str) -> list[tuple[int, str, float]]:
    shielding_data = []

    pattern = re.compile(
        r'^\s*(\d+)\s+([A-Za-z]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
    )

    with open(logfile, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if 'CHEMICAL SHIELDING SUMMARY' in line:
            j = i + 1

            # Skip header lines
            while j < len(lines) and not re.match(r'\s*\d+', lines[j]):
                j += 1

            # Read data
            while j < len(lines):
                line = lines[j]
                match = pattern.match(line)

                if match:
                    atom_idx = int(match.group(1))
                    element_symbol = match.group(2)
                    isotropic_value = float(match.group(3))

                    shielding_data.append(
                        (atom_idx, element_symbol, isotropic_value)
                    )
                else:
                    break  # stop when table ends

                j += 1

            break

    if not shielding_data:
        print("No NMR shielding data found in the log file.")

    return shielding_data

def parse_nmr_jcouplings(logfile: str) -> list[tuple[int, str, int, str, float]]:
    """Parse NMR spin-spin (J) coupling constants from an ORCA output file.

    Returns a list of tuples: (atom_a_idx, element_a, atom_b_idx, element_b, J_iso_Hz).
    """
    jcouplings = []

    # ORCA 6.x property block format, e.g.:
    #   NUCLEUS A = C    0 NUCLEUS B = C    1
    # (element symbol comes BEFORE the atom index for both A and B).
    pair_pattern = re.compile(
        r'NUCLEUS\s+A\s*=\s*([A-Za-z]+)\s+(\d+)\s+'
        r'NUCLEUS\s+B\s*=\s*([A-Za-z]+)\s+(\d+)',
        re.IGNORECASE
    )
    # The isotropic coupling is on the diagonalized total line, e.g.:
    #   J[0,1](Total)   37.220   46.839   48.366  iso=   44.142
    # Anchor on "J[" so the reduced-coupling "K[..](Total)" block is ignored.
    iso_pattern = re.compile(
        r'J\[\d+,\d+\]\(Total\).*?iso\s*=\s*(-?\d+\.\d+)',
        re.IGNORECASE
    )

    with open(logfile, 'r') as f:
        lines = f.readlines()

    current_pair = None
    for line in lines:
        m = pair_pattern.search(line)
        if m:
            # Store as (idx_a, elem_a, idx_b, elem_b) to match the return shape.
            current_pair = (int(m.group(2)), m.group(1), int(m.group(4)), m.group(3))
            continue
        if current_pair is not None:
            m_iso = iso_pattern.search(line)
            if m_iso:
                jcouplings.append((*current_pair, float(m_iso.group(1))))
                current_pair = None

    if not jcouplings:
        print("No NMR spin-spin coupling data found in the log file.")
    return jcouplings

def calculate_chemical_shifts(shielding_data: list[tuple[int, str, float]],
                              reference_shieldings: dict[str, float]) -> list[tuple[int, str, float]]:
    """Convert absolute shieldings to chemical shifts (delta = reference - shielding).

    Atoms whose element has no entry in ``reference_shieldings`` are skipped.
    """
    shifts_data = []
    for atom_idx, element, shielding in shielding_data:
        if element in reference_shieldings:
            delta = reference_shieldings[element] - shielding
            shifts_data.append((atom_idx, element, delta))
        else:
            print(f"No reference shielding constant for element {element}. Skipping atom {atom_idx}.")
    return shifts_data

def lorentzian_nmr(x, x0, gamma):
    """Unit-height Lorentzian line centered at ``x0`` with half-width ``gamma`` (ppm).

    ``x`` may be a scalar or numpy array (vectorized over the chemical-shift axis).
    """
    return (gamma**2) / ((x - x0)**2 + gamma**2)

# Cordero covalent radii (Angstrom) for the elements common in organic NMR.
# Used only for simple distance-based bond perception.
_COVALENT_RADII = {
    'H': 0.31, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39,
}

def parse_orca_geometry(logfile: str) -> tuple[list[str], np.ndarray]:
    """Parse the last 'CARTESIAN COORDINATES (ANGSTROEM)' block from an ORCA log.

    Returns (symbols, coords) in the SAME atom order as ORCA's NMR nucleus
    indices (which matches this block but NOT cclib's reordered atom list).
    coords is an (N, 3) numpy array in Angstrom.
    """
    line_pattern = re.compile(
        r'^\s*([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$'
    )
    with open(logfile, 'r') as f:
        lines = f.readlines()

    symbols, coords = [], []
    for i, line in enumerate(lines):
        if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
            block_symbols, block_coords = [], []
            j = i + 1
            while j < len(lines) and set(lines[j].strip()) <= set('-'):
                j += 1  # skip the dashed separator line(s)
            while j < len(lines):
                m = line_pattern.match(lines[j])
                if not m:
                    break
                block_symbols.append(m.group(1))
                block_coords.append([float(m.group(2)), float(m.group(3)), float(m.group(4))])
                j += 1
            if block_symbols:
                symbols, coords = block_symbols, block_coords  # keep the last block
    return symbols, np.asarray(coords, dtype=float).reshape(-1, 3)

def perceive_bonds(coords, symbols: list[str], tolerance: float = 1.3) -> list[set[int]]:
    """Infer bonds from interatomic distances and covalent radii.

    Returns a list of sets, where adjacency[i] holds the indices bonded to atom i.
    Two atoms are bonded when their distance is below (r_i + r_j) * tolerance.
    """
    coords = np.asarray(coords, dtype=float)
    n = len(symbols)
    adjacency = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            r_i = _COVALENT_RADII.get(symbols[i], 0.77)
            r_j = _COVALENT_RADII.get(symbols[j], 0.77)
            dist = float(np.linalg.norm(coords[i] - coords[j]))
            if 0.0 < dist < (r_i + r_j) * tolerance:
                adjacency[i].add(j)
                adjacency[j].add(i)
    return adjacency

def morgan_ranks(coords, symbols: list[str]) -> tuple[list[set[int]], list[int]]:
    """Assign each atom a topological-symmetry rank (Morgan-style refinement).

    Atoms start with an invariant of (element, degree) on the perceived bond graph
    and are iteratively re-ranked by (own rank, sorted neighbour ranks) until the set
    of distinct ranks stops growing. Two atoms sharing a final rank are equivalent by
    the symmetry of the whole molecular graph. Ranks are opaque grouping keys only.
    """
    adjacency = perceive_bonds(coords, symbols)
    n = len(symbols)
    ranks = [hash((symbols[i], len(adjacency[i]))) for i in range(n)]
    for _ in range(n):  # converges in at most N refinement rounds
        new_ranks = [
            hash((ranks[i], tuple(sorted(ranks[j] for j in adjacency[i]))))
            for i in range(n)
        ]
        if len(set(new_ranks)) == len(set(ranks)):
            break  # partition stable: no new classes distinguished
        ranks = new_ranks
    return adjacency, ranks

def compute_hydrogen_equivalence(coords, symbols: list[str]) -> dict[int, object]:
    """Group hydrogens that are bonded to symmetry-equivalent heavy atoms.

    Freely-rotating groups (CH3, CH2, NH2) average their protons to a single signal,
    and so do protons on heavy atoms that are equivalent by molecular symmetry (e.g.
    the two methyls of isopropanol or tert-butyl). Hydrogens are therefore keyed on the
    Morgan rank of their bonded heavy atom, which reduces to the bonded-atom identity
    when the molecule has no symmetry. ``coords``/``symbols`` must be in ORCA NMR
    nucleus order (use parse_orca_geometry, NOT cclib, which reorders atoms). Returns a
    dict {hydrogen_atom_idx: group_key}.

    NOTE: equivalence is topological only. Diastereotopic protons on the same carbon are
    treated as equivalent (a deliberate simplification, since a static geometry carries
    no stereochemistry in the graph).
    """
    adjacency, ranks = morgan_ranks(coords, symbols)
    equivalence = {}
    for idx, sym in enumerate(symbols):
        if sym != 'H':
            continue
        heavy_neighbors = [j for j in adjacency[idx] if symbols[j] != 'H']
        # Key on the bonded heavy atom's symmetry rank so protons on equivalent heavy
        # atoms merge. With no heavy neighbor (e.g. H2), the hydrogen forms its own group.
        equivalence[idx] = ('heavy', ranks[heavy_neighbors[0]]) if heavy_neighbors else f"H{idx}"
    return equivalence

def compute_carbon_equivalence(coords, symbols: list[str]) -> dict[int, int]:
    """Group carbons that are topologically equivalent by molecular symmetry.

    Unlike hydrogens (which also average within a freely-rotating group), equivalent
    carbons are related by the symmetry of the whole molecule (e.g. the two methyls of
    isopropanol, the three methyls of tert-butyl, or the para carbons of a symmetric
    ring). Carbons sharing a Morgan rank (see morgan_ranks) merge into a single averaged
    signal, independent of small DFT shift differences.

    ``coords``/``symbols`` must be in ORCA NMR nucleus order (use parse_orca_geometry,
    NOT cclib, which reorders atoms). Returns a dict {carbon_atom_idx: rank}; the rank
    is only used as an opaque grouping key.
    """
    _, ranks = morgan_ranks(coords, symbols)
    return {i: ranks[i] for i, sym in enumerate(symbols) if sym == 'C'}

# Multiplicity letters for first-order multiplets (number of equal coupling
# partners -> name). More than one coupling group concatenates letters (e.g. dd).
_MULTIPLICITY_LETTER = {1: 'd', 2: 't', 3: 'q', 4: 'p', 5: 'h', 6: 'hept'}

def build_nmr_peak_table(shifts_data: list[tuple[int, str, float]], element_symbol: str,
                         jcouplings: list[tuple[int, str, int, str, float]] | None = None,
                         equivalence: dict[int, object] | None = None,
                         j_threshold: float = 0.5) -> list[dict]:
    """Group nuclei of ``element_symbol`` into averaged multiplet peaks.

    Equivalent nuclei merge into one peak whose chemical shift is the average of
    the members and whose intensity is the member count. ``equivalence`` is an
    optional {atom_idx: group_key} map (e.g. from compute_hydrogen_equivalence);
    when omitted, nuclei are grouped by rounded chemical shift (~0.01 ppm).

    Coupling between two equivalent groups A and B is treated to first order:
    every nucleus of A is split by B into (size(B) + 1) lines with a single
    coupling constant equal to the AVERAGE of all pairwise A-B couplings. The
    multiplicity therefore comes from the partner group SIZE (not from how many
    couplings the calculation happened to print for one atom), which is what
    makes a CH3-CH2 pair read as the expected triplet/quartet even though the
    individual per-atom couplings differ in a single static geometry. Couplings
    whose averaged magnitude is below ``j_threshold`` (Hz) are dropped as
    unresolvable.

    Returns a list of peak dicts sorted by ascending shift::

        {'atom_idxs': [...], 'element': str, 'shift': float, 'count': int,
         'couplings': [{'partner_idxs': (...), 'n': int, 'J': float}, ...]}
    """
    shifts = [(i, e, s) for i, e, s in shifts_data if e == element_symbol]
    if not shifts:
        return []

    def key_of(atom_idx, shift):
        if equivalence is not None and atom_idx in equivalence:
            return ('grp', equivalence[atom_idx])
        return ('shift', round(shift, 2))

    groups = {}
    for atom_idx, element, shift in shifts:
        k = key_of(atom_idx, shift)
        grp = groups.setdefault(k, {'atom_idxs': [], 'shifts': []})
        grp['atom_idxs'].append(atom_idx)
        grp['shifts'].append(shift)

    # Pairwise homonuclear coupling lookup: (atom_a, atom_b) -> J_Hz (both ways).
    pair_j = {}
    if jcouplings:
        for idx_a, elem_a, idx_b, elem_b, j_hz in jcouplings:
            if elem_a == element_symbol and elem_b == element_symbol:
                pair_j[(idx_a, idx_b)] = j_hz
                pair_j[(idx_b, idx_a)] = j_hz

    group_items = list(groups.items())
    peaks = []
    for k, grp in group_items:
        atom_idxs = sorted(grp['atom_idxs'])
        avg_shift = sum(grp['shifts']) / len(grp['shifts'])
        couplings = []
        for k2, grp2 in group_items:
            if k2 == k:
                continue  # mutually equivalent nuclei do not split each other
            partner_idxs = sorted(grp2['atom_idxs'])
            js = [pair_j[(a, b)] for a in atom_idxs for b in partner_idxs if (a, b) in pair_j]
            if not js:
                continue
            j_avg = sum(js) / len(js)
            if abs(j_avg) < j_threshold:
                continue  # negligible / unresolvable coupling
            couplings.append({
                'partner_idxs': tuple(partner_idxs),
                'n': len(partner_idxs),  # multiplicity from partner GROUP size
                'J': j_avg,
            })
        couplings.sort(key=lambda c: -abs(c['J']))
        peaks.append({
            'atom_idxs': atom_idxs,
            'element': element_symbol,
            'shift': avg_shift,
            'count': len(atom_idxs),
            'couplings': couplings,
        })

    peaks.sort(key=lambda p: p['shift'])
    return peaks

def format_jcouplings_by_atom(jcouplings: list[tuple[int, str, int, str, float]] | None,
                              decimals: int = 4) -> dict[int, str]:
    """Summarise each nucleus's J-couplings as a readable string, keyed by atom index.

    Every coupling is listed from both partners' point of view, so a single ``H0-H1``
    entry in the log produces text for H0 naming H1 *and* text for H1 naming H0::

        {0: '1.0000 Hz (H1)', 1: '1.0000 Hz (H0)'}

    Partners are written as element symbol + atom index using ORCA's own 0-based nucleus
    numbering — the same numbering the structure viewer labels atoms with — so a value
    read off this table can be matched to an atom on screen. Couplings are listed
    strongest first; nuclei with no coupling are absent from the mapping.
    """
    by_atom: dict[int, list[tuple[float, str]]] = {}
    seen_pairs = set()

    for idx_a, elem_a, idx_b, elem_b, j_hz in jcouplings or []:
        # A log that prints a pair from both directions must not double-list it.
        pair = (min(idx_a, idx_b), max(idx_a, idx_b))
        if idx_a == idx_b or pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        by_atom.setdefault(idx_a, []).append((j_hz, f'{elem_b}{idx_b}'))
        by_atom.setdefault(idx_b, []).append((j_hz, f'{elem_a}{idx_a}'))

    formatted = {}
    for atom_idx, entries in by_atom.items():
        entries.sort(key=lambda entry: -abs(entry[0]))
        formatted[atom_idx] = ', '.join(
            f'{j_hz:.{decimals}f} Hz ({partner})' for j_hz, partner in entries)
    return formatted

def multiplicity_label(couplings: list[dict]) -> str:
    """Human-readable multiplicity (s/d/t/q/...) for a peak's coupling list."""
    if not couplings:
        return 's'
    return ''.join(_MULTIPLICITY_LETTER.get(c['n'], 'm') for c in couplings)

def generate_nmr_spectrum_interactive(
    shifts_data: list[tuple[int, str, float]], element_symbol: str,
    linewidth: float = 0.5, points: int = 10000, frequency: float = 400,
    plot_range: tuple[float, float] | None = None,
    jcouplings: list[tuple[int, str, int, str, float]] | None = None,
    equivalence: dict[int, object] | None = None,
) -> go.Figure | None:
    """Build an interactive NMR spectrum for ``element_symbol`` from chemical shifts.

    Nuclei are grouped into averaged multiplet peaks (:func:`build_nmr_peak_table`),
    each rendered as a first-order Pascal multiplet: coupling to ``n`` equivalent
    partners splits a line into ``n + 1`` sub-peaks separated by ``J / frequency``
    (Hz -> ppm), broadened with a Lorentzian (:func:`lorentzian_nmr`). The x-axis is
    reversed per NMR convention. Returns a Plotly figure, or ``None`` when no nuclei
    of ``element_symbol`` are present.
    """
    # Filter shifts for the desired element
    shifts = [(atom_idx, element, shift) for atom_idx, element, shift in shifts_data if element == element_symbol]

    if not shifts:
        print(f"No chemical shifts found for element {element_symbol}.")
        return

    # Group nuclei into averaged multiplet peaks. Equivalent nuclei (per the
    # `equivalence` map, or by rounded shift when absent) collapse to one peak
    # whose shift is averaged and whose height scales with the member count.
    peaks = build_nmr_peak_table(shifts_data, element_symbol,
                                 jcouplings=jcouplings, equivalence=equivalence)

    # Create the chemical shift axis
    if plot_range:
        min_shift, max_shift = plot_range
    else:
        min_shift = min(shift for _, _, shift in shifts) - 1
        max_shift = max(shift for _, _, shift in shifts) + 1

    x = np.linspace(min_shift, max_shift, points)

    # Initialize the spectrum
    spectrum = np.zeros_like(x)

    # Convert linewidth from Hz to ppm
    lw_ppm = linewidth / frequency

    # Generate Lorentzian multiplets per group and collect positions for annotations
    annotations = []
    for peak in peaks:
        shift = peak['shift']
        intensity = peak['count']
        # Use the n+1 distinct Pascal lines directly. Repeated binary splitting
        # materializes 2**n entries, most at duplicate positions, and can freeze
        # the result tab for large equivalent partner groups.
        sub_peaks = [(shift, 1.0)]
        for c in peak['couplings']:
            n = c['n']
            spacing = c['J'] / frequency
            sub_peaks = [
                (p + (k - n / 2) * spacing, w * comb(n, k) / 2**n)
                for p, w in sub_peaks for k in range(n + 1)
            ]
        for pos, weight in sub_peaks:
            spectrum += intensity * weight * lorentzian_nmr(x, pos, lw_ppm)
        label = f"{element_symbol}-" + ",".join(str(i) for i in peak['atom_idxs'])
        annotations.append({
            'shift': shift,
            'atom_label': label,
            'intensity': intensity
        })

    # Normalize the spectrum
    spectrum /= np.max(spectrum)

    # Create the Plotly figure
    fig = go.Figure()

    # Add the spectrum trace
    fig.add_trace(go.Scatter(
        x=x,
        y=spectrum,
        mode='lines',
        line=dict(color='black'),
        name=f'{element_symbol} NMR Spectrum'
    ))

    # Annotate at the centroid (chemical shift) using the actual spectrum value
    # so the arrow lands on the multiplet envelope rather than above it.
    for ann in annotations:
        idx = int(np.argmin(np.abs(x - ann['shift'])))
        fig.add_annotation(
            x=ann['shift'],
            y=spectrum[idx],
            text=ann['atom_label'],
            showarrow=True,
            arrowhead=1,
            ax=0,
            ay=-40,
            yshift=10,
            font=dict(color='red'),
        )

    # Invert x-axis for NMR spectrum convention
    fig.update_layout(
        xaxis=dict(
            autorange='reversed',
            title='Chemical Shift (ppm)',
        ),
        yaxis=dict(
            title='Intensity (a.u.)',
        ),
        title=f'{element_symbol} NMR Spectrum',
        showlegend=False,
    )

    # Adjust plot range if specified
    if plot_range:
        fig.update_xaxes(range=sorted(plot_range, reverse=True))

    return fig
