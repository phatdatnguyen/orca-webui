import os
import re
import cclib
from rdkit import Chem
import numpy as np
import plotly.graph_objects as go

def get_files_in_working_directory(working_directory_path):
    files = [f for f in os.listdir(working_directory_path) if not f.endswith('Zone.Identifier')]
    return files

def conformer_to_xyz_file(mol, conf_id, file_path, charge=0, multiplicity=1):
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

def add_bonds(mol, bond_factor = 1.25):
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

    # Convert to a regular Mol object
    mol_new = mol_new.GetMol()
    
    return mol_new

def mol_from_xyz_file(file_path, return_charge_and_multiplicity=False):
    with open(file_path, 'r') as file:
        xyz_string = file.read()

    lines = xyz_string.split('\n')
    if len(lines[0].split())==2:
        charge, multiplicity = map(int, lines[0].split())
    else:
        charge, multiplicity = 0, 1

    # Create a new empty molecule
    mol = Chem.RWMol()

    # Parse the atomic coordinates and add atoms
    elements = []
    for line in lines[1:]:
        if len(line.split())==4:
            parts = line.split()
            element = parts[0].capitalize()
            atom = Chem.Atom(element)
            mol.AddAtom(atom)
            elements.append(element)
        else:
            continue

    # Add 3D coordinates to the molecule
    conf = Chem.Conformer(mol.GetNumAtoms())
    atom_idx = 0
    for line in lines[1:]:
        if len(line.split())==4:
            parts = line.split()
            x, y, z = map(float, parts[1:4])
            conf.SetAtomPosition(atom_idx, (x, y, z))
            atom_idx += 1
        else:
            continue
    mol.AddConformer(conf)
                
    # Convert to a regular Mol object
    mol = mol.GetMol()
    
    if return_charge_and_multiplicity:
        return mol, charge, multiplicity
    else:
        return mol

def mol_from_orca_file(filename):
    # Parse the Gaussian file using cclib
    parser = cclib.io.ccopen(filename)
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

def get_orca_method_keyword(method_type, method_name, functional, basis):
    if method_type == "HF":
        return f'HF {basis}'
    elif method_type == "DFT":
        return f'{functional} {basis}'
    elif method_type in ["MP2", "CCSD", "CCSD(T)"]:
        return f'{method_type} {basis}'
    else:  # method_type == "Semi-empirical":
        return f'{method_name} {basis}'

def get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent, tddft=False, n_states=10, maxdim=5, root=0):
    f.write(f'%pal\n  nprocs {n_proc}\nend\n')
    if memory is not None:
        maxcore = int(memory * 1024)
        f.write(f'%maxcore {maxcore}\n')
    if tddft:
        f.write('%tddft\n')
        f.write(f'  nroots {n_states}\n')
        f.write(f'  maxdim {maxdim}\n')
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

def get_orca_coordinates(f, mol, charge, multiplicity):
    f.write(f'* xyz {charge} {multiplicity}\n')
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        pos = conf.GetAtomPosition(idx)
        f.write(f'{symbol} {pos.x:>12.6f} {pos.y:>12.6f} {pos.z:>12.6f}\n')
    f.write('*\n')

def write_sp_orca_input(mol, file_name, method_type, method_name, functional='B3LYP', basis='6-31G(d)', charge=0, multiplicity=1, solvation=False, solvation_model=None, solvent=None, n_proc=4, memory=2):
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, method_name, functional, basis)
        f.write(f'! {method} SP\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent)
        get_orca_coordinates(f, mol, charge, multiplicity)

def write_opt_orca_input(mol, file_name, method_type, method_name, functional='B3LYP', basis='6-31G(d,p)', charge=0, multiplicity=1, solvation=False, solvation_model=None, solvent=None, n_proc=4, memory=2):
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, method_name, functional, basis)
        f.write(f'! {method} Opt\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent)
        get_orca_coordinates(f, mol, charge, multiplicity)

def write_opt_freq_orca_input(mol, file_name, method_type, method_name, functional='B3LYP', basis='6-31G(d,p)', charge=0, multiplicity=1, solvation=False, solvation_model=None, solvent=None, n_proc=4, memory=2):
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, method_name, functional, basis)
        f.write(f'! {method} Opt Freq\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent)
        get_orca_coordinates(f, mol, charge, multiplicity)

def write_tddft_orca_input(mol, file_name, method_type, method_name, functional='B3LYP', basis='6-31G(d,p)', n_states=10, maxdim=5, root=0, charge=0, multiplicity=1, solvation=False, solvation_model=None, solvent=None, n_proc=4, memory=2):
    with open(file_name + '.inp', 'w') as f:
        method = get_orca_method_keyword(method_type, method_name, functional, basis)
        f.write(f'! {method} RIJCOSX TightSCF Opt\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent, tddft=True, n_states=n_states, maxdim=maxdim, root=root)
        get_orca_coordinates(f, mol, charge, multiplicity)

def write_nmr_orca_input(mol, file_name, method_type, functional='B3LYP', basis='6-31G(d,p)', spin_spin_coupling=False, charge=0, multiplicity=1, solvation=False, solvation_model=None, solvent=None, n_proc=4, memory=2):
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
        

def lorentzian_ir(wavenumber, position, intensity, width=10):
    return intensity * (width ** 2) / ((wavenumber - position) ** 2 + width ** 2)

def generate_ir_spectrum_interactive(
    frequencies,
    intensities,
    width=10.0,
    points=4000,
    plot_range=None,
    normalize=True,
    transmittance=True
):
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

def gaussian(x, x0, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-(x - x0)**2 / (2 * sigma**2))

def generate_absorption_emission_spectrum_interactive(wavelengths, oscs, points=10000, plot_range=None):
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

def parse_nmr_shielding_constants(logfile):
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

def parse_nmr_jcouplings(logfile):
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

def calculate_chemical_shifts(shielding_data, reference_shieldings):
    shifts_data = []
    for atom_idx, element, shielding in shielding_data:
        if element in reference_shieldings:
            delta = reference_shieldings[element] - shielding
            shifts_data.append((atom_idx, element, delta))
        else:
            print(f"No reference shielding constant for element {element}. Skipping atom {atom_idx}.")
    return shifts_data

def lorentzian_nmr(x, x0, gamma):
    return (gamma**2) / ((x - x0)**2 + gamma**2)

# Cordero covalent radii (Angstrom) for the elements common in organic NMR.
# Used only for simple distance-based bond perception.
_COVALENT_RADII = {
    'H': 0.31, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39,
}

def parse_orca_geometry(logfile):
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
    return symbols, np.asarray(coords, dtype=float)

def perceive_bonds(coords, symbols, tolerance=1.3):
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

def morgan_ranks(coords, symbols):
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

def compute_hydrogen_equivalence(coords, symbols):
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

def compute_carbon_equivalence(coords, symbols):
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

def build_nmr_peak_table(shifts_data, element_symbol, jcouplings=None, equivalence=None, j_threshold=0.5):
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

def multiplicity_label(couplings):
    """Human-readable multiplicity (s/d/t/q/...) for a peak's coupling list."""
    if not couplings:
        return 's'
    return ''.join(_MULTIPLICITY_LETTER.get(c['n'], 'm') for c in couplings)

def generate_nmr_spectrum_interactive(shifts_data, element_symbol, linewidth=0.5, points=10000, frequency=400, plot_range=None, jcouplings=None, equivalence=None):
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
        # First-order multiplet: each coupled partner halves the weight and
        # splits the peak by J/freq (Hz -> ppm). Coupling to n equivalent
        # partners splits n times, giving an (n+1)-line Pascal multiplet.
        sub_peaks = [(shift, 1.0)]
        for c in peak['couplings']:
            half_j_ppm = (c['J'] / 2.0) / frequency
            for _ in range(c['n']):
                sub_peaks = [(p - half_j_ppm, w / 2.0) for p, w in sub_peaks] + \
                            [(p + half_j_ppm, w / 2.0) for p, w in sub_peaks]
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
        fig.update_xaxes(range=plot_range)

    return fig