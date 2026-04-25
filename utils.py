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
        nmr_keywords = 'NMR'
        if spin_spin_coupling:
            nmr_keywords += ' SpinSpin'
        f.write(f'! {method} {nmr_keywords}\n\n')
        get_orca_header(f, n_proc, memory, solvation, solvation_model, solvent)
        get_orca_coordinates(f, mol, charge, multiplicity)

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

def generate_nmr_spectrum_interactive(shifts_data, element_symbol, linewidth=0.5, points=10000, frequency=400, plot_range=None):
    # Filter shifts for the desired element
    shifts = [(atom_idx, element, shift) for atom_idx, element, shift in shifts_data if element == element_symbol]

    if not shifts:
        print(f"No chemical shifts found for element {element_symbol}.")
        return

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

    # Count occurrences of each chemical shift for intensity scaling
    shift_counts = {}
    for atom_idx, element, shift in shifts:
        shift_counts[shift] = shift_counts.get(shift, 0) + 1

    # Generate Lorentzian peaks with intensities and collect peak positions for annotations
    annotations = []
    for atom_idx, element, shift in shifts:
        intensity = shift_counts[shift]
        spectrum += intensity * lorentzian_nmr(x, shift, lw_ppm)
        annotations.append({
            'shift': shift,
            'atom_label': f"{element}-{atom_idx}",
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

    # Add peak annotations
    for ann in annotations:
        fig.add_annotation(
            x=ann['shift'],
            y=lorentzian_nmr(ann['shift'], ann['shift'], lw_ppm) / np.max(spectrum),
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