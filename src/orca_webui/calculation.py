"""Calculation tab: turn a structure file into an ORCA input file and run ORCA.

The upper section builds an ``.inp`` file for the chosen calculation type (single
point, optimization, frequency, TD-DFT, NMR) and method/basis/solvation settings;
the lower section invokes the external ``orca`` binary on a selected ``.inp`` and
captures its stdout to a ``.log``.
"""
import os
import time
import multiprocessing
import psutil
import shutil
import subprocess
import threading
import math
import gradio as gr
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import cclib
from .utils import *

# Physical RAM in whole GB.
max_memory = math.floor(psutil.virtual_memory().total / (1024 ** 3))
# Upper bound / default of the core-count slider. ORCA fails to start when asked for
# more than half the logical CPU count (that half is the physical-core count on an
# SMT/hyperthreaded machine), so the logical count is halved here.
max_cores = max(1, multiprocessing.cpu_count() // 2)

def max_memory_per_core(n_cores):
    """Largest safe ``%maxcore`` value in whole GB for ``n_cores`` parallel processes.

    ORCA's ``%maxcore`` is the memory budget of *each* process, not of the calculation
    as a whole, so ORCA may claim ``nprocs x maxcore`` in total. Offering the full
    physical RAM per core lets a parallel run exhaust the machine, so the per-core
    ceiling is physical RAM divided by the core count, rounded down.
    """
    return max(1, math.floor(max_memory / max(1, int(n_cores))))
# Functional and basis-set menus. These live at module level because
# on_change_calculation_type rebuilds both dropdowns per calculation type: keeping one
# copy of each list stops the rebuilt menus from drifting from the initial ones.
FUNCTIONALS = ["BLYP", "B3LYP", "CAM-B3LYP", "B3PW91", "PBE", "PBE0", "revPBE", "B97-3C",
               "M06L", "M062X", "wB97X", "wB97X-D3", "wB97X-D3BJ", "wB97X-D4"]
# TD-DFT wants range-separated / GGA functionals; the hybrids above are dropped.
TDDFT_FUNCTIONALS = ["CAM-B3LYP", "PBE", "PBE0", "revPBE", "B97-3C",
                     "wB97X", "wB97X-D3", "wB97X-D3BJ", "wB97X-D4"]
_BASIS_SETS_GENERAL = ["STO-3G", "3-21G", "6-31G", "6-31G(d,p)", "6-31+G(d,p)", "6-31++G(d,p)",
                       "6-311G", "6-311G(d,p)", "cc-pVDZ", "cc-pVTZ", "cc-pVQZ",
                       "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]
# Jensen's segmented-contracted sets are built for NMR shielding, so they are offered
# only for NMR calculations.
_BASIS_SETS_NMR = ["pcSseg-0", "pcSseg-1", "pcSseg-2", "pcSseg-3"]
_BASIS_SETS_ECP = ["LanL2DZ", "SDD"]
BASIS_SETS = _BASIS_SETS_GENERAL + _BASIS_SETS_ECP
NMR_BASIS_SETS = _BASIS_SETS_GENERAL + _BASIS_SETS_NMR + _BASIS_SETS_ECP
DEFAULT_BASIS_SET = "3-21G"
# ORCA refuses to start a parallel run unless it is invoked by full pathname, so the
# binary is looked up on PATH at import time and offered as an editable default in the
# Settings accordion (empty when ORCA is not installed / not on PATH).
default_orca_path = shutil.which("orca") or ""

# Handle on the ORCA process launched by the Run button so the Stop button can kill it.
# The two clicks are served on different threads, hence the lock; _orca_stopped_by_user
# lets the run handler tell "user pressed Stop" apart from "ORCA failed".
_orca_process = None
_orca_stopped_by_user = False
_orca_lock = threading.Lock()

def on_working_directory_file_list_change(working_directory_file_list, input_file_name):
    """Repopulate the input-structure and input-file dropdowns from the file list.

    Keeps the currently-typed input file name selected if a matching ``.inp`` exists,
    otherwise falls back to the first available input file. Returns updates for the
    structure dropdown and the input-file dropdown. The file-list state is ``None``
    until a working directory is opened, so it is normalised to an empty list here.

    Both lists are sorted alphabetically (case-insensitively) rather than left in the
    arbitrary order ``os.listdir`` returns.
    """
    working_directory_file_list = working_directory_file_list or []
    input_file_name = input_file_name or ""
    structure_file_names = sorted((f for f in working_directory_file_list if f.endswith('.xyz') or f.endswith('.pdb') or f.endswith('.mol') or f.endswith('.mol2') or f.endswith('.log')), key=str.lower)
    input_file_names = sorted((f for f in working_directory_file_list if f.endswith('.inp')), key=str.lower)
    if input_file_name + ".inp" in input_file_names:
        input_file_name_value = input_file_name + ".inp"
    else:
        input_file_name_value = input_file_names[0] if len(input_file_names) > 0 else None
    
    return gr.update(choices=structure_file_names, value=structure_file_names[0] if len(structure_file_names) > 0 else None, interactive=True), \
           gr.update(choices=input_file_names, value=input_file_name_value, interactive=True)
   
def on_change_calculation_type(calculation_type, basis_set):
    """React to the calculation-type radio changing.

    Suggests a default input file name, restricts the available method types (NMR
    excludes semi-empirical/CC), swaps the functional list for TD-DFT, offers the
    NMR-only pcSseg basis sets for NMR, and toggles visibility of the TD-DFT sliders
    and the NMR spin-spin-coupling checkbox.
    """
    if calculation_type == "Single-Point":
        file_name = "single_point"
    elif calculation_type == "Geometry Optimization":
        file_name = "geometry_optimization"
    elif calculation_type == "Frequency":
        file_name = "frequency"
    elif calculation_type == "Time-Dependent Density Functional Theory":
        file_name = "tddft"
    else: # calculation_type == "NMR Spectrum"
        file_name = "nmr_spectrum"
    
    if calculation_type == "NMR Spectrum":
        method_types = ["HF", "DFT", "MP2"]
    else:
        method_types = ["HF", "DFT", "Semi-empirical", "MP2", "CCSD", "CCSD(T)"]

    if calculation_type == "Time-Dependent Density Functional Theory":
        functional_dropdown = gr.Dropdown(label="Functional", value="PBE", choices=TDDFT_FUNCTIONALS, allow_custom_value=True)
    else:
        functional_dropdown = gr.Dropdown(label="Functional", value="B3LYP", choices=FUNCTIONALS, allow_custom_value=True)

    basis_sets = NMR_BASIS_SETS if calculation_type == "NMR Spectrum" else BASIS_SETS
    # The dropdown allows custom values, so a pcSseg set left selected when switching
    # away from NMR would silently survive into a non-NMR input file. Fall back to the
    # default whenever the current choice is no longer on offer.
    basis_set_value = basis_set if basis_set in basis_sets else DEFAULT_BASIS_SET

    show_tddft = calculation_type == "Time-Dependent Density Functional Theory"

    show_spin_spin_couplings = (calculation_type == "NMR Spectrum")

    return file_name, gr.update(choices=method_types, value="DFT"), functional_dropdown, gr.update(choices=basis_sets, value=basis_set_value), gr.update(visible=show_tddft), gr.update(visible=show_tddft), gr.update(visible=show_tddft), gr.update(visible=show_spin_spin_couplings)

def on_mm_checkbox_change(mm_checkbox: gr.Checkbox):
    """Show/hide the force-field and max-iterations controls with the MM checkbox."""
    return gr.update(visible=mm_checkbox), gr.update(visible=mm_checkbox)

def on_solvation_checkbox_change(solvation_checkbox: gr.Checkbox):
    """Show/hide the solvation-model and solvent controls with the solvation checkbox."""
    return gr.update(visible=solvation_checkbox), gr.update(visible=solvation_checkbox)

def on_method_type_change(method_type):
    """Toggle method-name/functional/basis controls to match the chosen method type.

    Semi-empirical shows a method-name dropdown (AM1/PM3) and hides the basis set;
    DFT shows the functional; everything else shows only the basis set.
    """
    show_method_name = (method_type=="Semi-empirical")
    show_functional = method_type=="DFT"
    show_basis_set = not (method_type=="Semi-empirical")

    if method_type == "Semi-empirical":
        method_names = ["AM1", "PM3"]
        default_method_name = "AM1"
    else:
        method_names = [""]
        default_method_name = ""

    return gr.update(choices=method_names, value=default_method_name, visible=show_method_name), gr.update(visible=show_functional), gr.update(visible=show_basis_set)

def on_generate_input_file(working_directory_path, input_structure_file_dropdown, calculation_type,
                           use_mm, force_field, max_iters,
                           use_solvation, solvation_model, solvent,
                           method_type, method_name, functional, basis_set_textbox, n_states, maxdim, root, spin_spin_coupling, charge, multiplicity,
                           n_cores, memory, input_file_name):
    """Generate an ORCA ``.inp`` file for the selected structure and settings.

    Loads the structure (optionally pre-optimizing it with MMFF/UFF molecular
    mechanics), then dispatches to the matching ``write_*_orca_input`` helper based on
    ``calculation_type``. Returns ``(status_html, file_list)``; errors are reported in
    the status span rather than raised.
    """
    if input_structure_file_dropdown is None or input_structure_file_dropdown == "":
        gr.Warning("Please select an input structure")
        return "", get_files_in_working_directory(working_directory_path)
    
    try:
        # Get the molecule object. Same loader as the structure viewer, so the atom
        # order written below is the order the viewer labelled on screen.
        file_path = os.path.join(working_directory_path, input_structure_file_dropdown)
        mol = mol_from_structure_file(file_path)

        # Optimize geometry with molecular mechanics
        if use_mm:
            if force_field=="MMFF":
                AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
            else: # force_field=="UFF"
                AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)

        # Generate input file
        input_file_path = os.path.join(working_directory_path, input_file_name)
        if calculation_type == "Single-Point":
            write_sp_orca_input(mol, input_file_path, method_type, method_name, functional=functional, basis=basis_set_textbox,
                                charge=charge, multiplicity=multiplicity,
                                solvation=use_solvation, solvation_model=solvation_model, solvent=solvent, n_proc=n_cores, memory=memory)
        elif calculation_type == "Geometry Optimization":
            write_opt_orca_input(mol, input_file_path, method_type, method_name, functional=functional, basis=basis_set_textbox,
                                 charge=charge, multiplicity=multiplicity,
                                 solvation=use_solvation, solvation_model=solvation_model, solvent=solvent, n_proc=n_cores, memory=memory)
        elif calculation_type == "Frequency":
            write_opt_freq_orca_input(mol, input_file_path, method_type, method_name, functional=functional, basis=basis_set_textbox,
                                      charge=charge, multiplicity=multiplicity,
                                      solvation=use_solvation, solvation_model=solvation_model, solvent=solvent, n_proc=n_cores, memory=memory)
        elif calculation_type == "Time-Dependent Density Functional Theory":
            write_tddft_orca_input(mol, input_file_path, method_type, method_name, functional=functional, basis=basis_set_textbox,
                                     n_states=n_states, maxdim=maxdim, root=root, charge=charge, multiplicity=multiplicity,
                                     solvation=use_solvation, solvation_model=solvation_model, solvent=solvent, n_proc=n_cores, memory=memory)
        else: # calculation_type == "NMR Spectrum"
            write_nmr_orca_input(mol, input_file_path, method_type, functional=functional, basis=basis_set_textbox,
                                 spin_spin_coupling=spin_spin_coupling, charge=charge, multiplicity=multiplicity,
                                 solvation=use_solvation, solvation_model=solvation_model, solvent=solvent, n_proc=n_cores, memory=memory)

        status = "Input file generated."
        return f"<span style='color:green;'>{status}</span>", get_files_in_working_directory(working_directory_path)
    except Exception as exc:
        status = exc
        return f"<span style='color:red;'>{status}</span>", get_files_in_working_directory(working_directory_path)

def on_n_cores_change(n_cores, memory):
    """Re-cap the per-core memory slider when the core count changes.

    Keeps ``nprocs x maxcore`` within physical RAM by lowering the slider's maximum,
    clamping the current value down to it when the new ceiling is smaller.
    """
    per_core_max = max_memory_per_core(n_cores)
    return gr.update(maximum=per_core_max, value=min(memory, per_core_max))

def resolve_orca_executable(orca_executable_path):
    """Return an absolute path to the ORCA binary, or ``""`` if it cannot be located.

    ORCA aborts a parallel run (``%pal nprocs`` > 1) with *"ORCA has to be called with
    full pathname"* when it is invoked as a bare command name, so a blank entry falls
    back to the PATH lookup done at import and any relative entry is resolved and made
    absolute before it reaches ``subprocess``.
    """
    candidate = (orca_executable_path or "").strip()
    if not candidate:
        candidate = default_orca_path
    if not candidate:
        return ""
    if not os.path.isabs(candidate):
        candidate = shutil.which(candidate) or candidate
    return os.path.abspath(candidate)

def on_run_calculation(working_directory_path, input_file_name, orca_executable_path):
    """Run the external ORCA binary on ``input_file_name`` (a ``.inp`` in the wd).

    ORCA is invoked by absolute path (see :func:`resolve_orca_executable`) because it
    rejects parallel runs launched via a bare command name. stdout is captured to a
    same-named ``.log``; an existing ``.log`` from a previous run is first renamed to
    ``.log.bak`` so an accidental click cannot destroy it. Returns
    ``(status_html, file_list)`` with the elapsed time on success; a non-zero ORCA exit
    (or a missing binary) is reported in the status span rather than raised. The running
    process is published in ``_orca_process`` so :func:`on_stop_calculation` can kill it.
    """
    global _orca_process, _orca_stopped_by_user

    if input_file_name is None or input_file_name=="":
        gr.Warning("Please choose an input file.")
        return "", get_files_in_working_directory(working_directory_path)

    orca_executable = resolve_orca_executable(orca_executable_path)
    if not orca_executable or not os.path.isfile(orca_executable):
        gr.Warning("Please set the full path to the ORCA executable in the Settings accordion.")
        status = f"ORCA executable not found: {orca_executable or '(empty)'}"
        return f"<span style='color:red;'>{status}</span>", get_files_in_working_directory(working_directory_path)

    try:
        input_file_path = os.path.join(working_directory_path, input_file_name)
        output_file_path = os.path.join(working_directory_path, os.path.splitext(input_file_name)[0] + ".log")

        # Opening the .log below truncates it, so preserve any previous run's results
        # first: an accidental click on Run would otherwise destroy them. os.replace
        # overwrites an older backup atomically and works on Windows too.
        backed_up = False
        if os.path.isfile(output_file_path):
            os.replace(output_file_path, output_file_path + ".bak")
            backed_up = True

        # Run calculation with ORCA. The full pathname is required for parallel runs.
        # Popen rather than subprocess.run so the Stop button has something to kill;
        # the handle is published under the lock and cleared again on the way out.
        cmd = [orca_executable, input_file_path]
        start = time.time()
        print(f"Running command: {' '.join(cmd)} > {output_file_path}")
        with open(output_file_path, 'w') as outfile:
            process = subprocess.Popen(cmd, stdout=outfile, stderr=subprocess.STDOUT)
            with _orca_lock:
                _orca_process = process
                _orca_stopped_by_user = False
            try:
                return_code = process.wait()
            finally:
                with _orca_lock:
                    _orca_process = None
                    stopped = _orca_stopped_by_user
        end = time.time()
        duration = end - start

        log_name = os.path.basename(output_file_path)
        if stopped:
            status = f"Calculation stopped after {round(duration, 3)} s. Partial output in {log_name}."
            return f"<span style='color:red;'>{status}</span>", get_files_in_working_directory(working_directory_path)
        if return_code != 0:
            status = f"Error running calculation: ORCA exited with code {return_code}. See {log_name}."
            return f"<span style='color:red;'>{status}</span>", get_files_in_working_directory(working_directory_path)

        status = f"Calculation finished ({round(duration, 3)} s)."
        if backed_up:
            status += f" Previous log saved as {log_name}.bak."
        return f"<span style='color:green;'>{status}</span>", get_files_in_working_directory(working_directory_path)
    except Exception as exc:
        status = f"Error running calculation: {exc}"
        return f"<span style='color:red;'>{status}</span>", get_files_in_working_directory(working_directory_path)

def on_calculation_started():
    """Grey out the Run button for the duration of a calculation.

    Chained *ahead* of :func:`on_run_calculation` rather than folded into it: that
    handler blocks until ORCA exits, so anything it returned would only reach the
    browser once the run was already over.
    """
    return gr.update(interactive=False)

def on_calculation_finished():
    """Re-enable the Run button once no calculation is running.

    Chained with ``.then`` (not ``.success``) so the button comes back even when the
    run ended in an error or was stopped, which would otherwise leave it stuck.
    """
    return gr.update(interactive=True)

def on_stop_calculation():
    """Terminate the ORCA process started by the Run button, if one is running.

    ORCA's parallel driver spawns MPI workers, so terminating only the process we
    launched would orphan them and leave the cores busy; the whole process tree is
    signalled instead (SIGTERM first, SIGKILL for anything still alive after 5 s).
    The run handler reports the outcome once ``process.wait()`` returns, so this only
    needs to acknowledge the request.
    """
    global _orca_stopped_by_user

    with _orca_lock:
        process = _orca_process
        if process is None or process.poll() is not None:
            gr.Warning("No calculation is running.")
            return ""
        _orca_stopped_by_user = True

    try:
        parent = psutil.Process(process.pid)
        targets = parent.children(recursive=True) + [parent]
    except psutil.NoSuchProcess:
        targets = []

    for target in targets:
        try:
            target.terminate()
        except psutil.NoSuchProcess:
            pass
    _, alive = psutil.wait_procs(targets, timeout=5)
    for target in alive:
        try:
            target.kill()
        except psutil.NoSuchProcess:
            pass

    return "<span style='color:red;'>Stopping calculation...</span>"

def calculation_tab_content(working_directory_path_state, working_directory_file_list_state, status_markdown):
    """Build the "Calculation" tab (settings, input generation, run) and wire events.

    Takes the shared path/file-list states and status line; returns the tab component.
    """
    with gr.Tab("Calculation") as calculation_tab:
        with gr.Accordion("Settings", open=False):
            with gr.Row():
                orca_path_textbox = gr.Textbox(label="ORCA executable path", value=default_orca_path,
                                               placeholder="/full/path/to/orca",
                                               info="Full pathname is required — ORCA refuses parallel runs launched by bare command name.")
            with gr.Row():
                n_cores_slider = gr.Slider(label="Number of cores", value=max_cores, minimum=1, maximum=max_cores, step=1)
                memory_slider = gr.Slider(label="Memory per core (GB)", value=max_memory_per_core(max_cores),
                                          minimum=1, maximum=max_memory_per_core(max_cores), step=1,
                                          info=f"ORCA %maxcore is per process — cores x this must fit in {max_memory} GB of RAM.")
        with gr.Accordion("Generate Input File"):
            with gr.Row():
                with gr.Column(scale=1):
                    input_structure_file_dropdown = gr.Dropdown(label="Input structure", choices=[""], value="", interactive=False)
                with gr.Column(scale=4):
                    calculation_type_radio = gr.Radio(label="Type of calculation", value="Single-Point", choices=["Single-Point", "Geometry Optimization", "Frequency", "Time-Dependent Density Functional Theory", "NMR Spectrum"])
            with gr.Row():
                with gr.Column(scale=1):
                    with gr.Row():
                        mm_checkbox = gr.Checkbox(label="Optimize geometry with molecular mechanics", value=False)
                    with gr.Row():
                        with gr.Column(scale=1):
                            force_field_dropdown = gr.Dropdown(label="Force field", value="MMFF", choices=["MMFF", "UFF"], visible=False)
                        with gr.Column(scale=1):
                            max_iters_slider = gr.Slider(label="Max iterations", value=200, minimum=0, maximum=1000, step=1, visible=False)
                    with gr.Row():
                        solvation_checkbox = gr.Checkbox(label="Solvation", value=False)
                    with gr.Row():
                        with gr.Column(scale=1):
                            solvation_dropdown = gr.Dropdown(label="Solvation model", value="cpcm", choices=[("SMD", "smd"), ("CPCM", "cpcm")], visible=False)
                        with gr.Column(scale=1):
                            solvent_dropdown = gr.Dropdown(label="Solvent", value="water", choices=["water", ("DMSO", "dmso"),  "nitromethane", "acetonitrile", "methanol", "ethanol", "acetone", "dichloromethane",
                                                                                                    "dichloroethane", ("THF", "thf"), "aniline", "chlorobenzene", "chloroform", ("diethyl ether", "diethylether"),
                                                                                                    "toluene", "benzene", ("CCl4", "ccl4"), "cyclohexane", "heptane"], allow_custom_value=True, visible=False)
                with gr.Column(scale=1):
                    with gr.Row():
                        method_type_dropdown = gr.Dropdown(label="Type of method", value="DFT", choices=["HF", "DFT", "Semi-empirical", "MP2", "CCSD", "CCSD(T)"])
                        method_name_dropdown = gr.Dropdown(label="Method", choices=[], visible=False)
                    with gr.Row():
                        functional_dropdown = gr.Dropdown(label="Functional", value="B3LYP", choices=FUNCTIONALS, allow_custom_value=True)
                        basis_set_dropdown = gr.Dropdown(label="Basis set", value=DEFAULT_BASIS_SET, choices=BASIS_SETS, allow_custom_value=True)
                    with gr.Row():
                        n_states_slider = gr.Slider(label="Number of excited states", value=10, minimum=5, maximum=100, step=1, visible=False)
                        maxdim_slider = gr.Slider(label="Dimension of expansion space", value=5, minimum=1, maximum=10, step=1, visible=False)
                        root_slider = gr.Slider(label="Root state", value=0, minimum=0, maximum=99, step=1, visible=False)
                        spin_spin_coupling_checkbox = gr.Checkbox(label="Compute spin-spin couplings", value=False, visible=False)
                    with gr.Row():
                        charge_slider = gr.Slider(label="Charge", value=0, minimum=-2, maximum=2, step=1)
                        multiplicity_dropdown = gr.Dropdown(label="Multiplicity", value=1, choices=[("Singlet", 1), ("Doublet", 2),
                                                                                                    ("Triplet", 3), ("Quartet", 4),
                                                                                                    ("Quintet", 5), ("Sextet ", 6)])
                with gr.Column(scale=1):
                    input_file_name_textbox = gr.Textbox(label="File name", value="single_point")
                    generate_input_file_button = gr.Button(value="Generate input file")
        with gr.Accordion("Run calculation"):
            with gr.Row():
                with gr.Column(scale=1):
                    input_file_name_dropdown = gr.Dropdown(label="Input file", choices=[""], value="", interactive=False)
                with gr.Column(scale=2):
                    with gr.Row():
                        run_button = gr.Button("Run", variant="primary")
                        stop_button = gr.Button("Stop", variant="stop")

        working_directory_file_list_state.change(on_working_directory_file_list_change, [working_directory_file_list_state, input_file_name_textbox], [input_structure_file_dropdown, input_file_name_dropdown])
        calculation_type_radio.change(on_change_calculation_type, [calculation_type_radio, basis_set_dropdown], [input_file_name_textbox, method_type_dropdown, functional_dropdown, basis_set_dropdown, n_states_slider, maxdim_slider, root_slider, spin_spin_coupling_checkbox])
        mm_checkbox.change(on_mm_checkbox_change, mm_checkbox, [force_field_dropdown, max_iters_slider])
        solvation_checkbox.change(on_solvation_checkbox_change, solvation_checkbox, [solvation_dropdown, solvent_dropdown])
        method_type_dropdown.change(on_method_type_change, method_type_dropdown, [method_name_dropdown, functional_dropdown, basis_set_dropdown])
        n_cores_slider.change(on_n_cores_change, [n_cores_slider, memory_slider], memory_slider)
        generate_input_file_button.click(on_generate_input_file, [working_directory_path_state, input_structure_file_dropdown, calculation_type_radio,
                                                                  mm_checkbox, force_field_dropdown, max_iters_slider, solvation_checkbox, solvation_dropdown, solvent_dropdown,
                                                                  method_type_dropdown, method_name_dropdown, functional_dropdown, basis_set_dropdown, n_states_slider, maxdim_slider, root_slider, spin_spin_coupling_checkbox, charge_slider, multiplicity_dropdown,
                                                                  n_cores_slider, memory_slider, input_file_name_textbox],
                                                                 [status_markdown, working_directory_file_list_state])
        run_button.click(on_calculation_started, None, run_button) \
                  .then(on_run_calculation, [working_directory_path_state, input_file_name_dropdown, orca_path_textbox], [status_markdown, working_directory_file_list_state]) \
                  .then(on_calculation_finished, None, run_button)
        # concurrency_limit=None: without it this event would queue behind the running
        # calculation it is meant to interrupt, and only fire once ORCA had finished.
        stop_button.click(on_stop_calculation, None, status_markdown, concurrency_limit=None)
        
    return calculation_tab