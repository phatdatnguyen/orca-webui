# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A Gradio-based web UI that wraps the [ORCA](https://www.faccts.de/orca/) quantum chemistry program for simple computational chemistry workflows: conformer generation, then Single-Point / Geometry Optimization / Frequency / TDDFT (absorption/emission) / NMR calculations, then result visualization (energies, MOs, IR spectra, UV-Vis spectra, and simulated 1H/13C NMR spectra with multiplet analysis).

ORCA itself is an external binary invoked via `subprocess.run(["orca", input_file_path])` (see `on_run_calculation` in `src/orca_webui/calculation.py`). It is **not** bundled — the `orca` executable must be on `PATH` for calculations to actually run. Input-file generation and result parsing work without it.

## Packaging & layout

The project is a pip-installable package (`orca-webui` on PyPI, import name `orca_webui`) using a **src layout**: all modules live under `src/orca_webui/`. Metadata is in `pyproject.toml` (setuptools backend); `styles.css` ships as package data. The `orca-webui` console command maps to `orca_webui.app:main`.

## Commands

```bash
# Editable install with dev extras (pytest, build, twine) into a fresh venv
python -m venv orca-env
source orca-env/bin/activate          # Windows: orca-env\Scripts\activate.bat
pip install -e ".[dev]"

# Run the app (created data/ and static/ live in the current working directory)
orca-webui                            # or: python -m orca_webui.app

# Tests
pytest                                                        # whole suite
pytest tests/test_nmr_analysis.py                            # one file
pytest tests/test_nmr_analysis.py::test_peak_table_averages_equivalent_and_assigns_multiplicity  # one test

# Build & validate distributions
python -m build                       # -> dist/*.whl and dist/*.tar.gz
python -m twine check dist/*
```

`orca_webui.app.main()` auto-selects the first free port starting at 7860 and serves on `127.0.0.1`. Notable pins in `pyproject.toml`: `gradio==5.50.0` and `nglview==4.0` (the app targets those exact versions); `cclib>=1.8` from PyPI (the old README installed cclib from git — PyPI packages cannot declare git dependencies, so a released cclib is required).

The `orca-env/` venv lives inside the repo and is gitignored, as are build artifacts (`dist/`, `build/`, `*.egg-info/`).

### Tests

`tests/` is a pytest suite covering the `utils.py` engine (ORCA input generation, molecule I/O, spectra, and the full NMR parsing/equivalence/multiplet pipeline) plus regression tests for known fixes. Fixtures are self-contained (synthetic ORCA-log text in `conftest.py`) — they do **not** rely on the gitignored `data/` outputs. Tests needing the UI stack (`gradio`, `nglview`) or `cclib` use `pytest.importorskip`, so they skip cleanly in a minimal env and run fully inside the venv. `cclib` is imported lazily in `utils.mol_from_orca_file` so the rest of the engine (and most tests) import without it.

## Architecture

The app is a FastAPI host with a Gradio Blocks UI mounted at `/` and a `StaticFiles` mount at `/static`, assembled in `app.py`'s `build_app()` (entry point `main()`). `static/` holds transient HTML/cube/xyz artifacts (e.g. nglview structure viewers embedded via `<iframe>`); `build_app()` clears these on startup via `_cleanup_transient_files`. Both `data/` (user working directories) and `static/` are created **relative to the current working directory**, not the install location — so an installed package never writes into `site-packages`. `styles.css` is the exception: it's package data, read from beside the module.

### Working-directory model

Everything is organized around a **working directory** — a subfolder under `data/` selected in the left column. All file I/O (structure files, generated `.inp` inputs, ORCA `.log` outputs, exported `.csv`) happens inside it. `working_directory.py` owns this left column and the file browser/viewer.

State flows through Gradio `gr.State` objects created in `working_directory_blocks()` and threaded into every tab:
- `working_directory_path_state` — the current working directory path.
- `working_directory_file_list_state` — the file list; **this is the central event bus.** Any handler that changes files on disk returns a fresh `get_files_in_working_directory(...)` list into this state. Its `.change` event fires re-population handlers in each tab (`on_working_directory_file_list_change`), which is how the dropdowns of available structures / input files / result files stay in sync across tabs. When adding an operation that writes files, return the refreshed file list into this state to propagate the update.
- `status_markdown` — a shared status line; handlers return colored HTML spans (`<span style='color:green|red;'>...`) into it.

### UI structure (three tabs, one file each)

- `conformer_generation.py` — draw a molecule (`gradio_molecule2d`) or type SMILES → RDKit embeds N conformers → writes `.xyz`/`.pdb`/`.mol` files into the working directory.
- `calculation.py` — pick a structure file + calculation type + method/basis/solvation → generates an ORCA `.inp` file → optionally runs `orca` on it, capturing stdout to a `.log`. The calculation-type radio drives which UI controls are visible (TDDFT sliders, NMR spin-spin checkbox, etc.) via `on_change_calculation_type` / `on_method_type_change`.
- `result.py` — loads a `.log` via cclib, then conditionally reveals accordions (Energy, Geometry Optimization, Frequency, Absorption/Emission, NMR) based on what data is present in the parsed result. `on_load_result_file` returns a large positional tuple of `gr.update(...)` objects — the output list `nmr_result_outputs` and the return tuple must stay index-aligned; adding/removing an output means updating both.

### `utils.py` — the shared engine

Pure, UI-free functions used by all tabs. Two responsibilities:

1. **Molecule I/O & ORCA input writing.** `mol_from_xyz_file`, `mol_from_orca_file` (via cclib), `add_bonds` (distance/covalent-radius bond perception for RDKit mols that lack connectivity), and the `write_*_orca_input` family. All input writers funnel through `get_orca_method_keyword` + `get_orca_header` (writes `%pal`/`%maxcore`/`%tddft`/`%cpcm` blocks) + `get_orca_coordinates`.

2. **Spectrum simulation.** `generate_ir_spectrum_interactive` (Lorentzian), `generate_absorption_emission_spectrum_interactive` (Gaussian), and the NMR pipeline (see below). Spectra are Plotly figures.

### NMR pipeline (the most intricate part)

Simulating a realistic 1H/13C NMR spectrum from a static ORCA calculation requires care around **atom ordering** and **chemical equivalence**:

- **Atom order gotcha:** cclib *reorders* atoms, but ORCA's NMR nucleus indices follow the log's `CARTESIAN COORDINATES` block. So geometry for equivalence analysis is parsed directly from the log via `parse_orca_geometry` (NOT cclib). Mixing the two orderings will mis-assign shifts.
- `parse_nmr_shielding_constants` / `parse_nmr_jcouplings` — regex-scrape shielding and J-coupling values from the ORCA `.log` (ORCA 6.x property-block format).
- `perceive_bonds` → `morgan_ranks` → `compute_hydrogen_equivalence` / `compute_carbon_equivalence` — infer bonds from covalent radii, then Morgan-style graph refinement assigns topological-symmetry ranks so equivalent nuclei (e.g. the three protons of a CH3, or symmetric methyls) merge into one averaged peak. Equivalence is **topological only** — diastereotopic protons are deliberately treated as equivalent.
- `build_nmr_peak_table` — groups nuclei into averaged multiplet peaks; multiplicity comes from the *partner group size* with an averaged J, giving first-order (Pascal) multiplets. `generate_nmr_spectrum_interactive` renders them (1H with J-splitting; 13C proton-decoupled as singlets).

NMR reference shieldings (TMS defaults: H=32.9, C=203.4 ppm) are user-editable textboxes in the Result tab; the "Refresh" button re-runs `on_load_result_file` with new references without reloading the file.

## Conventions

- **Handler signature:** every UI callback is a top-level `on_*` function taking Gradio component values positionally and returning values positionally into the outputs list wired in the tab's `*_tab_content` / `*_blocks` builder. Keep the wiring (`.click`/`.change`/`.select` calls at the bottom of each builder) index-aligned with handler returns.
- **Errors surface to the user, never raise:** handlers wrap work in `try/except` and return a red status span + a refreshed file list rather than propagating exceptions. Match this pattern.
- Generated artifacts (`.inp`, `.log`, `.csv`, `.xyz`, `.html`, `.cube`, etc.) and `data/`, `static/`, `orca-env/` are all gitignored — never commit them.
