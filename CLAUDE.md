# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A Gradio-based web UI that wraps the [ORCA](https://www.faccts.de/orca/) quantum chemistry program for simple computational chemistry workflows: conformer generation, then Single-Point / Geometry Optimization / Frequency / TDDFT (absorption/emission/ECD) / NMR calculations, then result visualization (energies, MOs, IR spectra, UV-Vis/ECD spectra, and simulated 1H/13C NMR spectra with multiplet analysis).

ORCA itself is an external binary launched with `subprocess.Popen` using an absolute executable path and the selected working directory as `cwd` (see `on_run_calculation` in `src/orca_webui/calculation.py`). It is **not** bundled: install it separately and put it on `PATH` or set its executable path in the Calculation tab. Input-file generation and result parsing work without it.

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

# Build & validate the prepared 0.1.11 distributions
python -m build
python -m twine check dist/orca_webui-0.1.11-py3-none-any.whl dist/orca_webui-0.1.11.tar.gz
```

`orca_webui.app.main()` auto-selects the first free port starting at 7860 and serves on `127.0.0.1`. Notable pins in `pyproject.toml`: `gradio==5.50.0` and `nglview==4.0` (the app targets those exact versions); `cclib>=1.8` from PyPI (the old README installed cclib from git — PyPI packages cannot declare git dependencies, so a released cclib is required).

The `orca-env/` venv lives inside the repo and is gitignored, as are build artifacts (`dist/`, `build/`, `*.egg-info/`). Keep `pyproject.toml` and `src/orca_webui/__init__.py` version values aligned. The next release is 0.1.11; changes are recorded in `CHANGELOG.md`. See `Readme.md` for release preparation and publication steps. Building and checking distributions does not publish them.

### Tests

`tests/` covers ORCA input generation, molecule I/O, spectra, the NMR pipeline, process execution, file handling, and UI regressions. Most tests use synthetic logs and stand-in processes; none depend on gitignored `data/` outputs. `tests/test_orca_integration.py` also runs real single-point, optimization, TDDFT/ECD, and NMR calculations. These tests carry the `orca` marker and skip when ORCA is unavailable; set `ORCA_EXECUTABLE` to its full path or use `pytest -m "not orca"` to exclude them. Tests needing the UI stack or cclib use `pytest.importorskip`. `cclib` is imported lazily in `utils.mol_from_orca_file`, so most engine tests work without it.

## Architecture

The app is a FastAPI host with a Gradio Blocks UI mounted at `/` and a `StaticFiles` mount at `/static`, assembled in `app.py`'s `build_app()` (entry point `main()`). `static/` holds transient HTML/cube/xyz artifacts (e.g. nglview structure viewers embedded via `<iframe>`). `_cleanup_transient_files` removes only these artifacts on startup and preserves calculation logs; cleanup must not follow symlinks outside `static/`. Both `data/` (user working directories) and `static/` are created **relative to the current working directory**, or under the explicit `build_app(run_dir=...)` directory, not the install location — so an installed package never writes into `site-packages`. Viewers use unique HTML filenames. `styles.css` is the exception: it's package data, read from beside the module.

### Working-directory model

Everything is organized around a **working directory** — a subfolder under `data/` selected in the left column. All file I/O (structure files, generated `.inp` inputs, ORCA `.log` outputs, exported `.csv`) happens inside it. `working_directory.py` owns this left column and the file browser/viewer.

State flows through Gradio `gr.State` objects created in `working_directory_blocks()` and threaded into every tab:

- `working_directory_path_state` — the current working directory path.
- `working_directory_file_list_state` — the file list; **this is the central event bus.** Any handler that changes files on disk returns a fresh `get_files_in_working_directory(...)` list into this state. Its `.change` event fires re-population handlers in each tab (`on_working_directory_file_list_change`), which is how the dropdowns of available structures / input files / result files stay in sync across tabs. When adding an operation that writes files, return the refreshed file list into this state to propagate the update.
- `status_markdown` — the shared Calculation/Result status line created in `app.py`. `on_calculation_started` immediately displays orange `Running calculation...` and disables Run before the blocking callback. Completion, failure, or stopping replaces that text; re-enabling Run must preserve the final status. The conformer tab has its own status component.

### UI structure (three tabs, one file each)

- `conformer_generation.py` — draw a molecule (optional `gradio_molecule2d`) or type SMILES → RDKit embeds and optimizes conformers → writes `.xyz`/`.pdb`/`.mol` files. Without the editor dependency, the tab uses SMILES entry. Settings and output collisions are checked before embedding; unsuccessful optimizations are excluded.
- `calculation.py` — pick a structure file + calculation type + method/basis/solvation → generate an ORCA `.inp` file → optionally launch ORCA, capturing stdout/stderr to a `.log`. A process lock permits only one active calculation per app instance; Stop terminates the process tree. Existing logs rotate to `.log.bak` (replacing older backups), and launch failure restores the previous log. Validate paths and the executable before rotating files. The calculation-type radio drives control visibility via `on_change_calculation_type` / `on_method_type_change`.
- `result.py` — loads a `.log` via cclib, then reveals available Energy, Geometry Optimization, Frequency, Absorption/Emission/ECD, and NMR results. Convert cclib SCF and MO energies from eV to hartree; include both spin channels for unrestricted results. Optional fields may be absent in partial logs. `on_load_result_file` and `_empty_result_outputs` must stay index-aligned with `nmr_result_outputs`; update all three when changing outputs. Directory changes and load errors clear cached data and disable exports.

### `utils.py` — the shared engine

Pure, UI-free functions used by all tabs. Two responsibilities:

1. **Molecule I/O & ORCA input writing.** `mol_from_xyz_file` accepts standard atom-count/comment headers and legacy charge/multiplicity headers, validating every coordinate row. `mol_from_structure_file` is the shared viewer/input loader; `.log` geometry comes from `mol_from_orca_geometry` in ORCA nucleus order, falling back to `mol_from_orca_file` through cclib. `add_bonds` infers bonds from distances and covalent radii for molecules lacking connectivity. The `write_*_orca_input` family generates calculation inputs. All input writers funnel through `get_orca_method_keyword` + `get_orca_header` (writes `%pal`/`%maxcore`/`%tddft`/`%cpcm` blocks) + `get_orca_coordinates`.

2. **Spectrum simulation.** `generate_ir_spectrum_interactive` (Lorentzian), `generate_absorption_emission_spectrum_interactive` (Gaussian), `generate_ecd_spectrum_interactive` (signed rotatory strengths), and the NMR pipeline (see below). Spectra are Plotly figures.

### NMR pipeline (the most intricate part)

Simulating a realistic 1H/13C NMR spectrum from a static ORCA calculation requires care around **atom ordering** and **chemical equivalence**:

- **Atom order gotcha:** cclib *reorders* atoms, but ORCA's NMR nucleus indices follow the log's `CARTESIAN COORDINATES` block. So geometry for equivalence analysis is parsed directly from the log via `parse_orca_geometry` (NOT cclib). Mixing the two orderings will mis-assign shifts.
- `parse_nmr_shielding_constants` / `parse_nmr_jcouplings` — regex-scrape shielding and J-coupling values from the ORCA `.log` (ORCA 6.x property-block format).
- `perceive_bonds` → `morgan_ranks` → `compute_hydrogen_equivalence` / `compute_carbon_equivalence` — infer bonds from covalent radii, then Morgan-style graph refinement assigns topological-symmetry ranks so equivalent nuclei (e.g. the three protons of a CH3, or symmetric methyls) merge into one averaged peak. Equivalence is **topological only** — diastereotopic protons are deliberately treated as equivalent.
- `build_nmr_peak_table` — groups nuclei into averaged multiplet peaks; multiplicity comes from the *partner group size* with an averaged J, giving first-order (Pascal) multiplets. `generate_nmr_spectrum_interactive` renders them (1H with J-splitting; 13C proton-decoupled as singlets). Generate the `n+1` Pascal lines directly for `n` equivalent partners; repeated binary splitting creates exponentially many duplicates. Keep the chemical-shift axis descending.

NMR reference shieldings (TMS defaults: H=32.9, C=203.4 ppm) are user-editable textboxes in the Result tab; the "Refresh" button re-runs `on_load_result_file`, reading the selected log again with the new references. Blank, invalid, or nonfinite references fall back to the defaults.

## Conventions

- **Handler signature:** every UI callback is a top-level `on_*` function taking Gradio component values positionally and returning values positionally into the outputs list wired in the tab's `*_tab_content` / `*_blocks` builder. Keep the wiring (`.click`/`.change`/`.select` calls at the bottom of each builder) index-aligned with handler returns.
- **Errors surface to the user:** handlers catch failures and return an escaped red status span plus refreshed outputs, or show `gr.Warning` for file operations. Preserve return arity and avoid secondary errors while refreshing a missing directory.
- **File safety:** use `resolve_working_file` for user-selected files and editable output names. Reject paths and symlinks escaping the working directory. Uploads and conformer generation must reject output collisions. Text saves use an atomic replacement and the path of the file loaded into the editor, not the currently selected browser row. Clear file selections, viewers, editor contents, and cached results when the working directory changes.
- Generated artifacts (`.inp`, `.log`, `.csv`, `.xyz`, `.html`, `.cube`, etc.) and `data/`, `static/`, `orca-env/` are all gitignored — never commit them.
