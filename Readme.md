![Workflow](./images/workflow.png)

## Introduction
This web UI is for simple computational chemistry calculations with [ORCA](https://www.faccts.de/orca/):

* Single-Point Calculation

* Geometry Optimization (ground state or a chosen excited state)

* Frequency Analysis

* Absorption/Emission Spectrum Prediction

* ECD (Electronic Circular Dichroism) Spectrum Prediction

* NMR Prediction

<p style="display: flex; justify-content: space-between; gap: 10px;">
  <img src="./images/webui_1.png" alt="Description 1" style="width: 100%;" />
</p>
<p style="display: flex; justify-content: space-between; gap: 10px;">
  <img src="./images/webui_2.png" alt="Description 1" style="width: 49%;" />
  <img src="./images/webui_3.png" alt="Description 2" style="width: 49%;" />
</p>
<p style="display: flex; justify-content: space-between; gap: 10px;">
  <img src="./images/webui_4.png" alt="Description 1" style="width: 59%;" />
  <img src="./images/webui_5.png" alt="Description 2" style="width: 39%;" />
</p>

> **Note:** This web UI drives the external [ORCA](https://www.faccts.de/orca/)
> program but does **not** bundle it. Install ORCA separately (free for academic use)
> and make sure the `orca` executable is on your `PATH`; only then can calculations be
> run. Input-file generation and result parsing work without it.

## Installation

The recommended way to install is the [PyPI package](https://pypi.org/project/orca-webui/) (Python 3.10+).

First create and activate a virtual environment so the dependencies stay isolated from your system Python:

***For Linux and macOS***

```
python3 -m venv orca-env
source orca-env/bin/activate
```

***For Windows***

```
python -m venv orca-env
orca-env\Scripts\activate.bat
```

Then install the package from PyPI:

```
pip install orca-webui
```

This pulls in all Python dependencies (Gradio, RDKit, cclib, plotly, nglview, ...).

### Built-in molecular editor

The Conformer generation tab includes a molecular drawing editor based on Ketcher
3.18.0 and the app's own Gradio HTML integration. Its browser assets ship with the
package and are served locally, so drawing needs no CDN, Node installation, or
separate chemistry server. This release targets Gradio 6.30.0.

- Draw a molecule, then click **Use drawing** to send it to the **SMILES** field.
- To edit a SMILES string visually, enter it in the field and click **Load SMILES
  into editor**. Loading an empty field clears the drawing.
- Click **Generate** to create conformers from the current **SMILES** field. After
  changing the drawing, click **Use drawing** again before generating.

The app validates molecule text with RDKit before accepting it. You can also enter
SMILES and generate conformers directly without using the drawing editor.

### Newer ORCA versions (6.1+) and cclib

Results are parsed with [cclib](https://cclib.github.io/). The latest **released** cclib
(1.8.1) predates ORCA 6.1 and cannot parse its output — loading such a result file fails
with an error like *"cclib could not parse this file ... list index out of range."* The
fix is already in cclib's development version, so until the next cclib release, install it
from git:

```
pip install --upgrade "cclib @ git+https://github.com/cclib/cclib.git"
```

If your calculations come from ORCA 6.0 or earlier, the released cclib pulled in by
`pip install orca-webui` is sufficient.

## Start web UI

Run the console command from the directory where you want your `data/` and `static/`
folders to be created:

```
orca-webui
```

It starts a local server on the first free port at/after 7860 (open the printed URL in
your browser).

### Running a calculation

The **Settings** accordion on the Calculation tab controls how ORCA is invoked:

* **ORCA executable path** — defaults to whatever `which orca` finds. ORCA refuses to
  start a parallel run when it is launched by bare command name (*"For parallel runs
  ORCA has to be called with full pathname"*), so the full path is used. Set it by hand
  if ORCA is not on your `PATH`.
* **Number of cores** — capped at half the logical CPU count, which is the physical core
  count on a hyperthreaded machine; ORCA fails to start above that.
* **Memory per core (GB)** — this is ORCA's `%maxcore`, a **per-process** budget, so the
  calculation may claim `cores x maxcore` in total. The slider's ceiling therefore
  shrinks as the core count rises, keeping the product inside physical RAM.

**Run** displays **Running calculation...** in the shared status area and greys out
until the calculation finishes. The status then shows completion, failure, or a stop
message. **Stop** terminates ORCA and its child processes. Only one calculation can
run at a time in an app instance.

ORCA runs inside the selected working directory, so relative input references and
scratch files stay with the job. Before rerunning an input, its existing `.log` is
renamed to `.log.bak`, replacing any older backup. If ORCA fails to launch, the
previous log is restored. Restarting the app preserves calculation logs; startup
cleanup is limited to transient viewer files under `static/`.

### Working with files and results

Working-directory and output names must stay inside their selected directory.
Uploads and conformer generation reject existing output names; choose a new name to
keep both versions. **Save Text File** saves the file loaded into the editor, even if
you subsequently select another row. Switching directories clears the editor,
structure viewer, and loaded results.

XYZ imports accept both standard XYZ (atom count, comment, then coordinates) and
the app's legacy format (charge and multiplicity, then coordinates). Invalid elements,
nonfinite coordinates, and mismatched atom counts are reported as errors.

Result tables show orbital energies in hartree, including separate alpha and beta
orbitals for unrestricted calculations. Available energy and frequency data can be
displayed even when optional fields are missing. NMR multiplets use Pascal weights
directly to keep large equivalent partner groups responsive.

See [CHANGELOG.md](./CHANGELOG.md) for the changes prepared for version 0.1.12.

### Optimization and spectra are separate steps

**TD-DFT is a single point.** It computes excitation energies, oscillator strengths and
rotatory strengths at exactly the geometry you select — it does not optimize. Do the
optimization first with **Geometry Optimization**, then select the resulting `.log` as
the input structure for the TD-DFT step.

**Geometry Optimization** has an *Optimize* choice:

* **Ground state** — an ordinary optimization.
* **Excited state** — adds a `%tddft` block with `iroot N`, so ORCA follows the gradient
  of excited state *N* and returns that state's relaxed geometry.

Which one you want depends on the spectrum:

| Spectrum | Geometry to run TD-DFT on |
| --- | --- |
| Absorption, ECD | ground-state optimized |
| Emission (fluorescence) | excited-state optimized (`iroot 1` for S1) |

Emission happens *after* the excited state relaxes, so a TD-DFT single point on a
ground-state geometry gives a vertical excitation — the same number as the first
absorption band, with no Stokes shift. Optimize the excited state first, then load that
`.log` in the Result tab and press **Show emission spectrum**.

> **Note:** an excited-state optimization needs the third derivative of the
> exchange-correlation functional, which ORCA 6.1 does not provide natively for every
> functional. `B3LYP` and `wB97X-D4` abort with *"Invalid or unknown value for Exchange
> in DFT GXC. Please try using LibXC instead!"*, while `PBE` works. This limit is ORCA's,
> not the UI's, and it only affects excited-state optimization — ground-state
> optimization and TD-DFT single points work with any functional.

### ECD spectra

Any TD-DFT result offers **Show ECD spectrum** in the Result tab: ORCA 6.1 reports
rotatory strengths (its *CD SPECTRUM* blocks) for TD-DFT automatically, so no extra
keyword is needed and existing TD-DFT logs work too. The spectrum broadens the *signed*
rotatory strengths, so positive and negative Cotton effects appear above and below zero.
Broadening is done in wavenumber space and only displayed against wavelength.

Two enantiomers differ only by the sign of every rotatory strength, so their ECD spectra
are mirror images. One calculation therefore characterizes both.

### Reading atom indices off the structure viewer

The structure viewer labels each atom with its element and its index (`C0`, `C1`, `O2`,
...). The numbering is **0-based**, matching both the order of the `* xyz` coordinate
block in the generated `.inp` and the nucleus indices ORCA reports in its output — so a
label read off the viewer can be used directly to identify a nucleus.

The **NMR signals** table in the Result tab uses the same labels in its *J-Couplings*
column: a row for `H0` coupled to `H1` by 1 Hz reads `1.0000 Hz (H1)`, with multiple
partners listed strongest-first.

## Development (from source)

You never need to publish to PyPI to work on the app — that step is only for *other*
people to install it.

```
git clone https://github.com/phatdatnguyen/orca-webui
cd orca-webui

# create and activate a virtual environment
python -m venv orca-env
source orca-env/bin/activate          # Windows: orca-env\Scripts\activate.bat

# install in editable mode with the dev extras (pytest, build, twine)
pip install -e ".[dev]"

# run the app and the tests
orca-webui
pytest
```

### Tests

Most of the suite uses synthetic ORCA logs and stand-in processes, covering input
generation, parsing, spectra, execution, file handling, the molecular editor, and UI
regressions. It does not require files from `data/`.

`tests/test_orca_integration.py` additionally drives the **real** ORCA binary — it
generates an input with the app's writers, runs it, and parses the output with the app's
parsers, which is the only way to catch ORCA rejecting an input we generate or changing
an output format we scrape. Those tests are marked `orca` and **skip automatically when
ORCA is not installed**:

```
pytest                                          # skips the ORCA tests if it can't find orca
ORCA_EXECUTABLE=~/orca/orca_6_1_1/orca pytest   # ORCA installed but not on PATH
pytest -m "not orca"                            # skip them explicitly
```

The integration tests include water single points, methane NMR, and hydrogen peroxide
ground/excited-state optimizations and TD-DFT/ECD calculations. They use `nprocs 1`,
so they need no MPI. Runtime depends on the installed ORCA version and hardware;
parsing also requires a compatible cclib version.

### Maintaining the editor assets

The editor integration lives in `src/orca_webui/molecular_editor.py` and
`src/orca_webui/frontend/molecular_editor.js`. The pinned Ketcher standalone files,
including their license and dependency notices, live in
`src/orca_webui/assets/ketcher/` and are included in both distributions. These are
permanent package assets, separate from generated viewers in the runtime `static/`
directory.

Maintainers can restore a missing asset directory with:

```bash
python scripts/vendor_molecular_editor.py
```

The helper verifies the pinned download checksums and installs into a staging
directory first. It checks an existing installation instead of overwriting it.
Normal package installation and app startup do not download editor files.

### Editable installs

The `-e` (editable) flag is the important part: pip records a path hook pointing at
`src/`, so your edits are picked up on the next restart and you never reinstall. A plain
`pip install .` *copies* the package into `site-packages`, which is why changes appear to
be ignored until you reinstall. Note that `pip install --upgrade orca-webui` (e.g. to
check a published release) also replaces the editable link with a copy — use a throwaway
venv for that. To check which one you have:

```
python -m pip show orca-webui
```

Look for `Editable project location` pointing to this repository. If it is missing,
restore the editable install:

```
pip uninstall -y orca-webui
pip install -e .
```

### Running without installing at all

Because `orca_webui` is a regular package with relative imports, you can run it straight
from the source tree:

```
PYTHONPATH=src python -m orca_webui.app
```

Use `python -m orca_webui.app`, not `python src/orca_webui/app.py` — the latter fails on
the package-relative imports. Run it from the repository root, since `build_app()` creates
`data/` and `static/` relative to the current working directory.

### Auto-reload on file save

`build_app()` returns the FastAPI app, so it can be used directly as a uvicorn factory:

```
PYTHONPATH=src python -m uvicorn orca_webui.app:build_app --factory --reload --reload-dir src --port 7860
```

The server restarts whenever anything under `src/` changes. Note that `--reload` bypasses
`main()`, so the automatic "first free port at/after 7860" selection does not apply — pick
the port yourself.

### Copying the project to another directory

Virtual environments are **not** relocatable: `orca-env/pyvenv.cfg` and the activate
scripts contain the absolute path they were created at. After copying or moving the
project, delete the old environment and recreate it:

```
rm -rf orca-env
python -m venv orca-env
source orca-env/bin/activate
pip install -e ".[dev]"
```
