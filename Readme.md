## Introduction
This web UI is for simple computational chemistry calculations with [Orca](https://www.faccts.de/orca/):

* Single-Point Calculation

* Geometry Optimization

* Frequency Analysis

* Absorption/Emission Spectrum Prediction

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

### Optional: 2D structure editor

The conformer tab can show a 2D drawing editor powered by `gradio_molecule2d`. That
package pins `gradio<5.0`, which conflicts with the `gradio 5.x` this app uses, so it
is **not** installed automatically. To enable the editor, install it without its
dependencies:

```
pip install gradio_molecule2d --no-deps
```

Without it, the app runs fine — you simply enter a SMILES string directly instead of
drawing.

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

## Development (from source)

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