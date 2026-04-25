## Introduction
This web UI is for simple computational chemistry calculations with [Orca](https://www.faccts.de/orca/):

* Single-Point Calculation

* Geometry Optimization

* Frequency Analysis

* Absorption/Emission Spectrum Prediction

* NMR Prediction


## Installation
- Clone this repo: Open terminal

```
git clone https://github.com/phatdatnguyen/orca-webui
```

- Create and activate virtual environment:

***For Linux and MacOS***

```
cd orca-webui
python3 -m venv orca-env
source orca-env/bin/activate
```
***For Windows***

```
cd orca-webui
python -m venv orca-env
orca-env\Scripts\activate.bat
```

- Install packages:

```
pip install psutil
pip install rdkit
python -m pip install git+https://github.com/cclib/cclib.git
pip install plotly
pip install gradio_molecule2d
pip uninstall gradio
pip install gradio==5.50.0
pip install nglview==4.0
```

## Start web UI
To start the web UI:

***For Linux and MacOS***

```
source orca-env/bin/activate
python3 webui.py
```

***For Windows***
```
orca-env\Scripts\activate.bat
python webui.py
```