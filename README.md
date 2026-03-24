# QSIGA-CAL

QGIS plugin for configuration, calibration, and scenario simulation of hydrological, sedimentological, and water quality processes using the SIGA-CAL model.

Designed for watershed management and Nature-Based Solutions (NBS) prioritization workflows.

---

## Plugin Overview

```
Step 1 — Project Setup          →  Create project structure, load maps and time series
Step 2 — MODIS & Climatic       →  Download and process satellite LAI data (NASA MODIS)
Step 3 — Calibration            →  Calibrate and validate hydrological/water quality model
Step 4 — NBS Prioritization     →  Simulate NBS scenarios and prioritize intervention areas
Step 5 — Simulation             →  Run static/dynamic scenarios and export results
```

---

## Step-by-Step Description

### Step 1 — Project Setup

Creates the project folder structure and loads all base input files required by the SIGA-CAL model.

| | |
|---|---|
| **Input** | Topology file, raster maps (LULC, soil, DEM), Excel time series (precipitation, temperature, streamflow) |
| **Output** | `entradas/` folder with organized maps, series, calibration, and topology subfolders |
| **Key config** | Project name, basin name, output path, base scenario name |

---

### Step 2 — MODIS & Climatic Processing

Downloads MODIS MOD15A2H LAI product from NASA LAADS DAAC and processes it into phenology time series used by the model.

| | |
|---|---|
| **Input** | Date range, watershed boundary, NASA Earthdata Bearer token |
| **Output** | `03-Phenology/` folder with HDF4 files and processed LAI time series |
| **Key config** | NASA Earthdata token (valid 60 days — renew at [urs.earthdata.nasa.gov](https://urs.earthdata.nasa.gov/profile)) |

> A valid NASA Earthdata account is required. Generate a token at **Profile → Generate Token**.

---

### Step 3 — Calibration

Runs automated parameter calibration and validation using SPOTPY optimization algorithms. Covers four calibration targets independently.

| Calibration target | Available algorithms |
|---|---|
| Hydrological | DDS, LHS, SCE-UA |
| Hydrological validation | DDS, LHS, SCE-UA |
| Sediment | DDS, LHS, SCE-UA |
| Water quality | DDS, LHS, SCE-UA |

| | |
|---|---|
| **Input** | Calibration parameters file (`Parameters_*.txt`), observed time series |
| **Output** | Calibrated parameters, performance metrics table, time series plots |
| **Performance metrics** | NASH, IoAd, EAN, r, R², MSE, RMSE, RRMSE, COMPOS, t-Student |

---

### Step 4 — NBS Prioritization

Simulates a Nature-Based Solutions scenario and computes a spatial prioritization index to identify the most effective intervention areas within the watershed.

```
Step 4.1 — Configuration        →  Project setup and base simulation
Step 4.2 — Base simulation      →  Run Business-as-Usual scenario
Step 4.3 — NBS simulation       →  Run NBS scenario (with modified LULC or portfolio)
Step 4.4 — Prioritization       →  Compute delivery index (Dk) and rank polygons
```

| | |
|---|---|
| **Input** | Base and NBS simulation outputs, portfolio raster, polygon shapefile, variable weights |
| **Output** | Prioritization shapefile, Excel delivery table, coverage scenarios (20 %, 40 %, 80 %), `.tif` maps |
| **Key config** | Variable weights (Base Flow, Sediments, Nitrogen, Phosphorus, E. coli) |

---

### Step 5 — Simulation

Runs independent forward simulations for any configured scenario. Supports static LULC changes and dynamic multi-year scenarios. Exports results to Excel.

| Mode | Description |
|------|-------------|
| **Option 1 — Static** | Single simulation with a fixed LULC map |
| **Option 2A — Dynamic generation** | Generates annual LULC transition maps (scenarios A, B, C, D) |
| **Option 2B — Dynamic simulation** | Runs simulation using the generated dynamic scenario |
| **Option 3 — Export results** | Exports selected points and variables to Excel |

| | |
|---|---|
| **Input** | Execution file (`*.txt`), calibration file, sampler file, date range |
| **Output** | `salidas/{scenario}/` folder with simulation results, Excel files |

---

## Project File Structure

```
PluginQSigaCal/
├── __init__.py
├── QSIGA_CAL.py
├── metadata.txt
├── icon.png
├── logos/
└── Plugins/
    ├── common_functions/       ← Shared utilities (raster, files, SPOTPY setup)
    ├── requirements.txt        ← Python dependency list
    ├── Step_01/                ← Project setup module
    │   └── DATABASE.zip        ← Base model database (~4 GB) — downloaded separately
    ├── Step_02/                ← MODIS & climatic processing module
    ├── Step_03/                ← Calibration module
    ├── Step_04/                ← NBS prioritization module
    └── Step_05/                ← Simulation module
        └── bin/
            ├── SIGA.exe        ← SIGA-CAL model engine (requires Cygwin64)
            └── cnvrtr.exe
```

---

## Requirements

### Software

| Software | Version | Notes |
|----------|---------|-------|
| **QGIS** | 3.36 (64-bit) | Includes Python 3.12, GDAL, PyQt5 |
| **Windows** | 10 or 11 | 64-bit required |
| **PowerShell** | 5.1 or higher | Included in Windows 10/11 |
| **Cygwin64** | Latest | Required to run `SIGA.exe` — installed in `C:\cygwin64\` |

### Python packages

| Package | Version | Purpose |
|---------|---------|---------|
| `numpy` | **1.26.4** (exact) | Numerical operations |
| `scipy` | **1.13.0** (exact) | Scientific computing |
| `pandas` | ≥ 1.1 | Data manipulation |
| `spotpy` | ≥ 1.5.14 | Calibration algorithms (DDS, LHS, SCE) |
| `rasterio` | ≥ 1.2 | Raster I/O |
| `geopandas` | ≥ 0.9 | Vector spatial processing |
| `pysheds` | ≥ 0.3 | Watershed delineation |
| `numba` | ≥ 0.56 | Numerical acceleration |
| `matplotlib` | ≥ 3.3 | Plotting |
| `statsmodels` | ≥ 0.12 | Statistical analysis |
| `requests` | ≥ 2.25 | NASA MODIS downloads |
| `openpyxl`, `xlrd` | — | Excel read/write |
| `fiona`, `shapely`, `pyproj` | — | Geospatial support |

> `numpy==1.26.4` and `scipy==1.13.0` are pinned — other versions break `numba` and `pysheds`.
> GDAL and PyQt5 are bundled with QGIS and do not need separate installation.

---

## Installation

### 1 — Install Python dependencies

Open **OSGeo4W Shell** (comes with QGIS 3.36) as administrator:

> Start Menu → search `OSGeo4W Shell` → right-click → **Run as administrator**

```
python -m pip install "numpy==1.26.4" "scipy==1.13.0" pandas spotpy matplotlib openpyxl xlrd rasterio numba pysheds geopandas shapely pyproj fiona statsmodels requests seaborn scikit-learn
```

### 2 — Install Cygwin64

Download `cygwin64.zip` from the link below and save it to your Downloads folder:

> 📦 [https://tnc.box.com/s/e411u8hmda2qtqkw0nhp4ghbbmif5n04](https://tnc.box.com/s/e411u8hmda2qtqkw0nhp4ghbbmif5n04)

Open **PowerShell as administrator** and run:

```powershell
Expand-Archive -Path "$env:USERPROFILE\Downloads\cygwin64.zip" -DestinationPath "C:\" -Force
```

Verify that `C:\cygwin64\` exists before continuing.

### 3 — Copy the plugin to QGIS

Clone or download this repository, then run from the repository root:

```
xcopy /E /I /Y "PluginQSigaCal" "%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\QSIGA_CAL"
```

### 4 — Install the model database

Download `DATABASE.zip` from the same Box link and save it to your Downloads folder, then run:

```powershell
$dest = "$env:APPDATA\QGIS\QGIS3\profiles\default\python\plugins\QSIGA_CAL\Plugins\Step_01"
Move-Item -Path "$env:USERPROFILE\Downloads\DATABASE.zip" -Destination "$dest\DATABASE.zip"
```

> Step 3 must be completed first — the destination folder must already exist.

### 5 — Activate in QGIS

1. Open **QGIS 3.36**
2. Go to `Plugins → Manage and Install Plugins`
3. Click the **Installed** tab
4. Check **QSIGA_CAL**
5. Close the dialog

The **QSIGA-CAL** menu will appear in the toolbar.

---

## Updating

Pull the latest version and repeat step 3 only — no need to reinstall Cygwin or re-download the database:

```
git pull
xcopy /E /I /Y "PluginQSigaCal" "%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\QSIGA_CAL"
```

---

## Troubleshooting

### `No module named 'pysheds'` or similar

Repeat installation step 1.

### `numpy` version conflict

```
python -m pip install --force-reinstall "numpy==1.26.4" "scipy==1.13.0" --no-deps
```

### Plugin does not appear in QGIS

- Confirm the folder is named exactly `QSIGA_CAL` (uppercase)
- Confirm path: `C:\Users\<user>\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\QSIGA_CAL\`
- Restart QGIS

### Diagnose import errors

Open the QGIS Python Console (`Plugins → Python Console`) and run:

```python
import QSIGA_CAL
```

The console will show the exact missing module or error.

---

**Developed by:** The Nature Conservancy — miguel.canon@tnc.org
