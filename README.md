# BioSCART Minimal

BioSCART Minimal is a CT-only research prototype for biologically informed SCART/SFRT planning review. It imports DICOM CT and RTSTRUCT, rasterizes a selected GTV, creates transparent intratumor candidate regions, optionally extracts PyRadiomics features, and exports research artifacts that can be reviewed in 3D Slicer or a treatment planning system (TPS).

> **Research use only.** This is not a medical device, not treatment planning software, and must not be used for clinical treatment decisions.

---

## Quick start (for a new collaborator)

If you just received access to this private repo, here is the whole path from zero to a working run.

### 1. Prerequisites

| Tool | Version | Notes |
|------|---------|-------|
| Python | **3.9–3.11** | Use 3.9 if you want PyRadiomics (see note below). The core CT/RTSTRUCT workflow works on any 3.9+. |
| Git | any | To clone the repo. |
| (optional) 3D Slicer | any recent | To visually review the exported RTSTRUCT regions. |

Check what you have:

```powershell
python --version
git --version
```

### 2. Get the code

You should have received a GitHub collaborator invite by email — accept it first, then:

```powershell
git clone https://github.com/weisiyan-ops/bioscart-minimal.git
cd bioscart-minimal
```

### 3. Install

```powershell
python -m venv .venv
.\.venv\Scripts\pip.exe install -e .[dev,radiomics]
```

On macOS/Linux the activate/pip paths are `./.venv/bin/pip` instead.

**PyRadiomics note:** `pyradiomics` can be hard to build on Windows and on Python 3.12+. If the line above fails, install without it — the full CT/RTSTRUCT workflow still runs, only the radiomics feature extraction is skipped:

```powershell
.\.venv\Scripts\pip.exe install -e .[dev]
```

### 4. Verify the install

```powershell
.\.venv\Scripts\python.exe -m pytest tests -q -p no:cacheprovider
bioscart --help
```

If the tests pass and `--help` prints the options, you are ready to run a real case.

---

## What it does

1. Load a DICOM CT series and RTSTRUCT.
2. Select a GTV ROI by explicit name (`--gtv-name`) or simple GTV-name matching.
3. Convert the GTV contour to the CT image grid.
4. Generate five candidate BioSCART regions:
   - `BSCART_Rim_5mm` — inner rim shell (thickness set by `--rim-mm`)
   - `BSCART_Core`
   - `BSCART_CT_Low_Q25` — low-HU quartile
   - `BSCART_CT_High_Q75` — high-HU quartile
   - `BSCART_Texture_High_Q75` — high-texture quartile
5. Optionally extract PyRadiomics features.
6. Export JSON/CSV artifacts, RTSTRUCT regions, planning-objective templates, and (optionally) RTDOSE review metrics.

## Inputs

`--dicom-dir` should point to a folder containing **one CT series and its RTSTRUCT**. The RTSTRUCT must contain a GTV-like ROI. If the folder has multiple CT series, the largest is chosen automatically, or you can pin one with `--ct-series-uid`. A specific RTSTRUCT file can be forced with `--rtstruct`.

## Usage

Basic run:

```powershell
bioscart --dicom-dir C:\path\to\dicom --out-dir C:\path\to\output --gtv-name GTV
```

With optional TPS round-trip dose review:

```powershell
bioscart --dicom-dir C:\path\to\dicom --out-dir C:\path\to\output --gtv-name GTV --rtdose C:\path\to\RTDOSE.dcm
```

Useful options (`bioscart --help` for the full list):

| Option | Purpose |
|--------|---------|
| `--gtv-name` | Exact GTV ROI name (omit to auto-match a GTV-like name). |
| `--rim-mm` | Inner rim shell thickness in mm (default 5.0). |
| `--rtdose` | RTDOSE file for TPS round-trip dose evaluation. |
| `--ct-series-uid` | Pin a specific CT SeriesInstanceUID. |
| `--rtstruct` | Use a specific RTSTRUCT file. |
| `--pyradiomics-params` | Custom PyRadiomics YAML parameter file. |
| `--strict-tps-geometry` | Treat missing/mismatched DICOM geometry as a hard failure. |
| `--no-rtstruct-export` | Skip BioSCART RTSTRUCT export. |
| `--ice3-metrics` | ICE3/EDIC metrics JSON for combined review. |

## Outputs

Written to `--out-dir`:

- `bioscart_manifest.json`
- `region_summary.csv`
- `recommendations.json`
- `planning_objectives.json`
- `bioscart_regions_rtstruct.dcm`
- `dose_evaluation.json` — when `--rtdose` is provided
- `combined_review.json` — when ICE3/EDIC metrics are provided

To review the result, load the original CT plus `bioscart_regions_rtstruct.dcm` in 3D Slicer or your TPS.

## Development

```powershell
python -m pytest tests -q -p no:cacheprovider
```

## Data policy

Do not commit patient DICOM, RTSTRUCT, RTDOSE, NIfTI/NRRD volumes, generated output folders, screenshots with PHI, or local credential files. The `.gitignore` excludes common clinical-data and output patterns, but review every commit before pushing.

## Troubleshooting

- **`bioscart` not found** — the venv isn't active, or install didn't finish. Run via the full path: `.\.venv\Scripts\bioscart.exe --help`.
- **`pip install` fails on `pyradiomics`** — install with `.[dev]` only (see install note). Radiomics is optional.
- **"No GTV ROI found"** — pass the exact ROI name with `--gtv-name`, or check the RTSTRUCT ROI names.
- **Geometry / series mismatch errors** — confirm the CT and RTSTRUCT belong to the same study; use `--ct-series-uid` to disambiguate.

## Contact

Maintainer: Weisi Yan (weisiyan@gmail.com). For access requests or questions, open an issue or email the maintainer.
