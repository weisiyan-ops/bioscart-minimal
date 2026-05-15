# BioSCART Minimal

BioSCART Minimal is a CT-only research prototype for biologically informed SCART/SFRT planning review. It imports DICOM CT and RTSTRUCT, rasterizes a selected GTV, creates transparent intratumor candidate regions, optionally extracts PyRadiomics features, and exports research artifacts that can be reviewed in 3D Slicer or a treatment planning system.

This repository is for research use only. It is not a medical device, not treatment planning software, and must not be used for clinical treatment decisions.

## Current Workflow

1. Load a DICOM CT series and RTSTRUCT.
2. Select a GTV ROI by explicit name or simple GTV-name matching.
3. Convert the GTV contour to the CT image grid.
4. Generate five candidate BioSCART regions:
   - `BSCART_Rim_5mm`
   - `BSCART_Core`
   - `BSCART_CT_Low_Q25`
   - `BSCART_CT_High_Q75`
   - `BSCART_Texture_High_Q75`
5. Optionally extract PyRadiomics features.
6. Export JSON/CSV artifacts, RTSTRUCT regions, planning-objective templates, and optional RTDOSE review metrics.

## Installation

```powershell
python -m venv .venv
.\.venv\Scripts\pip.exe install -e .[dev,radiomics]
```

If PyRadiomics is not installed, the core CT/RTSTRUCT workflow still runs, but radiomics extraction is skipped.

## Example

```powershell
bioscart --dicom-dir C:\path\to\dicom --out-dir C:\path\to\output --gtv-name GTV
```

With optional TPS round-trip dose review:

```powershell
bioscart --dicom-dir C:\path\to\dicom --out-dir C:\path\to\output --gtv-name GTV --rtdose C:\path\to\RTDOSE.dcm
```

## Outputs

Typical outputs include:

- `bioscart_manifest.json`
- `region_summary.csv`
- `recommendations.json`
- `planning_objectives.json`
- `bioscart_regions_rtstruct.dcm`
- `dose_evaluation.json` when RTDOSE is provided
- `combined_review.json` when ICE3/EDIC metrics are provided

## Development

```powershell
python -m pytest tests -q -p no:cacheprovider
```

## Data Policy

Do not commit patient DICOM, RTSTRUCT, RTDOSE, NIfTI/NRRD volumes, generated output folders, screenshots with PHI, or local credential files. The `.gitignore` is configured to exclude common clinical data and output patterns, but every commit should still be reviewed before pushing.
