"""PyRadiomics integration for BioSCART minimal."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import SimpleITK as sitk

from .dicom_io import sitk_mask_from_array
from .regions import RegionMask


@dataclass
class RadiomicsRunResult:
    available: bool
    features_csv: str | None
    status_json: str
    message: str


def run_pyradiomics(
    image: sitk.Image,
    gtv_mask,
    regions: dict[str, RegionMask],
    output_dir: Path,
    params_path: Path | None = None,
) -> RadiomicsRunResult:
    """Run PyRadiomics for GTV and the five BioSCART regions if installed."""
    output_dir.mkdir(parents=True, exist_ok=True)
    status_path = output_dir / "radiomics_status.json"

    try:
        from radiomics import featureextractor  # type: ignore
    except ImportError:
        message = "PyRadiomics is not installed. Install with: pip install pyradiomics"
        status = {
            "available": False,
            "message": message,
            "features_csv": None,
        }
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        return RadiomicsRunResult(False, None, str(status_path), message)

    image_path = output_dir / "ct_image.nii.gz"
    sitk.WriteImage(image, str(image_path))

    extractor = (
        featureextractor.RadiomicsFeatureExtractor(str(params_path))
        if params_path
        else featureextractor.RadiomicsFeatureExtractor()
    )

    masks = {"GTV": gtv_mask}
    masks.update({name: region.mask for name, region in regions.items()})

    rows: list[dict[str, object]] = []
    for name, mask in masks.items():
        mask_img = sitk_mask_from_array(mask, image)
        mask_path = output_dir / f"{name}.nii.gz"
        sitk.WriteImage(mask_img, str(mask_path))

        result = extractor.execute(str(image_path), str(mask_path))
        row: dict[str, object] = {"region": name}
        for key, value in result.items():
            if key.startswith("diagnostics_"):
                continue
            try:
                row[key] = float(value)
            except Exception:
                row[key] = str(value)
        rows.append(row)

    csv_path = output_dir / "pyradiomics_features.csv"
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    status = {
        "available": True,
        "message": "PyRadiomics completed",
        "features_csv": str(csv_path),
        "params_path": str(params_path) if params_path else None,
        "regions": list(masks.keys()),
    }
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    return RadiomicsRunResult(True, str(csv_path), str(status_path), "PyRadiomics completed")

