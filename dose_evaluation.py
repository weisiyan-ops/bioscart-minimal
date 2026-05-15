"""RTDOSE round-trip evaluation for BioSCART Minimal."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pydicom
import SimpleITK as sitk

from .dicom_io import CTImage, sitk_mask_from_array
from .errors import DoseEvaluationError
from .regions import RegionMask


@dataclass
class RTDoseGrid:
    """Loaded RTDOSE grid in Gy."""

    array_gy: np.ndarray  # z, y, x
    sitk_image: sitk.Image
    path: Path
    dose_units: str
    frame_of_reference_uid: str | None


def load_rtdose(path: Path | str) -> RTDoseGrid:
    """Load DICOM RTDOSE into a SimpleITK image and NumPy array in Gy."""
    path = Path(path)
    ds = pydicom.dcmread(path, force=True)
    if str(getattr(ds, "Modality", "")).upper() != "RTDOSE":
        raise DoseEvaluationError(f"Expected RTDOSE, got {getattr(ds, 'Modality', None)}")

    scaling = float(getattr(ds, "DoseGridScaling", 1.0))
    array = ds.pixel_array.astype(np.float32) * scaling
    if array.ndim == 2:
        array = array[np.newaxis, :, :]

    dose_units = str(getattr(ds, "DoseUnits", "GY")).upper()
    if dose_units != "GY":
        raise DoseEvaluationError(f"Unsupported RTDOSE units: {dose_units}")

    row_spacing, col_spacing = [float(v) for v in getattr(ds, "PixelSpacing")]
    z_spacing = _dose_z_spacing(ds)
    image = sitk.GetImageFromArray(array)
    image.SetSpacing((col_spacing, row_spacing, z_spacing))
    image.SetOrigin(tuple(float(v) for v in getattr(ds, "ImagePositionPatient", [0, 0, 0])))

    frame_uid = str(getattr(ds, "FrameOfReferenceUID", "")) or None
    return RTDoseGrid(
        array_gy=array,
        sitk_image=image,
        path=path,
        dose_units=dose_units,
        frame_of_reference_uid=frame_uid,
    )


def evaluate_regions_on_rtdose(
    rtdose_path: Path | str,
    ct: CTImage,
    gtv_mask: np.ndarray,
    regions: dict[str, RegionMask],
    output_dir: Path | None = None,
) -> dict[str, Any]:
    """Compute dose statistics for GTV and BioSCART regions."""
    dose = load_rtdose(rtdose_path)

    masks: dict[str, np.ndarray] = {"GTV": gtv_mask}
    masks.update({name: region.mask for name, region in regions.items()})

    stats = []
    for name, mask in masks.items():
        dose_mask = _resample_mask_to_dose(mask, ct, dose)
        stats.append(_dose_stats(name, dose.array_gy, dose_mask))

    report = {
        "schema_version": "bioscart.dose_eval.v0.1",
        "clinical_status": "research_only_not_for_treatment_decisions",
        "rtdose_path": str(rtdose_path),
        "dose_units": dose.dose_units,
        "ct_frame_of_reference_uid": ct.frame_of_reference_uid,
        "dose_frame_of_reference_uid": dose.frame_of_reference_uid,
        "frame_of_reference_match": (
            ct.frame_of_reference_uid == dose.frame_of_reference_uid
            if ct.frame_of_reference_uid and dose.frame_of_reference_uid
            else None
        ),
        "region_dose_stats": stats,
    }

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        path = output_dir / "bioscart_dose_evaluation.json"
        path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        report["output_path"] = str(path)

    return report


def _resample_mask_to_dose(mask: np.ndarray, ct: CTImage, dose: RTDoseGrid) -> np.ndarray:
    mask_img = sitk_mask_from_array(mask, ct.sitk_image)
    same_grid = (
        mask_img.GetSize() == dose.sitk_image.GetSize()
        and np.allclose(mask_img.GetSpacing(), dose.sitk_image.GetSpacing(), atol=1e-4)
        and np.allclose(mask_img.GetOrigin(), dose.sitk_image.GetOrigin(), atol=1e-4)
    )
    if same_grid:
        return mask.astype(bool)

    resampled = sitk.Resample(
        mask_img,
        dose.sitk_image,
        sitk.Transform(),
        sitk.sitkNearestNeighbor,
        0,
        sitk.sitkUInt8,
    )
    return sitk.GetArrayFromImage(resampled).astype(bool)


def _dose_stats(name: str, dose_gy: np.ndarray, mask: np.ndarray) -> dict[str, Any]:
    values = dose_gy[mask]
    if values.size == 0:
        return {
            "region": name,
            "voxel_count": 0,
            "mean_gy": None,
            "min_gy": None,
            "max_gy": None,
            "d95_gy": None,
            "d50_gy": None,
            "d5_gy": None,
        }
    return {
        "region": name,
        "voxel_count": int(values.size),
        "mean_gy": round(float(np.mean(values)), 4),
        "min_gy": round(float(np.min(values)), 4),
        "max_gy": round(float(np.max(values)), 4),
        "d95_gy": round(float(np.percentile(values, 5)), 4),
        "d50_gy": round(float(np.percentile(values, 50)), 4),
        "d5_gy": round(float(np.percentile(values, 95)), 4),
    }


def _dose_z_spacing(ds) -> float:
    if hasattr(ds, "GridFrameOffsetVector") and len(ds.GridFrameOffsetVector) > 1:
        offsets = np.array([float(v) for v in ds.GridFrameOffsetVector], dtype=float)
        diffs = np.diff(offsets)
        diffs = diffs[np.abs(diffs) > 1e-6]
        if diffs.size:
            return float(np.median(np.abs(diffs)))
    return float(getattr(ds, "SliceThickness", 1.0))

