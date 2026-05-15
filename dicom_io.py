"""DICOM CT and RTSTRUCT loading for the minimal BioSCART prototype.

Scope is intentionally narrow for v0:
- CT image series
- one RTSTRUCT
- one GTV-like ROI
- axial image geometry without gantry tilt
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pydicom
import SimpleITK as sitk
from pydicom.dataset import Dataset
from skimage import draw


GTV_INCLUDE_RE = re.compile(r"(^|[^A-Z0-9])(GTV|GROSS|TUMOU?R|LESION|PRIMARY)([^A-Z0-9]|$)", re.I)
GTV_EXCLUDE_RE = re.compile(r"(PTV|CTV|ITV|PRV|BOOST|RING|BODY|BOLUS)", re.I)


@dataclass
class DICOMInventory:
    ct_series: dict[str, list[Path]] = field(default_factory=dict)
    rtstructs: list[Path] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)


@dataclass
class CTImage:
    array_hu: np.ndarray  # z, y, x
    sitk_image: sitk.Image
    slice_paths: list[Path]
    series_instance_uid: str
    origin_xyz_mm: tuple[float, float, float]
    spacing_zyx_mm: tuple[float, float, float]
    orientation: tuple[float, ...]
    frame_of_reference_uid: str | None = None
    warnings: list[str] = field(default_factory=list)


@dataclass
class GTVMask:
    roi_name: str
    roi_number: int
    mask: np.ndarray  # z, y, x boolean
    volume_cc: float
    warnings: list[str] = field(default_factory=list)


def scan_dicom_directory(dicom_dir: Path | str) -> DICOMInventory:
    """Scan a directory recursively for CT and RTSTRUCT files."""
    dicom_dir = Path(dicom_dir)
    inventory = DICOMInventory()

    for path in dicom_dir.rglob("*"):
        if not path.is_file():
            continue
        try:
            ds = pydicom.dcmread(path, stop_before_pixels=True, force=True)
        except Exception:
            inventory.skipped.append(str(path))
            continue

        modality = str(getattr(ds, "Modality", "")).upper()
        if modality == "CT":
            series_uid = str(getattr(ds, "SeriesInstanceUID", "UNKNOWN_SERIES"))
            inventory.ct_series.setdefault(series_uid, []).append(path)
        elif modality == "RTSTRUCT":
            inventory.rtstructs.append(path)

    for paths in inventory.ct_series.values():
        paths.sort()
    inventory.rtstructs.sort()
    return inventory


def select_largest_ct_series(inventory: DICOMInventory) -> tuple[str, list[Path]]:
    """Select the CT series with the most slices."""
    if not inventory.ct_series:
        raise ValueError("No CT series found in DICOM directory")
    series_uid, paths = max(inventory.ct_series.items(), key=lambda item: len(item[1]))
    return series_uid, paths


def load_ct_series(paths: Iterable[Path | str]) -> CTImage:
    """Load an axial CT series into HU array and SimpleITK image."""
    slice_items: list[tuple[float, int, Path, Dataset]] = []
    warnings: list[str] = []

    for raw_path in paths:
        path = Path(raw_path)
        ds = pydicom.dcmread(path, force=True)
        if str(getattr(ds, "Modality", "")).upper() != "CT":
            continue

        ipp = [float(v) for v in getattr(ds, "ImagePositionPatient", [0, 0, 0])]
        instance = int(getattr(ds, "InstanceNumber", 0))
        slice_items.append((ipp[2], instance, path, ds))

    if not slice_items:
        raise ValueError("No readable CT slices in selected series")

    slice_items.sort(key=lambda item: (item[0], item[1]))
    datasets = [item[3] for item in slice_items]
    slice_paths = [item[2] for item in slice_items]

    first = datasets[0]
    pixel_spacing = [float(v) for v in first.PixelSpacing]  # row, col
    row_spacing, col_spacing = pixel_spacing

    z_positions = np.array([item[0] for item in slice_items], dtype=float)
    if len(z_positions) > 1:
        z_diffs = np.diff(np.sort(z_positions))
        slice_spacing = float(np.median(np.abs(z_diffs[z_diffs != 0]))) if np.any(z_diffs != 0) else float(getattr(first, "SliceThickness", 1.0))
    else:
        slice_spacing = float(getattr(first, "SliceThickness", 1.0))

    arrays = []
    for ds in datasets:
        slope = float(getattr(ds, "RescaleSlope", 1.0))
        intercept = float(getattr(ds, "RescaleIntercept", 0.0))
        arrays.append(ds.pixel_array.astype(np.float32) * slope + intercept)

    volume = np.stack(arrays, axis=0)
    image = sitk.GetImageFromArray(volume)
    image.SetSpacing((col_spacing, row_spacing, slice_spacing))

    origin = tuple(float(v) for v in getattr(first, "ImagePositionPatient", [0, 0, 0]))
    image.SetOrigin(origin)

    orientation = tuple(float(v) for v in getattr(first, "ImageOrientationPatient", [1, 0, 0, 0, 1, 0]))
    if not _is_identity_axial(orientation):
        warnings.append(
            "Non-identity CT orientation detected. Minimal v0 assumes axial patient-coordinate alignment; visual QA is required."
        )

    frame_uid = str(getattr(first, "FrameOfReferenceUID", "")) or None
    frame_uids = {str(getattr(ds, "FrameOfReferenceUID", "")) for ds in datasets}
    frame_uids.discard("")
    if len(frame_uids) > 1:
        warnings.append("CT slices contain multiple FrameOfReferenceUID values; visual QA is required.")

    series_uid = str(getattr(first, "SeriesInstanceUID", "UNKNOWN_SERIES"))
    return CTImage(
        array_hu=volume,
        sitk_image=image,
        slice_paths=slice_paths,
        series_instance_uid=series_uid,
        origin_xyz_mm=origin,
        spacing_zyx_mm=(slice_spacing, row_spacing, col_spacing),
        orientation=orientation,
        frame_of_reference_uid=frame_uid,
        warnings=warnings,
    )


def load_rtstruct(path: Path | str) -> Dataset:
    """Load an RTSTRUCT dataset."""
    ds = pydicom.dcmread(Path(path), force=True)
    if str(getattr(ds, "Modality", "")).upper() != "RTSTRUCT":
        raise ValueError(f"Expected RTSTRUCT, got {getattr(ds, 'Modality', None)}")
    return ds


def find_gtv_roi(rtstruct: Dataset, preferred_name: str | None = None) -> tuple[str, int]:
    """Find a GTV-like ROI in an RTSTRUCT."""
    if not hasattr(rtstruct, "StructureSetROISequence"):
        raise ValueError("RTSTRUCT has no StructureSetROISequence")

    rois: list[tuple[str, int]] = []
    for item in rtstruct.StructureSetROISequence:
        name = str(getattr(item, "ROIName", "")).strip()
        number = int(getattr(item, "ROINumber"))
        rois.append((name, number))

    if preferred_name:
        for name, number in rois:
            if name.lower() == preferred_name.lower():
                return name, number
        raise ValueError(f"Requested GTV ROI '{preferred_name}' was not found")

    exact = [(name, number) for name, number in rois if name.upper() == "GTV"]
    if exact:
        return exact[0]

    candidates = [
        (name, number)
        for name, number in rois
        if GTV_INCLUDE_RE.search(name) and not GTV_EXCLUDE_RE.search(name)
    ]
    if not candidates:
        roi_list = ", ".join(name for name, _ in rois[:40])
        raise ValueError(f"No GTV-like ROI found. Available ROI names include: {roi_list}")

    candidates.sort(key=lambda item: (len(item[0]), item[0].lower()))
    return candidates[0]


def rasterize_roi_to_ct_grid(rtstruct: Dataset, roi_number: int, roi_name: str, ct: CTImage) -> GTVMask:
    """Rasterize RTSTRUCT contours for one ROI onto the CT image grid."""
    mask = np.zeros(ct.array_hu.shape, dtype=bool)
    warnings: list[str] = []

    contour_items = []
    for roi_contour in getattr(rtstruct, "ROIContourSequence", []):
        if int(getattr(roi_contour, "ReferencedROINumber", -1)) == roi_number:
            contour_items = list(getattr(roi_contour, "ContourSequence", []))
            break

    if not contour_items:
        raise ValueError(f"ROI '{roi_name}' has no contours")

    z_positions = np.array([
        float(pydicom.dcmread(path, stop_before_pixels=True, force=True).ImagePositionPatient[2])
        for path in ct.slice_paths
    ])

    origin_x, origin_y, _origin_z = ct.origin_xyz_mm
    slice_spacing, row_spacing, col_spacing = ct.spacing_zyx_mm

    slices: dict[int, list[tuple[np.ndarray, np.ndarray]]] = defaultdict(list)
    for contour in contour_items:
        data = list(getattr(contour, "ContourData", []))
        if len(data) < 9:
            continue
        points = np.array(data, dtype=float).reshape(-1, 3)
        contour_z = float(np.median(points[:, 2]))
        z_idx = int(np.argmin(np.abs(z_positions - contour_z)))
        z_error = abs(float(z_positions[z_idx]) - contour_z)
        if z_error > max(slice_spacing, 1.0):
            warnings.append(
                f"Skipped contour at z={contour_z:.2f} mm; nearest CT slice error {z_error:.2f} mm"
            )
            continue

        cols = (points[:, 0] - origin_x) / col_spacing
        rows = (points[:, 1] - origin_y) / row_spacing
        rr, cc_arr = draw.polygon(rows, cols, shape=mask.shape[1:])
        slices[z_idx].append((rr, cc_arr))

    for z_idx, polygons in slices.items():
        slice_mask = np.zeros(mask.shape[1:], dtype=bool)
        for rr, cc_arr in polygons:
            contour_fill = np.zeros_like(slice_mask)
            contour_fill[rr, cc_arr] = True
            slice_mask = np.logical_xor(slice_mask, contour_fill)
        mask[z_idx] = slice_mask

    if not np.any(mask):
        raise ValueError(f"ROI '{roi_name}' rasterized to an empty mask")

    voxel_cc = (slice_spacing * row_spacing * col_spacing) / 1000.0
    volume_cc = float(mask.sum() * voxel_cc)
    return GTVMask(
        roi_name=roi_name,
        roi_number=roi_number,
        mask=mask,
        volume_cc=volume_cc,
        warnings=warnings,
    )


def sitk_mask_from_array(mask: np.ndarray, reference: sitk.Image) -> sitk.Image:
    """Create a SimpleITK UInt8 mask image with reference geometry."""
    mask_img = sitk.GetImageFromArray(mask.astype(np.uint8))
    mask_img.CopyInformation(reference)
    return mask_img


def _is_identity_axial(orientation: tuple[float, ...], tolerance: float = 1e-3) -> bool:
    expected = np.array([1, 0, 0, 0, 1, 0], dtype=float)
    observed = np.array(orientation, dtype=float)
    return observed.shape == expected.shape and np.allclose(observed, expected, atol=tolerance)
