"""DICOM RTDOSE export for BioSCART SFRT dose prescription visualization.

Creates a synthetic RTDOSE file where each voxel inside a BioSCART region
carries the prescribed SFRT dose for that region. The planner imports this
into Eclipse alongside the BioSCART RTSTRUCT to see the spatial dose intent
as a color wash overlay.

This is a REFERENCE dose for visual guidance only -- it is NOT an
optimization input, NOT a calculated dose, and must not be used for
treatment delivery or plan evaluation.

Research use only.
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Any

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

from .dicom_io import CTImage
from .errors import DoseEvaluationError
from .regions import RegionMask
from .sfrt_dose_logic import RegionDoseAssignment


RTDOSE_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.2"


def export_prescription_rtdose(
    assignments: list[RegionDoseAssignment],
    regions: dict[str, RegionMask],
    ct: CTImage,
    reference_rtstruct: Dataset,
    output_dir: Path,
    dose_per_fraction: bool = True,
) -> Path:
    """Export a synthetic RTDOSE with prescribed SFRT doses painted onto regions.

    Parameters
    ----------
    assignments : list[RegionDoseAssignment]
        Region-to-dose mapping from sfrt_dose_logic.assign_doses_to_regions().
    regions : dict[str, RegionMask]
        BioSCART region masks (must match assignments by region_name).
    ct : CTImage
        The planning CT (defines the grid geometry).
    reference_rtstruct : Dataset
        Reference RTSTRUCT for copying patient/study tags.
    output_dir : Path
        Where to write the RTDOSE file.
    dose_per_fraction : bool
        If True (default), paint dose_per_fraction_gy values.
        If False, paint total_dose_gy values.

    Returns
    -------
    Path to the exported RTDOSE file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    dose_array = _build_dose_array(assignments, regions, ct, dose_per_fraction)
    ds = _build_rtdose_dataset(dose_array, ct, reference_rtstruct, dose_per_fraction)

    path = output_dir / "RD.BioSCART_Prescription.dcm"
    ds.save_as(str(path), write_like_original=False)
    return path


def _build_dose_array(
    assignments: list[RegionDoseAssignment],
    regions: dict[str, RegionMask],
    ct: CTImage,
    dose_per_fraction: bool,
) -> np.ndarray:
    """Paint prescribed doses onto the CT grid.

    When regions overlap, the HIGHEST dose wins (peak > intermediate > valley).
    Voxels outside all regions get 0.
    """
    dose_array = np.zeros(ct.array_hu.shape, dtype=np.float32)

    level_priority = {"valley": 0, "intermediate": 1, "peak": 2}
    sorted_assignments = sorted(
        assignments,
        key=lambda a: level_priority.get(a.dose_level, 1),
    )

    for assignment in sorted_assignments:
        region = regions.get(assignment.region_name)
        if region is None:
            continue
        dose_val = assignment.dose_per_fraction_gy if dose_per_fraction else assignment.total_dose_gy
        dose_array[region.mask] = dose_val

    return dose_array


def _build_rtdose_dataset(
    dose_array: np.ndarray,
    ct: CTImage,
    reference_rtstruct: Dataset,
    dose_per_fraction: bool,
) -> FileDataset:
    """Build a DICOM RTDOSE dataset from a dose array on the CT grid."""
    now = _dt.datetime.now()

    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = RTDOSE_SOP_CLASS_UID
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    file_meta.ImplementationClassUID = generate_uid(prefix="1.2.826.0.1.3680043.10.543.")

    ds = FileDataset(None, {}, file_meta=file_meta, preamble=b"\0" * 128)
    ds.is_little_endian = True
    ds.is_implicit_VR = False

    _copy_patient_study_tags(reference_rtstruct, ds)

    ds.SOPClassUID = RTDOSE_SOP_CLASS_UID
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.Modality = "RTDOSE"
    ds.SeriesInstanceUID = generate_uid()
    ds.SeriesNumber = "9902"
    ds.InstanceNumber = "1"
    ds.Manufacturer = "BioSCART Research"
    ds.ManufacturerModelName = "BioSCART Minimal SFRT"
    ds.SoftwareVersions = "0.1.0"
    ds.ContentDate = now.strftime("%Y%m%d")
    ds.ContentTime = now.strftime("%H%M%S")
    ds.ContentLabel = "BSCART_SFRT_RX"
    ds.ContentDescription = (
        "BioSCART SFRT prescription reference - RESEARCH ONLY - "
        "NOT a calculated dose - NOT for treatment delivery"
    )

    frame_uid = ct.frame_of_reference_uid or generate_uid()
    ds.FrameOfReferenceUID = frame_uid

    nz, ny, nx = dose_array.shape
    sz, sy, sx = ct.spacing_zyx_mm
    ox, oy, oz = ct.origin_xyz_mm

    ds.Rows = ny
    ds.Columns = nx
    ds.NumberOfFrames = nz
    ds.PixelSpacing = [f"{sy:.6f}", f"{sx:.6f}"]
    ds.ImagePositionPatient = [f"{ox:.6f}", f"{oy:.6f}", f"{oz:.6f}"]
    ds.ImageOrientationPatient = ["1", "0", "0", "0", "1", "0"]

    offsets = [round(float(i * sz), 6) for i in range(nz)]
    ds.GridFrameOffsetVector = offsets

    ds.DoseUnits = "GY"
    ds.DoseType = "PHYSICAL"
    ds.DoseSummationType = "FRACTION" if dose_per_fraction else "PLAN"
    ds.DoseComment = (
        "BioSCART SFRT prescribed dose reference. "
        "This is a visual guide for planning, NOT a calculated or deliverable dose. "
        "Research use only."
    )

    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 32
    ds.BitsStored = 32
    ds.HighBit = 31
    ds.PixelRepresentation = 0

    max_dose = float(np.max(dose_array)) if np.any(dose_array > 0) else 1.0
    scaling = max_dose / 4294967295.0
    ds.DoseGridScaling = f"{scaling:.10e}"

    if scaling > 0:
        int_array = np.clip(dose_array / scaling, 0, 4294967295).astype(np.uint32)
    else:
        int_array = np.zeros_like(dose_array, dtype=np.uint32)

    ds.PixelData = int_array.tobytes()
    ds[0x7FE0, 0x0010].VR = "OW"

    return ds


def _copy_patient_study_tags(source: Dataset, target: Dataset) -> None:
    for tag in [
        "PatientName", "PatientID", "PatientBirthDate", "PatientSex",
        "StudyInstanceUID", "StudyID", "StudyDate", "StudyTime",
        "AccessionNumber", "ReferringPhysicianName",
    ]:
        if hasattr(source, tag):
            setattr(target, tag, getattr(source, tag))
    if not hasattr(target, "StudyInstanceUID"):
        target.StudyInstanceUID = generate_uid()
    if not hasattr(target, "PatientID"):
        target.PatientID = "BIO-SCART-RESEARCH"
