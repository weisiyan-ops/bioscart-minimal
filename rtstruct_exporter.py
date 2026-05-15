"""DICOM RTSTRUCT export for BioSCART Minimal derived regions."""

from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Any

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid
from skimage import measure

from .dicom_io import CTImage
from .errors import StructureExportError
from .regions import RegionMask


RTSTRUCT_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.481.3"


def export_bioscart_rtstruct(
    regions: dict[str, RegionMask],
    ct: CTImage,
    reference_rtstruct: Dataset,
    output_dir: Path,
    structure_set_label: str = "BSCART_MIN",
) -> Path:
    """Export BioSCART masks as a DICOM RTSTRUCT for TPS/Slicer import."""
    if not regions:
        raise StructureExportError("No regions were provided for RTSTRUCT export")

    output_dir.mkdir(parents=True, exist_ok=True)
    ds = _build_rtstruct_dataset(regions, ct, reference_rtstruct, structure_set_label)
    path = output_dir / "RS.BioSCART_Minimal.dcm"
    ds.save_as(str(path), write_like_original=False)
    return path


def _build_rtstruct_dataset(
    regions: dict[str, RegionMask],
    ct: CTImage,
    reference_rtstruct: Dataset,
    structure_set_label: str,
) -> FileDataset:
    now = _dt.datetime.now()
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = RTSTRUCT_SOP_CLASS_UID
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    file_meta.ImplementationClassUID = generate_uid(prefix="1.2.826.0.1.3680043.10.543.")

    ds = FileDataset(None, {}, file_meta=file_meta, preamble=b"\0" * 128)
    ds.is_little_endian = True
    ds.is_implicit_VR = False

    _copy_patient_study_tags(reference_rtstruct, ds)

    ds.SOPClassUID = RTSTRUCT_SOP_CLASS_UID
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.Modality = "RTSTRUCT"
    ds.SeriesInstanceUID = generate_uid()
    ds.SeriesNumber = "9901"
    ds.InstanceNumber = "1"
    ds.Manufacturer = "BioSCART Research"
    ds.ManufacturerModelName = "BioSCART Minimal"
    ds.SoftwareVersions = "0.1.0"
    ds.StructureSetLabel = structure_set_label[:16]
    ds.StructureSetName = "BioSCART Minimal Research Structures"
    ds.StructureSetDescription = "Research only - not for treatment decision-making"
    ds.StructureSetDate = now.strftime("%Y%m%d")
    ds.StructureSetTime = now.strftime("%H%M%S")
    ds.ContentDate = ds.StructureSetDate
    ds.ContentTime = ds.StructureSetTime

    frame_uid = (
        ct.frame_of_reference_uid
        or _reference_rtstruct_frame_uid(reference_rtstruct)
        or generate_uid()
    )
    ds.ReferencedFrameOfReferenceSequence = Sequence([_referenced_frame_item(frame_uid)])

    roi_seq = []
    roi_contour_seq = []
    roi_obs_seq = []

    for roi_number, region in enumerate(regions.values(), start=1):
        roi_seq.append(_structure_set_roi_item(roi_number, region, frame_uid))
        roi_contour_seq.append(_roi_contour_item(roi_number, region, ct))
        roi_obs_seq.append(_roi_observation_item(roi_number, region))

    ds.StructureSetROISequence = Sequence(roi_seq)
    ds.ROIContourSequence = Sequence(roi_contour_seq)
    ds.RTROIObservationsSequence = Sequence(roi_obs_seq)
    return ds


def _copy_patient_study_tags(source: Dataset, target: Dataset) -> None:
    for tag in [
        "PatientName",
        "PatientID",
        "PatientBirthDate",
        "PatientSex",
        "StudyInstanceUID",
        "StudyID",
        "StudyDate",
        "StudyTime",
        "AccessionNumber",
        "ReferringPhysicianName",
    ]:
        if hasattr(source, tag):
            setattr(target, tag, getattr(source, tag))
    if not hasattr(target, "StudyInstanceUID"):
        target.StudyInstanceUID = generate_uid()
    if not hasattr(target, "PatientID"):
        target.PatientID = "BIO-SCART-RESEARCH"


def _referenced_frame_item(frame_uid: str) -> Dataset:
    item = Dataset()
    item.FrameOfReferenceUID = frame_uid
    return item


def _structure_set_roi_item(roi_number: int, region: RegionMask, frame_uid: str) -> Dataset:
    item = Dataset()
    item.ROINumber = roi_number
    item.ReferencedFrameOfReferenceUID = frame_uid
    item.ROIName = region.name[:64]
    item.ROIDescription = f"Research only - BioSCART Minimal - {region.role}"[:300]
    item.ROIGenerationAlgorithm = "AUTOMATIC"
    return item


def _roi_contour_item(roi_number: int, region: RegionMask, ct: CTImage) -> Dataset:
    item = Dataset()
    item.ReferencedROINumber = roi_number
    item.ROIDisplayColor = _roi_color(region.name)
    contours = _mask_to_contour_items(region.mask, ct)
    item.ContourSequence = Sequence(contours)
    return item


def _roi_observation_item(roi_number: int, region: RegionMask) -> Dataset:
    item = Dataset()
    item.ObservationNumber = roi_number
    item.ReferencedROINumber = roi_number
    item.RTROIInterpretedType = "CONTROL"
    item.ROIInterpreter = ""
    item.ROIObservationDescription = f"Research only: {region.role}"[:300]
    return item


def _mask_to_contour_items(mask: np.ndarray, ct: CTImage) -> list[Dataset]:
    slice_spacing, row_spacing, col_spacing = ct.spacing_zyx_mm
    origin_x, origin_y, origin_z = ct.origin_xyz_mm
    contour_items: list[Dataset] = []

    for z_idx in range(mask.shape[0]):
        slice_mask = mask[z_idx].astype(np.uint8)
        if not np.any(slice_mask):
            continue
        contours = measure.find_contours(slice_mask, 0.5)
        for contour in contours:
            if contour.shape[0] < 3:
                continue
            rows = contour[:, 0]
            cols = contour[:, 1]
            z = origin_z + z_idx * slice_spacing
            points: list[float] = []
            for row, col in zip(rows, cols):
                points.extend([
                    round(origin_x + float(col) * col_spacing, 4),
                    round(origin_y + float(row) * row_spacing, 4),
                    round(z, 4),
                ])
            if len(points) < 9:
                continue
            item = Dataset()
            item.ContourGeometricType = "CLOSED_PLANAR"
            item.NumberOfContourPoints = len(points) // 3
            item.ContourData = points
            contour_items.append(item)

    if not contour_items:
        raise StructureExportError("A BioSCART region produced no contours for RTSTRUCT export")
    return contour_items


def _reference_rtstruct_frame_uid(ds: Dataset) -> str | None:
    if hasattr(ds, "ReferencedFrameOfReferenceSequence") and ds.ReferencedFrameOfReferenceSequence:
        uid = getattr(ds.ReferencedFrameOfReferenceSequence[0], "FrameOfReferenceUID", None)
        if uid:
            return str(uid)
    if hasattr(ds, "StructureSetROISequence"):
        for roi in ds.StructureSetROISequence:
            uid = getattr(roi, "ReferencedFrameOfReferenceUID", None)
            if uid:
                return str(uid)
    return None


def _roi_color(name: str) -> list[int]:
    colors: dict[str, list[int]] = {
        "BSCART_Rim_5mm": [0, 220, 220],
        "BSCART_Core": [255, 80, 80],
        "BSCART_CT_Low_Q25": [80, 120, 255],
        "BSCART_CT_High_Q75": [255, 180, 0],
        "BSCART_Texture_High_Q75": [180, 80, 255],
        "BSCART_Uncertain": [150, 150, 150],
    }
    return colors.get(name, [128, 128, 128])

