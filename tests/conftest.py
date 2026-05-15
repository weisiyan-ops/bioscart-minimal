"""Shared fixtures for BioSCART Minimal tests. All synthetic — no real DICOM files."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pydicom
import pytest
import SimpleITK as sitk
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from bioscart_minimal.dicom_io import CTImage, GTVMask


SHAPE = (64, 128, 128)
SPACING_ZYX = (2.5, 0.9766, 0.9766)
ORIGIN_XYZ = (0.0, 0.0, 0.0)


def _sphere_mask(shape, center, radius_voxels):
    zz, yy, xx = np.ogrid[
        :shape[0], :shape[1], :shape[2]
    ]
    dist = np.sqrt(
        (zz - center[0]) ** 2 + (yy - center[1]) ** 2 + (xx - center[2]) ** 2
    )
    return dist <= radius_voxels


@pytest.fixture
def synthetic_ct_volume():
    vol = np.full(SHAPE, 30.0, dtype=np.float32)
    tumor = _sphere_mask(SHAPE, (32, 64, 64), 20)
    vol[tumor] += 30.0 * np.random.default_rng(42).random(tumor.sum()).astype(np.float32)
    return vol


@pytest.fixture
def synthetic_spacing_zyx():
    return SPACING_ZYX


@pytest.fixture
def synthetic_gtv_mask():
    return _sphere_mask(SHAPE, (32, 64, 64), 20)


@pytest.fixture
def small_gtv_mask():
    return _sphere_mask(SHAPE, (32, 64, 64), 3)


@pytest.fixture
def synthetic_ct_image(synthetic_ct_volume):
    img = sitk.GetImageFromArray(synthetic_ct_volume)
    col_sp, row_sp, slice_sp = SPACING_ZYX[2], SPACING_ZYX[1], SPACING_ZYX[0]
    img.SetSpacing((col_sp, row_sp, slice_sp))
    img.SetOrigin(ORIGIN_XYZ)
    return CTImage(
        array_hu=synthetic_ct_volume,
        sitk_image=img,
        slice_paths=[],
        series_instance_uid="1.2.3.4.5",
        origin_xyz_mm=ORIGIN_XYZ,
        spacing_zyx_mm=SPACING_ZYX,
        orientation=(1, 0, 0, 0, 1, 0),
        warnings=[],
    )


@pytest.fixture
def synthetic_gtv(synthetic_gtv_mask):
    voxel_cc = float(np.prod(SPACING_ZYX) / 1000.0)
    return GTVMask(
        roi_name="GTV",
        roi_number=1,
        mask=synthetic_gtv_mask,
        volume_cc=float(synthetic_gtv_mask.sum() * voxel_cc),
        warnings=[],
    )


@pytest.fixture
def small_gtv(small_gtv_mask):
    voxel_cc = float(np.prod(SPACING_ZYX) / 1000.0)
    return GTVMask(
        roi_name="GTV",
        roi_number=1,
        mask=small_gtv_mask,
        volume_cc=float(small_gtv_mask.sum() * voxel_cc),
        warnings=[],
    )


def _circle_contour_data(cx_mm, cy_mm, z_mm, radius_mm, n_points=36):
    """Generate flattened DICOM ContourData for a circle at given z."""
    angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
    xs = cx_mm + radius_mm * np.cos(angles)
    ys = cy_mm + radius_mm * np.sin(angles)
    zs = np.full_like(xs, z_mm)
    return list(np.column_stack([xs, ys, zs]).ravel())


def _build_rtstruct(roi_name, roi_number, contour_slices):
    """Build a minimal pydicom RTSTRUCT Dataset.

    contour_slices: list of (z_mm, list_of_contour_data_arrays)
    """
    ds = Dataset()
    ds.Modality = "RTSTRUCT"
    ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3"
    ds.SOPInstanceUID = generate_uid()

    roi_item = Dataset()
    roi_item.ROIName = roi_name
    roi_item.ROINumber = roi_number
    roi_item.ReferencedFrameOfReferenceUID = generate_uid()
    ds.StructureSetROISequence = Sequence([roi_item])

    contour_seq = []
    for _z_mm, contour_data_list in contour_slices:
        for cdata in contour_data_list:
            c_item = Dataset()
            c_item.ContourGeometricType = "CLOSED_PLANAR"
            c_item.NumberOfContourPoints = len(cdata) // 3
            c_item.ContourData = cdata
            contour_seq.append(c_item)

    roi_contour_item = Dataset()
    roi_contour_item.ReferencedROINumber = roi_number
    roi_contour_item.ContourSequence = Sequence(contour_seq)
    ds.ROIContourSequence = Sequence([roi_contour_item])

    return ds


@pytest.fixture
def mock_rtstruct_simple():
    """RTSTRUCT with a single circular GTV contour on 3 slices."""
    slices = []
    for z_mm in [78.0, 80.5, 83.0]:
        cdata = _circle_contour_data(62.5, 62.5, z_mm, radius_mm=15.0)
        slices.append((z_mm, [cdata]))
    return _build_rtstruct("GTV", 1, slices)


@pytest.fixture
def mock_rtstruct_with_hole():
    """RTSTRUCT with outer circle + inner hole on one slice."""
    outer = _circle_contour_data(62.5, 62.5, 80.0, radius_mm=20.0, n_points=72)
    inner = _circle_contour_data(62.5, 62.5, 80.0, radius_mm=5.0, n_points=36)
    return _build_rtstruct("GTV", 1, [(80.0, [outer, inner])])


@pytest.fixture
def mock_rtstruct_multi_roi():
    """RTSTRUCT with GTV, PTV, and CTV ROIs."""
    ds = Dataset()
    ds.Modality = "RTSTRUCT"
    ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3"
    ds.SOPInstanceUID = generate_uid()

    rois = []
    for name, number in [("PTV_Total", 1), ("CTV_60", 2), ("GTV", 3), ("BODY", 4)]:
        item = Dataset()
        item.ROIName = name
        item.ROINumber = number
        rois.append(item)
    ds.StructureSetROISequence = Sequence(rois)
    ds.ROIContourSequence = Sequence([])
    return ds
