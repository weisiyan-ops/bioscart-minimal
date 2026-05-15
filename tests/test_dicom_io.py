"""Tests for bioscart_minimal.dicom_io."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pydicom
import pytest
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

from bioscart_minimal.dicom_io import (
    CTImage,
    DICOMInventory,
    find_gtv_roi,
    rasterize_roi_to_ct_grid,
    scan_dicom_directory,
    select_largest_ct_series,
    sitk_mask_from_array,
)

from .conftest import ORIGIN_XYZ, SHAPE, SPACING_ZYX, _circle_contour_data


# ---------------------------------------------------------------------------
# scan_dicom_directory
# ---------------------------------------------------------------------------

class TestScanDICOMDirectory:

    def test_scan_empty_directory(self, tmp_path):
        inv = scan_dicom_directory(tmp_path)
        assert inv.ct_series == {}
        assert inv.rtstructs == []

    def test_scan_with_ct_file(self, tmp_path):
        ds = Dataset()
        ds.Modality = "CT"
        ds.SeriesInstanceUID = "1.2.3"
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        ds.SOPInstanceUID = generate_uid()
        ds.is_little_endian = True
        ds.is_implicit_VR = True
        path = tmp_path / "ct001.dcm"
        pydicom.dcmwrite(str(path), ds, write_like_original=False)

        inv = scan_dicom_directory(tmp_path)
        assert "1.2.3" in inv.ct_series
        assert len(inv.ct_series["1.2.3"]) == 1

    def test_scan_with_rtstruct_file(self, tmp_path):
        ds = Dataset()
        ds.Modality = "RTSTRUCT"
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3"
        ds.SOPInstanceUID = generate_uid()
        ds.is_little_endian = True
        ds.is_implicit_VR = True
        path = tmp_path / "rs001.dcm"
        pydicom.dcmwrite(str(path), ds, write_like_original=False)

        inv = scan_dicom_directory(tmp_path)
        assert len(inv.rtstructs) == 1


# ---------------------------------------------------------------------------
# select_largest_ct_series
# ---------------------------------------------------------------------------

class TestSelectLargestCTSeries:

    def test_selects_largest(self):
        inv = DICOMInventory(
            ct_series={
                "uid_a": [Path("a1.dcm"), Path("a2.dcm")],
                "uid_b": [Path("b1.dcm"), Path("b2.dcm"), Path("b3.dcm")],
            }
        )
        uid, paths = select_largest_ct_series(inv)
        assert uid == "uid_b"
        assert len(paths) == 3

    def test_raises_when_empty(self):
        inv = DICOMInventory()
        with pytest.raises(ValueError, match="No CT series"):
            select_largest_ct_series(inv)


# ---------------------------------------------------------------------------
# find_gtv_roi
# ---------------------------------------------------------------------------

class TestFindGTVROI:

    def test_find_exact_gtv(self, mock_rtstruct_multi_roi):
        name, number = find_gtv_roi(mock_rtstruct_multi_roi)
        assert name == "GTV"
        assert number == 3

    def test_find_preferred_name(self, mock_rtstruct_multi_roi):
        name, number = find_gtv_roi(mock_rtstruct_multi_roi, preferred_name="CTV_60")
        assert name == "CTV_60"
        assert number == 2

    def test_preferred_name_not_found_raises(self, mock_rtstruct_multi_roi):
        with pytest.raises(ValueError, match="not found"):
            find_gtv_roi(mock_rtstruct_multi_roi, preferred_name="NONEXISTENT")

    def test_regex_fallback(self):
        ds = Dataset()
        ds.Modality = "RTSTRUCT"
        item = Dataset()
        item.ROIName = "Gross Tumor Volume"
        item.ROINumber = 5
        ds.StructureSetROISequence = Sequence([item])
        name, number = find_gtv_roi(ds)
        assert name == "Gross Tumor Volume"
        assert number == 5

    def test_excludes_ptv(self):
        ds = Dataset()
        ds.Modality = "RTSTRUCT"
        rois = []
        for roi_name, roi_num in [("PTV_GTV_boost", 1), ("BODY", 2)]:
            item = Dataset()
            item.ROIName = roi_name
            item.ROINumber = roi_num
            rois.append(item)
        ds.StructureSetROISequence = Sequence(rois)
        with pytest.raises(ValueError, match="No GTV-like ROI"):
            find_gtv_roi(ds)

    def test_case_insensitive_match(self):
        ds = Dataset()
        ds.Modality = "RTSTRUCT"
        item = Dataset()
        item.ROIName = "gtv"
        item.ROINumber = 10
        ds.StructureSetROISequence = Sequence([item])
        name, number = find_gtv_roi(ds)
        assert number == 10

    def test_no_structure_set_raises(self):
        ds = Dataset()
        ds.Modality = "RTSTRUCT"
        with pytest.raises(ValueError, match="StructureSetROISequence"):
            find_gtv_roi(ds)


# ---------------------------------------------------------------------------
# rasterize_roi_to_ct_grid — XOR hole regression
# ---------------------------------------------------------------------------

class TestRasterizeROI:

    def _make_ct_for_rasterize(self, z_positions, tmp_path):
        """Build a minimal CTImage with real DICOM slice files for rasterization tests."""
        import SimpleITK as sitk

        nz = len(z_positions)
        shape = (nz, 128, 128)
        vol = np.zeros(shape, dtype=np.float32)
        img = sitk.GetImageFromArray(vol)
        col_sp = SPACING_ZYX[2]
        row_sp = SPACING_ZYX[1]
        slice_sp = float(np.median(np.diff(z_positions))) if nz > 1 else 2.5
        img.SetSpacing((col_sp, row_sp, slice_sp))
        img.SetOrigin(ORIGIN_XYZ)

        slice_paths = []
        for i, z in enumerate(z_positions):
            ds = Dataset()
            ds.Modality = "CT"
            ds.ImagePositionPatient = [0.0, 0.0, z]
            ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
            ds.SOPInstanceUID = generate_uid()
            ds.file_meta = pydicom.dataset.FileMetaDataset()
            ds.file_meta.MediaStorageSOPClassUID = ds.SOPClassUID
            ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
            ds.file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
            path = tmp_path / f"ct_slice_{i:03d}.dcm"
            pydicom.dcmwrite(str(path), ds, write_like_original=False)
            slice_paths.append(path)

        return CTImage(
            array_hu=vol,
            sitk_image=img,
            slice_paths=slice_paths,
            series_instance_uid="1.2.3",
            origin_xyz_mm=ORIGIN_XYZ,
            spacing_zyx_mm=(slice_sp, row_sp, col_sp),
            orientation=(1, 0, 0, 0, 1, 0),
            warnings=[],
        )

    def test_simple_contour_nonempty(self, mock_rtstruct_simple, tmp_path):
        ct = self._make_ct_for_rasterize([78.0, 80.5, 83.0], tmp_path)
        gtv = rasterize_roi_to_ct_grid(mock_rtstruct_simple, 1, "GTV", ct)
        assert np.any(gtv.mask)
        assert gtv.volume_cc > 0

    def test_xor_hole_handling(self, mock_rtstruct_with_hole, tmp_path):
        """Regression: inner contour must create a hole via XOR, not fill."""
        ct = self._make_ct_for_rasterize([80.0], tmp_path)
        gtv = rasterize_roi_to_ct_grid(mock_rtstruct_with_hole, 1, "GTV", ct)

        center_row = int(round(62.5 / SPACING_ZYX[1]))
        center_col = int(round(62.5 / SPACING_ZYX[2]))
        assert not gtv.mask[0, center_row, center_col], (
            "Center voxel should be False (hole from inner contour XOR)"
        )

        rim_row = int(round((62.5 + 15.0) / SPACING_ZYX[1]))
        rim_col = int(round(62.5 / SPACING_ZYX[2]))
        assert gtv.mask[0, rim_row, rim_col], (
            "Rim voxel should be True (outside inner hole, inside outer contour)"
        )

    def test_empty_contour_raises(self, tmp_path):
        ds = Dataset()
        ds.Modality = "RTSTRUCT"
        item = Dataset()
        item.ROIName = "GTV"
        item.ROINumber = 1
        ds.StructureSetROISequence = Sequence([item])
        ds.ROIContourSequence = Sequence([])

        ct = self._make_ct_for_rasterize([0.0], tmp_path)
        with pytest.raises(ValueError, match="no contours"):
            rasterize_roi_to_ct_grid(ds, 1, "GTV", ct)


# ---------------------------------------------------------------------------
# sitk_mask_from_array
# ---------------------------------------------------------------------------

class TestSITKMaskFromArray:

    def test_preserves_geometry(self, synthetic_ct_image, synthetic_gtv_mask):
        mask_img = sitk_mask_from_array(synthetic_gtv_mask, synthetic_ct_image.sitk_image)
        assert mask_img.GetSpacing() == synthetic_ct_image.sitk_image.GetSpacing()
        assert mask_img.GetOrigin() == synthetic_ct_image.sitk_image.GetOrigin()
        assert mask_img.GetSize() == synthetic_ct_image.sitk_image.GetSize()
