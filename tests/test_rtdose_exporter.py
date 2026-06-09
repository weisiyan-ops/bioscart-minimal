"""Tests for the RTDOSE prescription exporter."""

import numpy as np
import pydicom
import pytest

from bioscart_minimal.dicom_io import CTImage
from bioscart_minimal.regions import RegionMask
from bioscart_minimal.rtdose_exporter import export_prescription_rtdose, _build_dose_array
from bioscart_minimal.sfrt_dose_logic import (
    RegionDoseAssignment,
    assign_doses_to_regions,
    scart_lung_protocol,
)


@pytest.fixture
def mock_ct():
    import SimpleITK as sitk
    from pathlib import Path
    array = np.zeros((10, 20, 20), dtype=np.float32)
    img = sitk.GetImageFromArray(array)
    img.SetSpacing((1.0, 1.0, 2.5))
    img.SetOrigin((-100.0, -100.0, -50.0))
    return CTImage(
        array_hu=array,
        sitk_image=img,
        slice_paths=[Path("dummy.dcm")],
        series_instance_uid="1.2.3.4.5.6.8",
        origin_xyz_mm=(-100.0, -100.0, -50.0),
        spacing_zyx_mm=(2.5, 1.0, 1.0),
        orientation=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        frame_of_reference_uid="1.2.3.4.5.6.7",
        warnings=[],
    )


@pytest.fixture
def mock_regions():
    mask_core = np.zeros((10, 20, 20), dtype=bool)
    mask_core[3:7, 8:12, 8:12] = True

    mask_rim = np.zeros((10, 20, 20), dtype=bool)
    mask_rim[2:8, 7:13, 7:13] = True
    mask_rim[3:7, 8:12, 8:12] = False

    mask_texture = np.zeros((10, 20, 20), dtype=bool)
    mask_texture[4:6, 9:11, 9:11] = True

    return {
        "BSCART_Core": RegionMask(
            "BSCART_Core", mask_core, "peak_candidate",
            "geometric_ct_only", 5.0, ["core"],
        ),
        "BSCART_Rim_5mm": RegionMask(
            "BSCART_Rim_5mm", mask_rim, "border_or_valley_candidate",
            "geometric_ct_only", 8.0, ["rim"],
        ),
        "BSCART_Texture_High_Q75": RegionMask(
            "BSCART_Texture_High_Q75", mask_texture, "heterogeneity_candidate",
            "ct_texture_surrogate", 1.0, ["texture"],
        ),
    }


@pytest.fixture
def mock_rtstruct():
    ds = pydicom.Dataset()
    ds.PatientName = "TEST"
    ds.PatientID = "BSCART001"
    ds.StudyInstanceUID = "1.2.3.4.5"
    ds.StudyID = "1"
    return ds


class TestBuildDoseArray:
    def test_peak_higher_than_valley(self, mock_regions, mock_ct):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        dose = _build_dose_array(assignments, mock_regions, mock_ct, dose_per_fraction=True)

        core_mask = mock_regions["BSCART_Core"].mask
        rim_mask = mock_regions["BSCART_Rim_5mm"].mask
        core_dose = dose[core_mask].mean()
        rim_dose = dose[rim_mask].mean()
        assert core_dose > rim_dose

    def test_outside_regions_is_zero(self, mock_regions, mock_ct):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        dose = _build_dose_array(assignments, mock_regions, mock_ct, dose_per_fraction=True)

        all_masks = np.zeros_like(dose, dtype=bool)
        for r in mock_regions.values():
            all_masks |= r.mask
        assert np.all(dose[~all_masks] == 0.0)

    def test_overlap_peak_wins(self, mock_regions, mock_ct):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        dose = _build_dose_array(assignments, mock_regions, mock_ct, dose_per_fraction=True)

        texture_mask = mock_regions["BSCART_Texture_High_Q75"].mask
        core_mask = mock_regions["BSCART_Core"].mask
        overlap = texture_mask & core_mask
        if np.any(overlap):
            peak_dose = proto.levels["peak"].dose_per_fraction_gy
            assert np.allclose(dose[overlap], peak_dose)

    def test_dose_per_fraction_vs_total(self, mock_regions, mock_ct):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)

        dose_fx = _build_dose_array(assignments, mock_regions, mock_ct, dose_per_fraction=True)
        dose_total = _build_dose_array(assignments, mock_regions, mock_ct, dose_per_fraction=False)

        core_mask = mock_regions["BSCART_Core"].mask
        fx_val = dose_fx[core_mask].mean()
        total_val = dose_total[core_mask].mean()
        assert total_val > fx_val


class TestExportPrescriptionRTDOSE:
    def test_exports_file(self, mock_regions, mock_ct, mock_rtstruct, tmp_path):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        path = export_prescription_rtdose(
            assignments, mock_regions, mock_ct, mock_rtstruct, tmp_path,
        )
        assert path.exists()
        assert path.name == "RD.BioSCART_Prescription.dcm"

    def test_dicom_readable(self, mock_regions, mock_ct, mock_rtstruct, tmp_path):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        path = export_prescription_rtdose(
            assignments, mock_regions, mock_ct, mock_rtstruct, tmp_path,
        )
        ds = pydicom.dcmread(str(path))
        assert ds.Modality == "RTDOSE"
        assert ds.DoseUnits == "GY"
        assert ds.Rows == 20
        assert ds.Columns == 20
        assert ds.NumberOfFrames == 10

    def test_dose_grid_scaling_nonzero(self, mock_regions, mock_ct, mock_rtstruct, tmp_path):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        path = export_prescription_rtdose(
            assignments, mock_regions, mock_ct, mock_rtstruct, tmp_path,
        )
        ds = pydicom.dcmread(str(path))
        scaling = float(ds.DoseGridScaling)
        assert scaling > 0

    def test_pixel_data_has_dose(self, mock_regions, mock_ct, mock_rtstruct, tmp_path):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        path = export_prescription_rtdose(
            assignments, mock_regions, mock_ct, mock_rtstruct, tmp_path,
        )
        ds = pydicom.dcmread(str(path))
        scaling = float(ds.DoseGridScaling)
        array = ds.pixel_array.astype(np.float64) * scaling
        assert np.max(array) > 0

    def test_research_only_label(self, mock_regions, mock_ct, mock_rtstruct, tmp_path):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        path = export_prescription_rtdose(
            assignments, mock_regions, mock_ct, mock_rtstruct, tmp_path,
        )
        ds = pydicom.dcmread(str(path))
        assert "RESEARCH ONLY" in ds.ContentDescription
        assert "NOT" in ds.DoseComment

    def test_frame_of_reference_matches_ct(self, mock_regions, mock_ct, mock_rtstruct, tmp_path):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        path = export_prescription_rtdose(
            assignments, mock_regions, mock_ct, mock_rtstruct, tmp_path,
        )
        ds = pydicom.dcmread(str(path))
        assert ds.FrameOfReferenceUID == mock_ct.frame_of_reference_uid

    def test_geometry_matches_ct(self, mock_regions, mock_ct, mock_rtstruct, tmp_path):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        path = export_prescription_rtdose(
            assignments, mock_regions, mock_ct, mock_rtstruct, tmp_path,
        )
        ds = pydicom.dcmread(str(path))
        ipp = [float(v) for v in ds.ImagePositionPatient]
        assert ipp == list(mock_ct.origin_xyz_mm)

    def test_total_dose_mode(self, mock_regions, mock_ct, mock_rtstruct, tmp_path):
        proto = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, proto)
        path = export_prescription_rtdose(
            assignments, mock_regions, mock_ct, mock_rtstruct, tmp_path,
            dose_per_fraction=False,
        )
        ds = pydicom.dcmread(str(path))
        assert ds.DoseSummationType == "PLAN"
