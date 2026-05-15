"""Tests for RTDOSE region evaluation."""

from __future__ import annotations

import numpy as np
import pydicom
import SimpleITK as sitk
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

from bioscart_minimal.dicom_io import CTImage
from bioscart_minimal.dose_evaluation import evaluate_regions_on_rtdose
from bioscart_minimal.regions import RegionMask


def _write_rtdose(path, dose_counts):
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.2"
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds = FileDataset(str(path), {}, file_meta=file_meta, preamble=b"\0" * 128)
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    ds.SOPClassUID = file_meta.MediaStorageSOPClassUID
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.Modality = "RTDOSE"
    ds.FrameOfReferenceUID = "1.2.3"
    ds.Rows = dose_counts.shape[1]
    ds.Columns = dose_counts.shape[2]
    ds.NumberOfFrames = dose_counts.shape[0]
    ds.PixelSpacing = [1.0, 1.0]
    ds.ImagePositionPatient = [0.0, 0.0, 0.0]
    ds.GridFrameOffsetVector = [float(i) for i in range(dose_counts.shape[0])]
    ds.DoseUnits = "GY"
    ds.DoseType = "PHYSICAL"
    ds.DoseSummationType = "PLAN"
    ds.DoseGridScaling = 0.1
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0
    ds.PixelData = dose_counts.astype(np.uint16).tobytes()
    ds.save_as(str(path), write_like_original=False)


class TestDoseEvaluation:

    def test_evaluates_region_stats(self, tmp_path):
        dose_counts = np.arange(2 * 5 * 5, dtype=np.uint16).reshape(2, 5, 5)
        rtdose = tmp_path / "dose.dcm"
        _write_rtdose(rtdose, dose_counts)

        image = sitk.GetImageFromArray(np.zeros((2, 5, 5), dtype=np.float32))
        image.SetSpacing((1.0, 1.0, 1.0))
        image.SetOrigin((0.0, 0.0, 0.0))
        ct = CTImage(
            array_hu=np.zeros((2, 5, 5), dtype=np.float32),
            sitk_image=image,
            slice_paths=[],
            series_instance_uid="ct",
            origin_xyz_mm=(0.0, 0.0, 0.0),
            spacing_zyx_mm=(1.0, 1.0, 1.0),
            orientation=(1, 0, 0, 0, 1, 0),
            frame_of_reference_uid="1.2.3",
        )
        gtv = np.zeros((2, 5, 5), dtype=bool)
        gtv[:, 1:4, 1:4] = True
        region_mask = np.zeros_like(gtv)
        region_mask[:, 2:4, 2:4] = True
        regions = {
            "BSCART_Test": RegionMask(
                name="BSCART_Test",
                mask=region_mask,
                role="test",
                evidence_level="test",
                volume_cc=1.0,
                notes=[],
            )
        }
        report = evaluate_regions_on_rtdose(rtdose, ct, gtv, regions, tmp_path)
        assert report["frame_of_reference_match"] is True
        stats = {item["region"]: item for item in report["region_dose_stats"]}
        assert stats["GTV"]["voxel_count"] == 18
        assert stats["BSCART_Test"]["mean_gy"] is not None

