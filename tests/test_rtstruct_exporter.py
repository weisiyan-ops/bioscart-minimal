"""Tests for BioSCART RTSTRUCT export."""

from __future__ import annotations

import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence

from bioscart_minimal.regions import build_five_gtv_regions
from bioscart_minimal.rtstruct_exporter import export_bioscart_rtstruct


def _reference_rtstruct() -> Dataset:
    ds = Dataset()
    ds.Modality = "RTSTRUCT"
    ds.PatientID = "TEST"
    ds.StudyInstanceUID = "1.2.840.113619.2.55.3"
    frame = Dataset()
    frame.FrameOfReferenceUID = "1.2.3.4"
    ds.ReferencedFrameOfReferenceSequence = Sequence([frame])
    return ds


class TestRTSTRUCTExporter:

    def test_exports_importable_rtstruct(self, tmp_path, synthetic_ct_image, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        synthetic_ct_image.frame_of_reference_uid = "1.2.3.4"
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        out = export_bioscart_rtstruct(
            regions=regions,
            ct=synthetic_ct_image,
            reference_rtstruct=_reference_rtstruct(),
            output_dir=tmp_path,
        )
        ds = pydicom.dcmread(out)
        assert ds.Modality == "RTSTRUCT"
        names = {str(item.ROIName) for item in ds.StructureSetROISequence}
        assert "BSCART_Core" in names
        assert "BSCART_Rim_5mm" in names
        assert ds.StructureSetDescription.startswith("Research only")

