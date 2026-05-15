"""Tests for DICOM geometry validation."""

from __future__ import annotations

import pytest
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence

from bioscart_minimal.errors import DICOMGeometryError
from bioscart_minimal.geometry import validate_ct_rtstruct_geometry


def _rtstruct_with_frame(frame_uid: str) -> Dataset:
    ds = Dataset()
    item = Dataset()
    item.FrameOfReferenceUID = frame_uid
    ds.ReferencedFrameOfReferenceSequence = Sequence([item])
    return ds


class TestGeometryValidation:

    def test_matching_frame_passes(self, synthetic_ct_image):
        synthetic_ct_image.frame_of_reference_uid = "1.2.3"
        result = validate_ct_rtstruct_geometry(
            synthetic_ct_image,
            _rtstruct_with_frame("1.2.3"),
            strict=True,
        )
        assert result.passed

    def test_mismatched_frame_fails(self, synthetic_ct_image):
        synthetic_ct_image.frame_of_reference_uid = "1.2.3"
        result = validate_ct_rtstruct_geometry(
            synthetic_ct_image,
            _rtstruct_with_frame("9.9.9"),
            strict=False,
        )
        assert not result.passed
        assert result.errors

    def test_strict_missing_frame_raises(self, synthetic_ct_image):
        synthetic_ct_image.frame_of_reference_uid = None
        with pytest.raises(DICOMGeometryError):
            validate_ct_rtstruct_geometry(synthetic_ct_image, Dataset(), strict=True)

