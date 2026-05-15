"""DICOM geometry validation for BioSCART Minimal."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from pydicom.dataset import Dataset

from .dicom_io import CTImage
from .errors import DICOMGeometryError


@dataclass
class GeometryValidationResult:
    """Structured result of CT/RTSTRUCT geometry validation."""

    passed: bool
    ct_frame_of_reference_uid: str | None = None
    rtstruct_frame_of_reference_uid: str | None = None
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "ct_frame_of_reference_uid": self.ct_frame_of_reference_uid,
            "rtstruct_frame_of_reference_uid": self.rtstruct_frame_of_reference_uid,
            "warnings": self.warnings,
            "errors": self.errors,
        }


def validate_ct_rtstruct_geometry(
    ct: CTImage,
    rtstruct: Dataset,
    strict: bool = False,
) -> GeometryValidationResult:
    """Validate CT and RTSTRUCT linkage before TPS-facing export.

    In non-strict mode this returns warnings for missing metadata so the
    retrospective research workflow can still run. In strict TPS mode, missing
    or mismatched FrameOfReferenceUID is a hard failure.
    """
    warnings: list[str] = []
    errors: list[str] = []

    ct_for_uid = ct.frame_of_reference_uid
    rs_for_uid = _rtstruct_frame_of_reference_uid(rtstruct)

    if not ct_for_uid:
        _add_issue(
            "CT FrameOfReferenceUID is missing; TPS-facing export requires explicit geometry linkage.",
            strict,
            warnings,
            errors,
        )
    if not rs_for_uid:
        _add_issue(
            "RTSTRUCT FrameOfReferenceUID is missing; TPS-facing export requires explicit geometry linkage.",
            strict,
            warnings,
            errors,
        )
    if ct_for_uid and rs_for_uid and ct_for_uid != rs_for_uid:
        errors.append(
            f"CT FrameOfReferenceUID ({ct_for_uid}) does not match RTSTRUCT FrameOfReferenceUID ({rs_for_uid})."
        )

    if not _is_identity_axial(ct.orientation):
        _add_issue(
            "Non-identity CT orientation is not supported for TPS-facing BioSCART export in minimal v0.",
            strict,
            warnings,
            errors,
        )

    if any(spacing <= 0 for spacing in ct.spacing_zyx_mm):
        errors.append(f"Invalid CT spacing: {ct.spacing_zyx_mm}")

    result = GeometryValidationResult(
        passed=not errors,
        ct_frame_of_reference_uid=ct_for_uid,
        rtstruct_frame_of_reference_uid=rs_for_uid,
        warnings=warnings,
        errors=errors,
    )
    if strict and errors:
        raise DICOMGeometryError("; ".join(errors))
    return result


def _rtstruct_frame_of_reference_uid(rtstruct: Dataset) -> str | None:
    if hasattr(rtstruct, "ReferencedFrameOfReferenceSequence") and rtstruct.ReferencedFrameOfReferenceSequence:
        uid = getattr(rtstruct.ReferencedFrameOfReferenceSequence[0], "FrameOfReferenceUID", None)
        if uid:
            return str(uid)
    if hasattr(rtstruct, "StructureSetROISequence"):
        for roi in rtstruct.StructureSetROISequence:
            uid = getattr(roi, "ReferencedFrameOfReferenceUID", None)
            if uid:
                return str(uid)
    return None


def _add_issue(
    message: str,
    strict: bool,
    warnings: list[str],
    errors: list[str],
) -> None:
    if strict:
        errors.append(message)
    else:
        warnings.append(message)


def _is_identity_axial(orientation: tuple[float, ...], tolerance: float = 1e-3) -> bool:
    expected = np.array([1, 0, 0, 0, 1, 0], dtype=float)
    observed = np.array(orientation, dtype=float)
    return observed.shape == expected.shape and np.allclose(observed, expected, atol=tolerance)

