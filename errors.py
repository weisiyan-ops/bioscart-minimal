"""Error taxonomy for BioSCART Minimal."""

from __future__ import annotations


class BioSCARTError(Exception):
    """Base exception for BioSCART errors."""


class BioSCARTInputError(BioSCARTError):
    """Input data are missing, ambiguous, or unsupported."""


class DICOMGeometryError(BioSCARTError):
    """DICOM geometry or frame-of-reference validation failed."""


class StructureExportError(BioSCARTError):
    """Derived structure export failed."""


class ObjectiveTemplateError(BioSCARTError):
    """Objective template creation failed."""


class DoseEvaluationError(BioSCARTError):
    """RTDOSE loading or region dose evaluation failed."""

