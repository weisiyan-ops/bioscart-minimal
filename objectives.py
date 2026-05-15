"""Vendor-neutral BioSCART research objective templates."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .dicom_io import GTVMask
from .regions import RegionMask


@dataclass
class ObjectiveSuggestion:
    """One research-only objective suggestion for TPS review."""

    structure: str
    suggested_role: str
    objective_class: str
    dose_value: float | None = None
    dose_units: str = "Gy"
    requires_physician_approval: bool = True
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "structure": self.structure,
            "suggested_role": self.suggested_role,
            "objective_class": self.objective_class,
            "dose_value": self.dose_value,
            "dose_units": self.dose_units,
            "requires_physician_approval": self.requires_physician_approval,
            "notes": self.notes,
        }


def build_objective_template(
    gtv: GTVMask,
    regions: dict[str, RegionMask],
    protocol_id: str = "BIO-SCART-MINIMAL-v0",
) -> dict[str, Any]:
    """Create a vendor-neutral research objective template.

    The template intentionally avoids prescription doses unless a protocol later
    supplies them.
    """
    suggestions = [_gtv_suggestion(gtv)]
    for region in regions.values():
        suggestions.append(_region_suggestion(region))

    return {
        "schema_version": "bioscart.objectives.v0.2",
        "clinical_status": "research_only_not_for_treatment_decisions",
        "protocol_id": protocol_id,
        "source_gtv": {
            "roi_name": gtv.roi_name,
            "roi_number": gtv.roi_number,
            "volume_cc": round(gtv.volume_cc, 3),
        },
        "objective_policy": {
            "prescription_doses_included": False,
            "final_tps_planning_required": True,
            "physician_physicist_review_required": True,
        },
        "objectives": [item.to_dict() for item in suggestions],
    }


def write_objective_template(template: dict[str, Any], output_dir: Path) -> Path:
    """Write objective template JSON."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "bioscart_objectives.json"
    path.write_text(json.dumps(template, indent=2), encoding="utf-8")
    return path


def _gtv_suggestion(gtv: GTVMask) -> ObjectiveSuggestion:
    return ObjectiveSuggestion(
        structure="GTV",
        suggested_role="source_tumor_volume",
        objective_class="maintain_protocol_gtv_coverage",
        notes=[
            f"Source clinical GTV ROI: {gtv.roi_name}",
            "Use the clinical protocol for target coverage; BioSCART does not define prescription.",
        ],
    )


def _region_suggestion(region: RegionMask) -> ObjectiveSuggestion:
    role_to_class = {
        "border_or_valley_candidate": "review_for_low_intermediate_valley_or_border_dose",
        "peak_candidate": "review_for_high_dose_vertex_candidate",
        "necrosis_or_hypoxia_surrogate_review": "review_as_uncertain_surrogate_not_prescription",
        "viable_or_enhancing_surrogate_review": "review_for_viable_surrogate_overlap",
        "heterogeneity_candidate": "review_for_heterogeneity_candidate",
    }
    return ObjectiveSuggestion(
        structure=region.name,
        suggested_role=region.role,
        objective_class=role_to_class.get(region.role, "manual_review"),
        notes=region.notes,
    )

