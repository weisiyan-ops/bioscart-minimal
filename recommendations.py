"""Research-only recommendation generator for BioSCART minimal."""

from __future__ import annotations

import json
from pathlib import Path

from .dicom_io import CTImage, GTVMask
from .regions import MIN_GTV_VOLUME_CC_FOR_RADIOMICS, RegionMask


def build_recommendations(
    ct: CTImage,
    gtv: GTVMask,
    regions: dict[str, RegionMask],
    radiomics_available: bool,
) -> dict:
    """Create transparent, research-only BioSCART recommendations."""
    recs = {
        "clinical_status": "research_only_not_for_treatment_decisions",
        "case_summary": {
            "ct_series_instance_uid": ct.series_instance_uid,
            "gtv_roi_name": gtv.roi_name,
            "gtv_volume_cc": round(gtv.volume_cc, 3),
            "radiomics_available": radiomics_available,
            "gtv_volume_adequate_for_radiomics": gtv.volume_cc >= MIN_GTV_VOLUME_CC_FOR_RADIOMICS,
        },
        "global_recommendations": [
            "Use CT-only BioSCART regions as review candidates, not validated biologic truth.",
            "Do not label hot/cold/excluded immune phenotype without pathology/IHC or spatial-omics validation.",
            "Do not label confirmed hypoxia without hypoxia PET, validated MRI surrogate, or tissue validation.",
            "Keep final planning and prescription decisions inside the clinical TPS and protocol review workflow.",
        ],
        "region_recommendations": [],
        "warnings": ct.warnings + gtv.warnings,
    }

    for name, region in regions.items():
        action = _action_for_region(name)
        recs["region_recommendations"].append({
            "region": name,
            "role": region.role,
            "volume_cc": round(region.volume_cc, 3),
            "suggested_research_action": action,
            "notes": region.notes,
        })

    return recs


def write_recommendations(recommendations: dict, output_dir: Path) -> tuple[Path, Path]:
    """Write JSON and Markdown recommendations."""
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "bioscart_recommendations.json"
    md_path = output_dir / "bioscart_recommendations.md"

    json_path.write_text(json.dumps(recommendations, indent=2), encoding="utf-8")

    lines = [
        "# BioSCART Minimal Recommendations",
        "",
        "**Status:** Research use only. Not for clinical treatment decision-making.",
        "",
        "## Case Summary",
        "",
        f"- GTV ROI: `{recommendations['case_summary']['gtv_roi_name']}`",
        f"- GTV volume: `{recommendations['case_summary']['gtv_volume_cc']} cc`",
        f"- PyRadiomics available: `{recommendations['case_summary']['radiomics_available']}`",
        "",
        "## Global Recommendations",
        "",
    ]
    lines.extend(f"- {item}" for item in recommendations["global_recommendations"])
    lines.extend(["", "## Region Recommendations", ""])

    for item in recommendations["region_recommendations"]:
        lines.extend([
            f"### {item['region']}",
            "",
            f"- Role: `{item['role']}`",
            f"- Volume: `{item['volume_cc']} cc`",
            f"- Suggested research action: {item['suggested_research_action']}",
            "",
            "Notes:",
        ])
        lines.extend(f"- {note}" for note in item["notes"])
        lines.append("")

    if recommendations["warnings"]:
        lines.extend(["## Warnings", ""])
        lines.extend(f"- {warning}" for warning in recommendations["warnings"])
        lines.append("")

    md_path.write_text("\n".join(lines), encoding="utf-8")
    return json_path, md_path


def _action_for_region(name: str) -> str:
    if name == "BSCART_Rim_5mm":
        return "Review as border/valley candidate; preserve as immune-access zone when protocol and OAR constraints permit."
    if name == "BSCART_Core":
        return "Review as geometric peak candidate; use only after OAR and dose-gradient feasibility checks."
    if name == "BSCART_CT_Low_Q25":
        return "Classify as uncertain necrosis/hypoxia surrogate; prioritize for validation, not automatic escalation."
    if name == "BSCART_CT_High_Q75":
        return "Review as viable/enhancing surrogate; consider overlap with core/texture regions for peak candidate ranking."
    if name == "BSCART_Texture_High_Q75":
        return "Review as heterogeneous/resistant candidate; validate against pathology before biologic dose painting claims."
    return "Review manually."

