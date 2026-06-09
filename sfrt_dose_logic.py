"""SFRT/SCART SIB dose prescription logic for BioSCART regions.

Maps CT-based tumor subregions to spatially varying SCART dose levels
(peak / intermediate / valley) and generates Eclipse-compatible SIB
optimization objectives. Integrates with ICE3 for immune safety checking.

Research use only. Not for clinical treatment decisions without physician
and physicist review.

References:
    - Wu et al., "GRID/Lattice/SCART dosimetric comparison" (Yan co-author)
    - Green Journal SFRT review: 187 studies, 3842 patients, PVDR >= 3-5:1
    - Eclipse VMAT SIB workflow (Varian RapidArc Operations pp.58-62)
    - SCART spatial transcriptomics: periphery 67.2% immune infiltration
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .regions import RegionMask


# ─── SCART Dose Level Presets ───────────────────────────────────────

@dataclass
class SCARTDoseLevel:
    """One dose tier in the SFRT/SCART SIB prescription."""
    name: str
    dose_per_fraction_gy: float
    total_dose_gy: float
    fractions: int
    role: str
    rationale: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "dose_per_fraction_gy": self.dose_per_fraction_gy,
            "total_dose_gy": self.total_dose_gy,
            "fractions": self.fractions,
            "role": self.role,
            "rationale": self.rationale,
        }


@dataclass
class SCARTProtocol:
    """A complete SCART/SFRT dose painting protocol."""
    protocol_id: str
    site: str
    description: str
    fractions: int
    levels: dict[str, SCARTDoseLevel]
    pvdr_target: float
    ice3_alc_nadir_threshold: float | None = None
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "protocol_id": self.protocol_id,
            "site": self.site,
            "description": self.description,
            "fractions": self.fractions,
            "levels": {k: v.to_dict() for k, v in self.levels.items()},
            "pvdr_target": self.pvdr_target,
            "ice3_alc_nadir_threshold": self.ice3_alc_nadir_threshold,
            "notes": self.notes,
        }


# ─── Built-in Protocols ─────────────────────────────────────────────

def scart_lung_protocol(fractions: int = 5) -> SCARTProtocol:
    """SCART protocol for bulky lung tumors (single-fraction or hypofractionated)."""
    return SCARTProtocol(
        protocol_id="SCART-LUNG-v1",
        site="lung",
        description="SCART dose painting for bulky NSCLC: ablative core + immune-preserving periphery",
        fractions=fractions,
        levels={
            "peak": SCARTDoseLevel(
                name="peak",
                dose_per_fraction_gy=15.0 / fractions * fractions if fractions == 1 else 66.7 / fractions,
                total_dose_gy=15.0 if fractions == 1 else 66.7,
                fractions=fractions,
                role="ablative_core",
                rationale="Ablative dose to radioresistant hypoxic core; overcome hypoxia-mediated resistance",
            ),
            "intermediate": SCARTDoseLevel(
                name="intermediate",
                dose_per_fraction_gy=8.0 / fractions * fractions if fractions == 1 else 40.0 / fractions,
                total_dose_gy=8.0 if fractions == 1 else 40.0,
                fractions=fractions,
                role="transition_zone",
                rationale="Moderate dose to heterogeneous/viable tumor; balance cell kill and vascular preservation",
            ),
            "valley": SCARTDoseLevel(
                name="valley",
                dose_per_fraction_gy=3.0 / fractions * fractions if fractions == 1 else 20.0 / fractions,
                total_dose_gy=3.0 if fractions == 1 else 20.0,
                fractions=fractions,
                role="immune_sanctuary",
                rationale="Low dose preserves tumor vasculature and immune cell access; enables bystander/abscopal effect",
            ),
        },
        pvdr_target=5.0,
        ice3_alc_nadir_threshold=0.5,
        notes=[
            "Based on LATTICE SIB: 66.7 Gy to vertices / 20 Gy to PTV over 5 fx",
            "SCART replaces geometric vertices with CT-guided biologic regions",
            "Spatial transcriptomics: SCART periphery = 67.2% immune infiltration vs SBRT 44.1%",
            "Adjust doses per institutional protocol and OAR constraints",
        ],
    )


def scart_esophagus_protocol(fractions: int = 28) -> SCARTProtocol:
    """SCART protocol for esophageal cancer (conventional fractionation SIB).

    CRITICAL LUMEN SAFETY: Peak dose is applied ONLY to BSCART_Core
    (geometrically interior, away from luminal surface). The Rim region
    (which includes the luminal mucosal surface) receives valley/standard
    dose to prevent mucosal breakdown, fistula, and perforation.

    ARTDECO (JCO 2021) boosted the ENTIRE GTV to 61.6 Gy and FAILED due to
    luminal toxicity. BioSCART-E escalates only the Core (away from lumen)
    while de-escalating the Rim (luminal surface) — biologically guided
    spatial dose painting instead of uniform escalation.
    """
    return SCARTProtocol(
        protocol_id="SCART-ESOPH-v2",
        site="esophagus",
        description="BioSCART-E: CT-guided SIB for esophageal CRT — lumen-sparing core escalation + ICE3 immune protection",
        fractions=fractions,
        levels={
            "peak": SCARTDoseLevel(
                name="peak",
                dose_per_fraction_gy=2.1,
                total_dose_gy=round(2.1 * fractions, 1),
                fractions=fractions,
                role="ablative_core",
                rationale="Modest escalation (58.8 Gy) to geometric Core ONLY — away from luminal surface. "
                          "Matches German SIB study (safe at 58.8 Gy). "
                          "ARTDECO failed at 61.6 Gy because it escalated EVERYWHERE including lumen.",
            ),
            "intermediate": SCARTDoseLevel(
                name="intermediate",
                dose_per_fraction_gy=1.8,
                total_dose_gy=round(1.8 * fractions, 1),
                fractions=fractions,
                role="standard_coverage",
                rationale="Standard CRT dose (50.4 Gy) to viable tumor regions and heterogeneous zones",
            ),
            "valley": SCARTDoseLevel(
                name="valley",
                dose_per_fraction_gy=1.6,
                total_dose_gy=round(1.6 * fractions, 1),
                fractions=fractions,
                role="lumen_sparing_immune_sanctuary",
                rationale="De-escalated dose to Rim (luminal surface) and low-density regions. "
                          "Protects esophageal mucosa from fistula/perforation while preserving "
                          "immune cell access. Also reduces blood dose (esophageal CRT causes "
                          "severe lymphopenia due to large mediastinal fields).",
            ),
        },
        pvdr_target=1.3,
        ice3_alc_nadir_threshold=0.5,
        notes=[
            "LUMEN SAFETY: Peak dose ONLY to BSCART_Core (interior, away from lumen)",
            "Rim region includes luminal mucosal surface -- receives valley dose (lumen protection)",
            "ARTDECO (JCO 2021): 61.6 Gy unselected boost → fistula, bleeding, no benefit",
            "German SIB study: 58.8 Gy SIB was safe and feasible (peak matches 58.8 Gy here)",
            "BioSCART-E innovation: escalate WHERE safe (Core) + protect WHERE vulnerable (Rim/lumen)",
            "ICE3 checks that escalation does not cause unacceptable lymphopenia",
            "Peak 58.8 Gy is below the 61.6 Gy ARTDECO failure threshold",
            "Valley 44.8 Gy provides lumen protection while maintaining tumor coverage",
        ],
    )


def scart_single_fraction_protocol() -> SCARTProtocol:
    """Classic single-fraction SCART (Lattice-style) for bulky tumors."""
    return SCARTProtocol(
        protocol_id="SCART-1FX-v1",
        site="any_bulky",
        description="Single-fraction SCART: 15 Gy peak / 3 Gy valley for bulky tumors",
        fractions=1,
        levels={
            "peak": SCARTDoseLevel(
                name="peak", dose_per_fraction_gy=15.0, total_dose_gy=15.0,
                fractions=1, role="ablative_core",
                rationale="Ablative dose to vertices/core; standard LATTICE peak dose",
            ),
            "intermediate": SCARTDoseLevel(
                name="intermediate", dose_per_fraction_gy=7.5, total_dose_gy=7.5,
                fractions=1, role="transition_zone",
                rationale="50% of peak dose to high-texture/viable regions",
            ),
            "valley": SCARTDoseLevel(
                name="valley", dose_per_fraction_gy=3.0, total_dose_gy=3.0,
                fractions=1, role="immune_sanctuary",
                rationale="Valley dose preserves vasculature and immune infiltration",
            ),
        },
        pvdr_target=5.0,
        ice3_alc_nadir_threshold=0.8,
        notes=[
            "Eclipse implementation: 15 Gy to lattice spheres / 2-3 Gy valley in single fraction",
            "SCART replaces geometric spheres with CT-guided biologic regions",
        ],
    )


PROTOCOLS = {
    "lung": scart_lung_protocol,
    "esophagus": scart_esophagus_protocol,
    "single_fraction": scart_single_fraction_protocol,
}


# ─── Region-to-Dose Mapping ────────────────────────────────────────

@dataclass
class RegionDoseAssignment:
    """Maps a BioSCART region to a SCART dose level."""
    region_name: str
    dose_level: str
    dose_per_fraction_gy: float
    total_dose_gy: float
    confidence: str
    rationale: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "region_name": self.region_name,
            "dose_level": self.dose_level,
            "dose_per_fraction_gy": self.dose_per_fraction_gy,
            "total_dose_gy": self.total_dose_gy,
            "confidence": self.confidence,
            "rationale": self.rationale,
        }


def assign_doses_to_regions(
    regions: dict[str, RegionMask],
    protocol: SCARTProtocol,
) -> list[RegionDoseAssignment]:
    """Map each BioSCART region to a SCART dose level based on its role."""

    role_to_level = {
        "peak_candidate": "peak",
        "viable_or_enhancing_surrogate_review": "peak",
        "heterogeneity_candidate": "intermediate",
        "border_or_valley_candidate": "valley",
        "necrosis_or_hypoxia_surrogate_review": "valley",
    }

    assignments = []
    for name, region in regions.items():
        level_key = role_to_level.get(region.role, "intermediate")
        level = protocol.levels[level_key]

        confidence = "geometric" if region.evidence_level == "geometric_ct_only" else "ct_surrogate"

        assignments.append(RegionDoseAssignment(
            region_name=name,
            dose_level=level_key,
            dose_per_fraction_gy=level.dose_per_fraction_gy,
            total_dose_gy=level.total_dose_gy,
            confidence=confidence,
            rationale=f"{region.role} -> {level_key}: {level.rationale}",
        ))

    return assignments


# ─── Eclipse SIB Objective Generator ───────────────────────────────

@dataclass
class SIBObjective:
    """One Eclipse-compatible SIB optimization objective."""
    structure: str
    objective_type: str
    volume_pct: float | None
    dose_cgy: float
    priority: int
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "structure": self.structure,
            "objective_type": self.objective_type,
            "volume_pct": self.volume_pct,
            "dose_cgy": self.dose_cgy,
            "priority": self.priority,
            "notes": self.notes,
        }


def generate_sib_objectives(
    assignments: list[RegionDoseAssignment],
    protocol: SCARTProtocol,
    oar_constraints: dict[str, list[dict]] | None = None,
) -> list[SIBObjective]:
    """Generate Eclipse-compatible SIB optimization objectives.

    Produces Lower (minimum) and Upper (maximum) objectives for each
    region, plus standard OAR constraints. These can be loaded into
    Eclipse's Optimization dialog.
    """
    objectives = []

    for a in assignments:
        total_cgy = a.total_dose_gy * 100

        # Lower objective: minimum dose to this region
        lower_dose = total_cgy * 0.95
        objectives.append(SIBObjective(
            structure=a.region_name,
            objective_type="Lower",
            volume_pct=100.0,
            dose_cgy=round(lower_dose),
            priority=100 if a.dose_level == "peak" else 90 if a.dose_level == "intermediate" else 80,
            notes=f"SIB {a.dose_level}: minimum coverage",
        ))

        # Upper objective: maximum dose
        upper_dose = total_cgy * 1.07
        objectives.append(SIBObjective(
            structure=a.region_name,
            objective_type="Upper",
            volume_pct=0.0,
            dose_cgy=round(upper_dose),
            priority=100 if a.dose_level == "peak" else 90 if a.dose_level == "intermediate" else 80,
            notes=f"SIB {a.dose_level}: hot spot limit",
        ))

    # Standard OAR constraints
    default_oars = {
        "Cord": [{"type": "Upper", "vol": 0.0, "dose_cgy": 4500, "priority": 85}],
        "Total_Lung": [
            {"type": "Upper", "vol": 20.0, "dose_cgy": 2000, "priority": 70},
            {"type": "Mean", "vol": None, "dose_cgy": 2000, "priority": 70},
        ],
        "Heart": [
            {"type": "Upper", "vol": 33.0, "dose_cgy": 3000, "priority": 70},
            {"type": "Mean", "vol": None, "dose_cgy": 2000, "priority": 70},
        ],
    }

    oar_specs = oar_constraints if oar_constraints else default_oars
    for structure, constraints in oar_specs.items():
        for c in constraints:
            objectives.append(SIBObjective(
                structure=structure,
                objective_type=c["type"],
                volume_pct=c.get("vol"),
                dose_cgy=c["dose_cgy"],
                priority=c.get("priority", 70),
                notes="Standard OAR constraint",
            ))

    # NTO
    objectives.append(SIBObjective(
        structure="NTO",
        objective_type="Automatic",
        volume_pct=None,
        dose_cgy=0,
        priority=100,
        notes="Normal Tissue Objective: Automatic mode, priority matches peak target",
    ))

    return objectives


# ─── ICE3 Immune Safety Check ──────────────────────────────────────

@dataclass
class ImmuneSafetyResult:
    """Result of checking SFRT dose plan against ICE3 immune metrics."""
    safe: bool
    ice3_available: bool
    alc_nadir_predicted: float | None
    blood_dose_gy: float | None
    threshold: float | None
    recommendation: str
    details: list[str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "safe": self.safe,
            "ice3_available": self.ice3_available,
            "alc_nadir_predicted": self.alc_nadir_predicted,
            "blood_dose_gy": self.blood_dose_gy,
            "threshold": self.threshold,
            "recommendation": self.recommendation,
            "details": self.details,
        }


def check_immune_safety(
    ice3_metrics: dict[str, Any] | None,
    protocol: SCARTProtocol,
) -> ImmuneSafetyResult:
    """Check whether the SFRT plan is immune-safe per ICE3 metrics."""

    if ice3_metrics is None:
        return ImmuneSafetyResult(
            safe=True,
            ice3_available=False,
            alc_nadir_predicted=None,
            blood_dose_gy=None,
            threshold=protocol.ice3_alc_nadir_threshold,
            recommendation="ICE3 metrics not provided; immune safety not assessed. "
                           "Run ICE3 on the RTDOSE to evaluate blood dose.",
            details=["No ICE3 metrics JSON provided via --ice3-metrics"],
        )

    alc_nadir = ice3_metrics.get("alc_predicted_nadir") or ice3_metrics.get("alc_nadir")
    blood_dose = ice3_metrics.get("mean_blood_dose_gy") or ice3_metrics.get("edic_gy")
    threshold = protocol.ice3_alc_nadir_threshold

    details = []
    if alc_nadir is not None:
        details.append(f"ICE3 predicted ALC nadir: {alc_nadir:.2f} x10^9/L")
    if blood_dose is not None:
        details.append(f"ICE3 mean blood dose: {blood_dose:.3f} Gy/fx")
    if threshold is not None:
        details.append(f"Protocol ALC nadir threshold: {threshold:.2f} x10^9/L")

    if alc_nadir is not None and threshold is not None:
        safe = alc_nadir >= threshold
        if safe:
            recommendation = (
                f"IMMUNE SAFE: Predicted ALC nadir ({alc_nadir:.2f}) >= threshold ({threshold:.2f}). "
                f"SFRT dose escalation acceptable from immune perspective."
            )
        else:
            recommendation = (
                f"IMMUNE RISK: Predicted ALC nadir ({alc_nadir:.2f}) < threshold ({threshold:.2f}). "
                f"Consider: (1) reduce peak dose, (2) add avoidance sectors to protect major vessels, "
                f"(3) switch to fewer fractions, (4) accept elevated lymphopenia risk with documentation."
            )
    else:
        safe = True
        recommendation = (
            "ICE3 metrics provided but ALC nadir prediction not found. "
            "Cannot assess immune safety. Check ICE3 output format."
        )

    return ImmuneSafetyResult(
        safe=safe,
        ice3_available=True,
        alc_nadir_predicted=alc_nadir,
        blood_dose_gy=blood_dose,
        threshold=threshold,
        recommendation=recommendation,
        details=details,
    )


# ─── Full SFRT Planning Report ─────────────────────────────────────

def build_sfrt_plan(
    regions: dict[str, RegionMask],
    protocol: SCARTProtocol,
    ice3_metrics: dict[str, Any] | None = None,
    oar_constraints: dict[str, list[dict]] | None = None,
) -> dict[str, Any]:
    """Build complete SFRT/SCART dose plan from BioSCART regions.

    Returns a JSON-serializable dict with:
    - Protocol specification
    - Region-to-dose assignments
    - Eclipse SIB objectives
    - ICE3 immune safety check
    - PVDR (peak-to-valley dose ratio)
    """
    assignments = assign_doses_to_regions(regions, protocol)
    objectives = generate_sib_objectives(assignments, protocol, oar_constraints)
    immune_check = check_immune_safety(ice3_metrics, protocol)

    peak_doses = [a.total_dose_gy for a in assignments if a.dose_level == "peak"]
    valley_doses = [a.total_dose_gy for a in assignments if a.dose_level == "valley"]
    actual_pvdr = None
    if peak_doses and valley_doses:
        actual_pvdr = round(max(peak_doses) / min(valley_doses), 2)

    return {
        "schema_version": "bioscart.sfrt_plan.v0.1",
        "clinical_status": "research_only_not_for_treatment_decisions",
        "protocol": protocol.to_dict(),
        "region_dose_assignments": [a.to_dict() for a in assignments],
        "eclipse_sib_objectives": [o.to_dict() for o in objectives],
        "dose_summary": {
            "peak_total_gy": max(peak_doses) if peak_doses else None,
            "valley_total_gy": min(valley_doses) if valley_doses else None,
            "actual_pvdr": actual_pvdr,
            "target_pvdr": protocol.pvdr_target,
            "pvdr_meets_target": actual_pvdr >= protocol.pvdr_target if actual_pvdr else None,
        },
        "immune_safety": immune_check.to_dict(),
    }


def write_sfrt_plan(plan: dict[str, Any], output_dir: Path) -> Path:
    """Write SFRT plan to JSON."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "bioscart_sfrt_plan.json"
    path.write_text(json.dumps(plan, indent=2), encoding="utf-8")
    return path
