"""ICE3-ready combined plan review scaffolding.

This module does not run ICE3 itself. It defines the bridge contract for
combining BioSCART region dose metrics with externally produced EDIC/ICE3
metrics.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_ice3_metrics(path: Path | str | None) -> dict[str, Any] | None:
    """Load external ICE3/EDIC metrics JSON if provided."""
    if path is None:
        return None
    return json.loads(Path(path).read_text(encoding="utf-8"))


def build_combined_plan_review(
    dose_evaluation: dict[str, Any] | None,
    ice3_metrics: dict[str, Any] | None = None,
    comparator: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a combined BioSCART + ICE3 research plan review.

    The output supports null-result reporting. It does not compute clinical
    risk or make treatment recommendations.
    """
    review = {
        "schema_version": "bioscart.combined_review.v0.1",
        "clinical_status": "research_only_not_for_treatment_decisions",
        "axes": {
            "bioscart_region_dose": dose_evaluation is not None,
            "ice3_edic": ice3_metrics is not None,
            "comparator_available": comparator is not None,
        },
        "interpretation_policy": [
            "No clinical risk score is produced.",
            "No plan is selected automatically.",
            "Null or no-benefit results are valid research outputs.",
        ],
        "bioscart_summary": _summarize_dose_eval(dose_evaluation),
        "ice3_summary": _summarize_ice3(ice3_metrics),
        "null_result_status": _null_result_status(dose_evaluation, ice3_metrics, comparator),
    }
    return review


def write_combined_plan_review(review: dict[str, Any], output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "bioscart_ice3_combined_review.json"
    path.write_text(json.dumps(review, indent=2), encoding="utf-8")
    return path


def _summarize_dose_eval(dose_evaluation: dict[str, Any] | None) -> dict[str, Any]:
    if dose_evaluation is None:
        return {"available": False}
    stats = dose_evaluation.get("region_dose_stats", [])
    return {
        "available": True,
        "region_count": len(stats),
        "regions": [
            {
                "region": item.get("region"),
                "mean_gy": item.get("mean_gy"),
                "d95_gy": item.get("d95_gy"),
                "d5_gy": item.get("d5_gy"),
            }
            for item in stats
        ],
    }


def _summarize_ice3(ice3_metrics: dict[str, Any] | None) -> dict[str, Any]:
    if ice3_metrics is None:
        return {"available": False}
    keys_of_interest = [
        "edic_gy",
        "edric_gy",
        "mean_blood_dose_gy",
        "lymphocyte_survival_fraction",
        "alc_predicted_change",
    ]
    summary = {key: ice3_metrics.get(key) for key in keys_of_interest if key in ice3_metrics}
    return {
        "available": True,
        "provided_metric_count": len(ice3_metrics),
        "selected_metrics": summary,
    }


def _null_result_status(
    dose_evaluation: dict[str, Any] | None,
    ice3_metrics: dict[str, Any] | None,
    comparator: dict[str, Any] | None,
) -> dict[str, Any]:
    if comparator is None:
        return {
            "can_assess_benefit": False,
            "reason": "No comparator plan metrics were provided.",
        }
    return {
        "can_assess_benefit": True,
        "allowed_outcomes": [
            "candidate_better",
            "candidate_worse",
            "no_significant_difference",
            "inconclusive",
        ],
        "note": "Statistical comparison is protocol-dependent and not performed by this bridge.",
    }

