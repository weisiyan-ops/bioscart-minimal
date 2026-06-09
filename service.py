"""Service API for BioSCART Minimal.

CLI, future GUI, and future TPS wrappers should use this service instead of
duplicating orchestration logic.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import SimpleITK as sitk

from .dicom_io import (
    find_gtv_roi,
    load_ct_series,
    load_rtstruct,
    rasterize_roi_to_ct_grid,
    scan_dicom_directory,
    select_largest_ct_series,
    sitk_mask_from_array,
)
from .dose_evaluation import evaluate_regions_on_rtdose
from .geometry import validate_ct_rtstruct_geometry
from .ice3_bridge import build_combined_plan_review, load_ice3_metrics, write_combined_plan_review
from .objectives import build_objective_template, write_objective_template
from .radiomics_runner import run_pyradiomics
from .recommendations import build_recommendations, write_recommendations
from .regions import build_five_gtv_regions, check_gtv_volume_for_radiomics
from .rtstruct_exporter import export_bioscart_rtstruct
from .rtdose_exporter import export_prescription_rtdose
from .sfrt_dose_logic import PROTOCOLS, assign_doses_to_regions, build_sfrt_plan, write_sfrt_plan


@dataclass
class BioSCARTRunConfig:
    """Configuration for a BioSCART Minimal case run."""

    dicom_dir: Path
    out_dir: Path
    gtv_name: str | None = None
    ct_series_uid: str | None = None
    rtstruct_path: Path | None = None
    rtdose_path: Path | None = None
    ice3_metrics_path: Path | None = None
    rim_mm: float = 5.0
    pyradiomics_params: Path | None = None
    strict_tps_geometry: bool = False
    export_rtstruct: bool = True
    sfrt_protocol: str | None = None


@dataclass
class BioSCARTRunResult:
    """Important output paths for a BioSCART Minimal run."""

    manifest_path: Path
    recommendations_path: Path
    objective_template_path: Path
    rtstruct_export_path: Path | None
    dose_evaluation_path: Path | None
    combined_review_path: Path | None
    sfrt_plan_path: Path | None = None
    sfrt_rtdose_path: Path | None = None


class BioSCARTService:
    """TPS-adjacent BioSCART Minimal service API."""

    def run_case(self, config: BioSCARTRunConfig) -> BioSCARTRunResult:
        out_dir = Path(config.out_dir)
        masks_dir = out_dir / "masks"
        radiomics_dir = out_dir / "radiomics"
        out_dir.mkdir(parents=True, exist_ok=True)
        masks_dir.mkdir(parents=True, exist_ok=True)

        inventory = scan_dicom_directory(config.dicom_dir)
        if config.ct_series_uid:
            ct_paths = inventory.ct_series.get(config.ct_series_uid)
            if not ct_paths:
                raise ValueError(f"Requested CT series not found: {config.ct_series_uid}")
            ct_series_uid = config.ct_series_uid
        else:
            ct_series_uid, ct_paths = select_largest_ct_series(inventory)

        rtstruct_path = config.rtstruct_path or _select_rtstruct(inventory.rtstructs)
        ct = load_ct_series(ct_paths)
        rtstruct = load_rtstruct(rtstruct_path)
        geometry = validate_ct_rtstruct_geometry(
            ct,
            rtstruct,
            strict=config.strict_tps_geometry,
        )

        roi_name, roi_number = find_gtv_roi(rtstruct, config.gtv_name)
        gtv = rasterize_roi_to_ct_grid(rtstruct, roi_number, roi_name, ct)

        volume_warnings = check_gtv_volume_for_radiomics(gtv.volume_cc)
        gtv.warnings.extend(volume_warnings)

        regions = build_five_gtv_regions(
            ct_hu=ct.array_hu,
            gtv_mask=gtv.mask,
            spacing_zyx_mm=ct.spacing_zyx_mm,
            rim_mm=config.rim_mm,
        )

        sitk.WriteImage(ct.sitk_image, str(out_dir / "ct_image.nii.gz"))
        sitk.WriteImage(sitk_mask_from_array(gtv.mask, ct.sitk_image), str(masks_dir / "GTV.nii.gz"))
        for name, region in regions.items():
            sitk.WriteImage(sitk_mask_from_array(region.mask, ct.sitk_image), str(masks_dir / f"{name}.nii.gz"))

        region_summary_path = out_dir / "bioscart_region_summary.csv"
        _write_region_summary(region_summary_path, gtv, regions)

        radiomics_result = run_pyradiomics(
            image=ct.sitk_image,
            gtv_mask=gtv.mask,
            regions=regions,
            output_dir=radiomics_dir,
            params_path=config.pyradiomics_params,
        )

        recs = build_recommendations(
            ct=ct,
            gtv=gtv,
            regions=regions,
            radiomics_available=radiomics_result.available,
        )
        rec_json, rec_md = write_recommendations(recs, out_dir)

        objective_template = build_objective_template(gtv, regions)
        objective_path = write_objective_template(objective_template, out_dir)

        rtstruct_export_path = None
        if config.export_rtstruct:
            rtstruct_export_path = export_bioscart_rtstruct(
                regions=regions,
                ct=ct,
                reference_rtstruct=rtstruct,
                output_dir=out_dir / "dicom_export",
            )

        dose_eval = None
        dose_eval_path = None
        if config.rtdose_path:
            dose_eval = evaluate_regions_on_rtdose(
                rtdose_path=config.rtdose_path,
                ct=ct,
                gtv_mask=gtv.mask,
                regions=regions,
                output_dir=out_dir,
            )
            dose_eval_path = Path(dose_eval["output_path"])

        combined_review_path = None
        ice3_metrics = None
        if dose_eval is not None or config.ice3_metrics_path:
            ice3_metrics = load_ice3_metrics(config.ice3_metrics_path)
            combined = build_combined_plan_review(dose_eval, ice3_metrics)
            combined_review_path = write_combined_plan_review(combined, out_dir)

        sfrt_plan_path = None
        sfrt_rtdose_path = None
        if config.sfrt_protocol:
            protocol_factory = PROTOCOLS.get(config.sfrt_protocol)
            if protocol_factory is None:
                raise ValueError(
                    f"Unknown SFRT protocol: {config.sfrt_protocol!r}. "
                    f"Available: {', '.join(PROTOCOLS.keys())}"
                )
            protocol = protocol_factory()
            assignments = assign_doses_to_regions(regions, protocol)
            sfrt_plan = build_sfrt_plan(
                regions=regions,
                protocol=protocol,
                ice3_metrics=ice3_metrics,
            )
            sfrt_plan_path = write_sfrt_plan(sfrt_plan, out_dir)
            sfrt_rtdose_path = export_prescription_rtdose(
                assignments=assignments,
                regions=regions,
                ct=ct,
                reference_rtstruct=rtstruct,
                output_dir=out_dir / "dicom_export",
            )

        manifest = {
            "schema_version": "bioscart.minimal.v0.2",
            "clinical_status": "research_only_not_for_treatment_decisions",
            "input_dicom_dir": str(config.dicom_dir),
            "ct_series_instance_uid": ct_series_uid,
            "ct_slice_count": len(ct_paths),
            "rtstruct_path": str(rtstruct_path),
            "rtdose_path": str(config.rtdose_path) if config.rtdose_path else None,
            "gtv_roi_name": roi_name,
            "gtv_roi_number": roi_number,
            "gtv_volume_cc": round(gtv.volume_cc, 3),
            "geometry_validation": geometry.to_dict(),
            "outputs": {
                "ct_image": str(out_dir / "ct_image.nii.gz"),
                "masks_dir": str(masks_dir),
                "region_summary_csv": str(region_summary_path),
                "radiomics_status_json": radiomics_result.status_json,
                "radiomics_features_csv": radiomics_result.features_csv,
                "recommendations_json": str(rec_json),
                "recommendations_md": str(rec_md),
                "objectives_json": str(objective_path),
                "bioscart_rtstruct": str(rtstruct_export_path) if rtstruct_export_path else None,
                "dose_evaluation_json": str(dose_eval_path) if dose_eval_path else None,
                "combined_review_json": str(combined_review_path) if combined_review_path else None,
                "sfrt_plan_json": str(sfrt_plan_path) if sfrt_plan_path else None,
                "sfrt_rtdose": str(sfrt_rtdose_path) if sfrt_rtdose_path else None,
            },
            "warnings": ct.warnings + gtv.warnings + geometry.warnings,
        }
        manifest_path = out_dir / "bioscart_minimal_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

        return BioSCARTRunResult(
            manifest_path=manifest_path,
            recommendations_path=rec_md,
            objective_template_path=objective_path,
            rtstruct_export_path=rtstruct_export_path,
            dose_evaluation_path=dose_eval_path,
            combined_review_path=combined_review_path,
            sfrt_plan_path=sfrt_plan_path,
            sfrt_rtdose_path=sfrt_rtdose_path,
        )


def _select_rtstruct(paths: list[Path]) -> Path:
    if not paths:
        raise ValueError("No RTSTRUCT found in DICOM directory")
    return paths[0]


def _write_region_summary(path: Path, gtv, regions) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["region", "role", "evidence_level", "volume_cc", "voxel_count", "notes"],
        )
        writer.writeheader()
        writer.writerow({
            "region": "GTV",
            "role": "source_tumor_mask",
            "evidence_level": "clinical_rtstruct",
            "volume_cc": round(gtv.volume_cc, 3),
            "voxel_count": int(gtv.mask.sum()),
            "notes": "Clinical source contour selected by GTV name matching.",
        })
        for region in regions.values():
            writer.writerow({
                "region": region.name,
                "role": region.role,
                "evidence_level": region.evidence_level,
                "volume_cc": round(region.volume_cc, 3),
                "voxel_count": int(region.mask.sum()),
                "notes": " | ".join(region.notes),
            })

