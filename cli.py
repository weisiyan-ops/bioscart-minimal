"""Command line entrypoint for the minimal BioSCART prototype."""

from __future__ import annotations

import argparse
from pathlib import Path

from .service import BioSCARTRunConfig, BioSCARTService


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="BioSCART minimal CT-only research prototype"
    )
    parser.add_argument("--dicom-dir", required=True, help="Input DICOM directory")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--gtv-name", default=None, help="Exact GTV ROI name to use")
    parser.add_argument("--ct-series-uid", default=None, help="CT SeriesInstanceUID to use")
    parser.add_argument("--rtstruct", default=None, help="Specific RTSTRUCT path")
    parser.add_argument("--rim-mm", type=float, default=5.0, help="Inner rim shell thickness in mm")
    parser.add_argument("--pyradiomics-params", default=None, help="Optional PyRadiomics YAML parameter file")
    parser.add_argument("--rtdose", default=None, help="Optional RTDOSE path for TPS round-trip evaluation")
    parser.add_argument("--ice3-metrics", default=None, help="Optional ICE3/EDIC metrics JSON for combined review")
    parser.add_argument("--strict-tps-geometry", action="store_true", help="Treat missing/mismatched DICOM geometry as hard failure")
    parser.add_argument("--no-rtstruct-export", action="store_true", help="Disable BioSCART RTSTRUCT export")
    args = parser.parse_args(argv)

    config = BioSCARTRunConfig(
        dicom_dir=Path(args.dicom_dir),
        out_dir=Path(args.out_dir),
        gtv_name=args.gtv_name,
        ct_series_uid=args.ct_series_uid,
        rtstruct_path=Path(args.rtstruct) if args.rtstruct else None,
        rtdose_path=Path(args.rtdose) if args.rtdose else None,
        ice3_metrics_path=Path(args.ice3_metrics) if args.ice3_metrics else None,
        rim_mm=args.rim_mm,
        pyradiomics_params=Path(args.pyradiomics_params) if args.pyradiomics_params else None,
        strict_tps_geometry=args.strict_tps_geometry,
        export_rtstruct=not args.no_rtstruct_export,
    )
    result = BioSCARTService().run_case(config)

    print(f"BioSCART minimal complete: {result.manifest_path}")
    print(f"Recommendations: {result.recommendations_path}")
    print(f"Objectives: {result.objective_template_path}")
    if result.rtstruct_export_path:
        print(f"BioSCART RTSTRUCT: {result.rtstruct_export_path}")
    if result.dose_evaluation_path:
        print(f"Dose evaluation: {result.dose_evaluation_path}")
    if result.combined_review_path:
        print(f"Combined review: {result.combined_review_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
