"""Batch runner for BioSCART Minimal v0 — runs multiple cases and produces a summary."""

from __future__ import annotations

import csv
import json
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bioscart_minimal.service import BioSCARTRunConfig, BioSCARTService


def find_rtdose(case_dir: Path) -> Path | None:
    """Find the first RTDOSE file in a case directory."""
    candidates = sorted(case_dir.glob("RD.*.dcm"))
    return candidates[0] if candidates else None


def run_batch(
    data_root: Path,
    output_root: Path,
    max_cases: int = 10,
    include_rtdose: bool = True,
) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    service = BioSCARTService()

    case_dirs = sorted(
        [d for d in data_root.iterdir() if d.is_dir()],
        key=lambda d: d.name,
    )[:max_cases]

    summary_rows: list[dict] = []
    summary_path = output_root / "batch_summary.csv"

    print(f"BioSCART Batch: {len(case_dirs)} cases from {data_root.name}")
    print(f"Output: {output_root}")
    print("-" * 70)

    for i, case_dir in enumerate(case_dirs, 1):
        case_id = case_dir.name
        case_out = output_root / case_id
        print(f"[{i}/{len(case_dirs)}] {case_id} ... ", end="", flush=True)

        t0 = time.time()
        row: dict = {"case_id": case_id, "status": "FAILED", "error": ""}

        try:
            rtdose_path = find_rtdose(case_dir) if include_rtdose else None

            config = BioSCARTRunConfig(
                dicom_dir=case_dir,
                out_dir=case_out,
                rtdose_path=rtdose_path,
            )
            result = service.run_case(config)

            manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
            region_csv = case_out / "bioscart_region_summary.csv"
            region_data = {}
            if region_csv.exists():
                with region_csv.open(encoding="utf-8") as f:
                    for r in csv.DictReader(f):
                        region_data[r["region"]] = float(r["volume_cc"])

            row.update({
                "status": "OK",
                "gtv_roi_name": manifest.get("gtv_roi_name", ""),
                "gtv_volume_cc": manifest.get("gtv_volume_cc", 0),
                "ct_slices": manifest.get("ct_slice_count", 0),
                "geometry_passed": manifest.get("geometry_validation", {}).get("passed", False),
                "rim_cc": region_data.get("BSCART_Rim_5mm", 0),
                "core_cc": region_data.get("BSCART_Core", 0),
                "ct_low_cc": region_data.get("BSCART_CT_Low_Q25", 0),
                "ct_high_cc": region_data.get("BSCART_CT_High_Q75", 0),
                "texture_cc": region_data.get("BSCART_Texture_High_Q75", 0),
                "rtdose_used": rtdose_path is not None,
                "dose_eval": result.dose_evaluation_path is not None,
                "rtstruct_exported": result.rtstruct_export_path is not None,
                "warnings": "; ".join(manifest.get("warnings", [])),
                "elapsed_sec": round(time.time() - t0, 1),
            })
            print(f"OK  GTV={row['gtv_volume_cc']:.1f}cc  {row['elapsed_sec']}s")

        except Exception as e:
            row["error"] = str(e)
            row["elapsed_sec"] = round(time.time() - t0, 1)
            print(f"FAILED  {e}")
            traceback.print_exc()

        summary_rows.append(row)

    all_fields = []
    for r in summary_rows:
        for k in r:
            if k not in all_fields:
                all_fields.append(k)

    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    ok = sum(1 for r in summary_rows if r["status"] == "OK")
    fail = len(summary_rows) - ok
    print("-" * 70)
    print(f"Batch complete: {ok} OK, {fail} FAILED")
    print(f"Summary: {summary_path}")
    return summary_path


if __name__ == "__main__":
    data_root = Path(
        r"C:\Users\wya245\OneDrive - University of Kentucky\Desktop\SaMD\Data"
        r"\Thoracic DICOM_Export Cross Ref V2.0"
    )
    output_root = Path(r"C:\Users\wya245\SaMD\bioscart_output")

    max_cases = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    run_batch(data_root, output_root, max_cases=max_cases, include_rtdose=True)
