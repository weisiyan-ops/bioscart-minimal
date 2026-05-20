"""Aggregate dose evaluation across BioSCART batch cases.

Answers: does the standard plan already differentiate dose by BioSCART region?
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np


def aggregate(output_root: Path) -> Path:
    case_dirs = sorted(
        [d for d in output_root.iterdir() if d.is_dir() and (d / "bioscart_dose_evaluation.json").exists()]
    )

    rows: list[dict] = []
    for case_dir in case_dirs:
        dose_path = case_dir / "bioscart_dose_evaluation.json"
        manifest_path = case_dir / "bioscart_minimal_manifest.json"
        dose_data = json.loads(dose_path.read_text(encoding="utf-8"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

        case_id = case_dir.name
        gtv_vol = manifest.get("gtv_volume_cc", 0)

        for region_stat in dose_data.get("region_dose_stats", []):
            rows.append({
                "case_id": case_id,
                "gtv_volume_cc": gtv_vol,
                "region": region_stat["region"],
                "voxel_count": region_stat["voxel_count"],
                "mean_gy": region_stat["mean_gy"],
                "min_gy": region_stat["min_gy"],
                "max_gy": region_stat["max_gy"],
                "d95_gy": region_stat["d95_gy"],
                "d50_gy": region_stat["d50_gy"],
                "d5_gy": region_stat["d5_gy"],
            })

    long_csv = output_root / "dose_aggregate_long.csv"
    fields = ["case_id", "gtv_volume_cc", "region", "voxel_count",
              "mean_gy", "min_gy", "max_gy", "d95_gy", "d50_gy", "d5_gy"]
    with long_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    regions_of_interest = ["GTV", "BSCART_Core", "BSCART_Rim_5mm",
                           "BSCART_CT_Low_Q25", "BSCART_CT_High_Q75",
                           "BSCART_Texture_High_Q75"]

    print("=" * 90)
    print("DOSE AGGREGATION: Core vs Rim vs CT-Intensity vs Texture")
    print("=" * 90)

    cases_with_core = [r for r in rows if r["region"] == "BSCART_Core" and r["voxel_count"] > 10]
    if not cases_with_core:
        print("No cases with meaningful Core region for comparison.")
        return long_csv

    eligible_cases = {r["case_id"] for r in cases_with_core}
    print(f"\nCases with meaningful Core region (>10 voxels): {len(eligible_cases)}/{len(case_dirs)}")

    print(f"\n{'Case':<12} {'GTV(cc)':>8} {'GTV Mean':>9} {'Core Mean':>10} {'Rim Mean':>9} {'Core-Rim':>9} {'CT-Low':>8} {'CT-High':>9} {'Texture':>8}")
    print("-" * 90)

    core_rim_diffs = []
    ctlow_cthigh_diffs = []

    for case_id in sorted(eligible_cases):
        case_rows = {r["region"]: r for r in rows if r["case_id"] == case_id}
        gtv = case_rows.get("GTV", {})
        core = case_rows.get("BSCART_Core", {})
        rim = case_rows.get("BSCART_Rim_5mm", {})
        ct_low = case_rows.get("BSCART_CT_Low_Q25", {})
        ct_high = case_rows.get("BSCART_CT_High_Q75", {})
        texture = case_rows.get("BSCART_Texture_High_Q75", {})

        core_mean = core.get("mean_gy", 0)
        rim_mean = rim.get("mean_gy", 0)
        diff = core_mean - rim_mean
        core_rim_diffs.append(diff)

        ctlow_mean = ct_low.get("mean_gy", 0)
        cthigh_mean = ct_high.get("mean_gy", 0)
        ctlow_cthigh_diffs.append(cthigh_mean - ctlow_mean)

        print(f"{case_id:<12} {gtv.get('gtv_volume_cc', 0):>8.1f} "
              f"{gtv.get('mean_gy', 0):>9.2f} "
              f"{core_mean:>10.2f} {rim_mean:>9.2f} "
              f"{diff:>+9.2f} "
              f"{ctlow_mean:>8.2f} {cthigh_mean:>9.2f} "
              f"{texture.get('mean_gy', 0):>8.2f}")

    print("-" * 90)
    print(f"\nCore - Rim dose difference (mean Gy):")
    print(f"  Mean:   {np.mean(core_rim_diffs):+.3f} Gy")
    print(f"  Median: {np.median(core_rim_diffs):+.3f} Gy")
    print(f"  Range:  [{min(core_rim_diffs):+.3f}, {max(core_rim_diffs):+.3f}] Gy")

    print(f"\nCT-High - CT-Low dose difference (mean Gy):")
    print(f"  Mean:   {np.mean(ctlow_cthigh_diffs):+.3f} Gy")
    print(f"  Median: {np.median(ctlow_cthigh_diffs):+.3f} Gy")
    print(f"  Range:  [{min(ctlow_cthigh_diffs):+.3f}, {max(ctlow_cthigh_diffs):+.3f}] Gy")

    print(f"\nInterpretation:")
    mean_diff = np.mean(core_rim_diffs)
    if abs(mean_diff) < 0.5:
        print(f"  Standard plans deliver near-uniform dose across Core vs Rim (diff {mean_diff:+.2f} Gy).")
        print(f"  BioSCART dose painting would represent a departure from current practice.")
    else:
        direction = "higher" if mean_diff > 0 else "lower"
        print(f"  Standard plans already deliver {direction} dose to Core vs Rim ({mean_diff:+.2f} Gy).")

    ct_diff = np.mean(ctlow_cthigh_diffs)
    if abs(ct_diff) < 0.5:
        print(f"  Standard plans do NOT differentiate by CT intensity (diff {ct_diff:+.2f} Gy).")
    else:
        print(f"  Standard plans show {ct_diff:+.2f} Gy difference between CT-High and CT-Low regions.")

    print(f"\nFull data: {long_csv}")
    return long_csv


if __name__ == "__main__":
    output_root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(r"C:\Users\wya245\SaMD\bioscart_output")
    aggregate(output_root)
