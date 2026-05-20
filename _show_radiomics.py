"""Print key radiomics features per BioSCART region for one case."""
import csv, sys
from pathlib import Path

csv_path = Path(sys.argv[1])
with csv_path.open(encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

keys = [
    "original_firstorder_Mean",
    "original_firstorder_Median",
    "original_firstorder_Entropy",
    "original_firstorder_Kurtosis",
    "original_firstorder_Skewness",
    "original_firstorder_Variance",
    "original_firstorder_10Percentile",
    "original_firstorder_90Percentile",
    "original_glcm_Contrast",
    "original_glcm_Correlation",
    "original_glcm_JointEntropy",
    "original_glcm_ClusterShade",
    "original_ngtdm_Coarseness",
    "original_ngtdm_Complexity",
    "original_ngtdm_Contrast",
    "original_shape_Sphericity",
    "original_shape_Elongation",
    "original_shape_MeshVolume",
]

labels = []
for r in rows:
    lab = r["region"].replace("BSCART_","B:").replace("_Q25","<").replace("_Q75",">").replace("_5mm","")
    labels.append(lab)

print(f"{'Feature':<32}", end="")
for lab in labels:
    print(f"{lab:>14}", end="")
print()
print("-" * (32 + 14 * len(rows)))

for k in keys:
    short = k.replace("original_","").replace("firstorder_","1st:").replace("glcm_","glcm:").replace("ngtdm_","ngtdm:").replace("shape_","shp:")
    print(f"{short:<32}", end="")
    for r in rows:
        val = r.get(k, "")
        try:
            v = float(val)
            if abs(v) > 10000:
                print(f"{v:>14.0f}", end="")
            elif abs(v) > 10:
                print(f"{v:>14.2f}", end="")
            elif abs(v) > 0.01:
                print(f"{v:>14.4f}", end="")
            else:
                print(f"{v:>14.6f}", end="")
        except:
            print(f"{str(val):>14}", end="")
    print()
