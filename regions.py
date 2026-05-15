"""GTV subregion generation for minimal BioSCART."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage


@dataclass
class RegionMask:
    name: str
    mask: np.ndarray
    role: str
    evidence_level: str
    volume_cc: float
    notes: list[str]

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "role": self.role,
            "evidence_level": self.evidence_level,
            "volume_cc": round(self.volume_cc, 3),
            "voxel_count": int(self.mask.sum()),
            "notes": self.notes,
        }


MIN_GTV_VOLUME_CC_FOR_RADIOMICS = 10.0


def check_gtv_volume_for_radiomics(
    gtv_volume_cc: float,
    min_cc: float = MIN_GTV_VOLUME_CC_FOR_RADIOMICS,
) -> list[str]:
    """Return warnings if GTV volume is too small for reliable subregion radiomics."""
    if gtv_volume_cc < min_cc:
        return [
            f"GTV volume ({gtv_volume_cc:.1f} cc) is below the {min_cc:.0f} cc minimum "
            f"for reliable five-region radiomics. Features may be statistically unreliable."
        ]
    return []


def build_five_gtv_regions(
    ct_hu: np.ndarray,
    gtv_mask: np.ndarray,
    spacing_zyx_mm: tuple[float, float, float],
    rim_mm: float = 5.0,
) -> dict[str, RegionMask]:
    """Create five transparent, CT-only research regions inside the GTV.

    The masks are allowed to overlap. They are candidate biologic habitats, not
    validated prescription structures.
    """
    if ct_hu.shape != gtv_mask.shape:
        raise ValueError("CT image and GTV mask shapes do not match")
    if not np.any(gtv_mask):
        raise ValueError("GTV mask is empty")

    voxel_cc = float(np.prod(spacing_zyx_mm) / 1000.0)
    gtv_values = ct_hu[gtv_mask]
    q25 = float(np.percentile(gtv_values, 25))
    q75 = float(np.percentile(gtv_values, 75))

    distance_inside = ndimage.distance_transform_edt(gtv_mask, sampling=spacing_zyx_mm)
    rim = gtv_mask & (distance_inside <= rim_mm)
    core = gtv_mask & (distance_inside > rim_mm)
    if not np.any(core):
        core = gtv_mask & (distance_inside >= np.percentile(distance_inside[gtv_mask], 60))

    ct_low = gtv_mask & (ct_hu <= q25)
    ct_high = gtv_mask & (ct_hu >= q75)

    texture = local_std(ct_hu, size=5)
    texture_q75 = float(np.percentile(texture[gtv_mask], 75))
    texture_high = gtv_mask & (texture >= texture_q75)

    return {
        "BSCART_Rim_5mm": RegionMask(
            name="BSCART_Rim_5mm",
            mask=rim,
            role="border_or_valley_candidate",
            evidence_level="geometric_ct_only",
            volume_cc=float(rim.sum() * voxel_cc),
            notes=[
                "Outer GTV shell. Candidate immune-access/transition region for low-to-intermediate dose review.",
                "Not an immune-excluded label without pathology/IHC validation.",
            ],
        ),
        "BSCART_Core": RegionMask(
            name="BSCART_Core",
            mask=core,
            role="peak_candidate",
            evidence_level="geometric_ct_only",
            volume_cc=float(core.sum() * voxel_cc),
            notes=[
                "Inner GTV region after rim contraction. Candidate high-dose vertex review if OAR geometry allows.",
            ],
        ),
        "BSCART_CT_Low_Q25": RegionMask(
            name="BSCART_CT_Low_Q25",
            mask=ct_low,
            role="necrosis_or_hypoxia_surrogate_review",
            evidence_level="ct_intensity_surrogate",
            volume_cc=float(ct_low.sum() * voxel_cc),
            notes=[
                f"Lowest CT-intensity quartile inside GTV (threshold <= {q25:.1f} HU).",
                "May represent necrosis, low cellularity, artifact, or hypoxia surrogate depending on tumor/site/contrast.",
                "Do not call this confirmed hypoxia without hypoxia PET, MRI, or pathology support.",
            ],
        ),
        "BSCART_CT_High_Q75": RegionMask(
            name="BSCART_CT_High_Q75",
            mask=ct_high,
            role="viable_or_enhancing_surrogate_review",
            evidence_level="ct_intensity_surrogate",
            volume_cc=float(ct_high.sum() * voxel_cc),
            notes=[
                f"Highest CT-intensity quartile inside GTV (threshold >= {q75:.1f} HU).",
                "Candidate viable/enhancing/cellular region depending on acquisition and contrast phase.",
            ],
        ),
        "BSCART_Texture_High_Q75": RegionMask(
            name="BSCART_Texture_High_Q75",
            mask=texture_high,
            role="heterogeneity_candidate",
            evidence_level="ct_texture_surrogate",
            volume_cc=float(texture_high.sum() * voxel_cc),
            notes=[
                f"Highest local CT heterogeneity quartile inside GTV (local std threshold >= {texture_q75:.1f} HU).",
                "Candidate resistant/heterogeneous habitat for review; validate against pathology before escalation claims.",
            ],
        ),
    }


def local_std(image: np.ndarray, size: int = 5) -> np.ndarray:
    """Compute local standard deviation with a cubic window."""
    image = image.astype(np.float32)
    mean = ndimage.uniform_filter(image, size=size)
    mean_sq = ndimage.uniform_filter(image * image, size=size)
    var = np.maximum(mean_sq - mean * mean, 0.0)
    return np.sqrt(var)

