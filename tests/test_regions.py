"""Tests for bioscart_minimal.regions."""

from __future__ import annotations

import numpy as np
import pytest

from bioscart_minimal.regions import (
    MIN_GTV_VOLUME_CC_FOR_RADIOMICS,
    RegionMask,
    build_five_gtv_regions,
    check_gtv_volume_for_radiomics,
    local_std,
)


EXPECTED_REGION_NAMES = {
    "BSCART_Rim_5mm",
    "BSCART_Core",
    "BSCART_CT_Low_Q25",
    "BSCART_CT_High_Q75",
    "BSCART_Texture_High_Q75",
}


class TestBuildFiveGTVRegions:

    def test_returns_five_regions(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        assert len(regions) == 5

    def test_all_region_names_correct(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        assert set(regions.keys()) == EXPECTED_REGION_NAMES

    def test_all_regions_within_gtv(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        for name, region in regions.items():
            outside = region.mask & ~synthetic_gtv_mask
            assert not np.any(outside), f"{name} has voxels outside GTV"

    def test_rim_core_cover_gtv(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        combined = regions["BSCART_Rim_5mm"].mask | regions["BSCART_Core"].mask
        assert np.array_equal(combined, synthetic_gtv_mask)

    def test_rim_core_minimal_overlap(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        overlap = regions["BSCART_Rim_5mm"].mask & regions["BSCART_Core"].mask
        assert not np.any(overlap), "Rim and core should not overlap"

    def test_ct_low_threshold(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        q25 = float(np.percentile(synthetic_ct_volume[synthetic_gtv_mask], 25))
        low_mask = regions["BSCART_CT_Low_Q25"].mask
        assert np.all(synthetic_ct_volume[low_mask] <= q25 + 1e-6)

    def test_ct_high_threshold(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        q75 = float(np.percentile(synthetic_ct_volume[synthetic_gtv_mask], 75))
        high_mask = regions["BSCART_CT_High_Q75"].mask
        assert np.all(synthetic_ct_volume[high_mask] >= q75 - 1e-6)

    def test_small_gtv_fallback_core(self, synthetic_ct_volume, small_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, small_gtv_mask, synthetic_spacing_zyx)
        assert np.any(regions["BSCART_Core"].mask), "Core must exist even for small GTV"

    def test_shape_mismatch_raises(self, synthetic_gtv_mask, synthetic_spacing_zyx):
        bad_ct = np.zeros((10, 10, 10), dtype=np.float32)
        with pytest.raises(ValueError, match="shapes do not match"):
            build_five_gtv_regions(bad_ct, synthetic_gtv_mask, synthetic_spacing_zyx)

    def test_empty_mask_raises(self, synthetic_ct_volume, synthetic_spacing_zyx):
        empty = np.zeros(synthetic_ct_volume.shape, dtype=bool)
        with pytest.raises(ValueError, match="empty"):
            build_five_gtv_regions(synthetic_ct_volume, empty, synthetic_spacing_zyx)

    def test_region_volumes_positive(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        for name, region in regions.items():
            assert region.volume_cc > 0, f"{name} has zero volume"

    def test_to_dict_serialization(self, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        for region in regions.values():
            d = region.to_dict()
            assert "name" in d
            assert "role" in d
            assert "evidence_level" in d
            assert "volume_cc" in d
            assert "voxel_count" in d
            assert isinstance(d["notes"], list)


class TestCheckGTVVolumeForRadiomics:

    def test_adequate_volume_no_warnings(self):
        assert check_gtv_volume_for_radiomics(25.0) == []

    def test_small_volume_produces_warning(self):
        warnings = check_gtv_volume_for_radiomics(3.0)
        assert len(warnings) == 1
        assert "below" in warnings[0]

    def test_custom_threshold(self):
        assert check_gtv_volume_for_radiomics(8.0, min_cc=5.0) == []
        assert len(check_gtv_volume_for_radiomics(3.0, min_cc=5.0)) == 1

    def test_boundary_value(self):
        assert check_gtv_volume_for_radiomics(MIN_GTV_VOLUME_CC_FOR_RADIOMICS) == []
        assert len(check_gtv_volume_for_radiomics(MIN_GTV_VOLUME_CC_FOR_RADIOMICS - 0.01)) == 1


class TestLocalStd:

    def test_constant_image_zero_std(self):
        img = np.ones((10, 10, 10), dtype=np.float32) * 42.0
        result = local_std(img, size=3)
        assert np.allclose(result, 0.0, atol=1e-5)

    def test_varying_image_nonzero_std(self):
        rng = np.random.default_rng(0)
        img = rng.standard_normal((20, 20, 20)).astype(np.float32) * 100
        result = local_std(img, size=5)
        assert np.any(result > 0)
