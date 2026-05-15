"""Tests for bioscart_minimal.recommendations."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from bioscart_minimal.recommendations import build_recommendations, write_recommendations
from bioscart_minimal.regions import build_five_gtv_regions


class TestBuildRecommendations:

    def _build(self, synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        return build_recommendations(
            ct=synthetic_ct_image,
            gtv=synthetic_gtv,
            regions=regions,
            radiomics_available=False,
        )

    def test_research_only_status(self, synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        recs = self._build(synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        assert recs["clinical_status"] == "research_only_not_for_treatment_decisions"

    def test_all_regions_present(self, synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        recs = self._build(synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        region_names = {r["region"] for r in recs["region_recommendations"]}
        assert "BSCART_Rim_5mm" in region_names
        assert "BSCART_Core" in region_names
        assert "BSCART_CT_Low_Q25" in region_names
        assert "BSCART_CT_High_Q75" in region_names
        assert "BSCART_Texture_High_Q75" in region_names

    def test_warnings_propagated(self, synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        synthetic_ct_image.warnings.append("test warning from CT")
        synthetic_gtv.warnings.append("test warning from GTV")
        recs = self._build(synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        assert "test warning from CT" in recs["warnings"]
        assert "test warning from GTV" in recs["warnings"]

    def test_volume_adequate_flag_true(self, synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        recs = self._build(synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        assert recs["case_summary"]["gtv_volume_adequate_for_radiomics"] is True

    def test_volume_adequate_flag_false(self, synthetic_ct_image, small_gtv, synthetic_ct_volume, small_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, small_gtv_mask, synthetic_spacing_zyx)
        recs = build_recommendations(
            ct=synthetic_ct_image,
            gtv=small_gtv,
            regions=regions,
            radiomics_available=False,
        )
        assert recs["case_summary"]["gtv_volume_adequate_for_radiomics"] is False


class TestWriteRecommendations:

    def test_json_output_valid(self, tmp_path, synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        recs = build_recommendations(synthetic_ct_image, synthetic_gtv, regions, False)
        json_path, md_path = write_recommendations(recs, tmp_path)
        data = json.loads(Path(json_path).read_text(encoding="utf-8"))
        assert "clinical_status" in data
        assert "region_recommendations" in data

    def test_markdown_has_header(self, tmp_path, synthetic_ct_image, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        recs = build_recommendations(synthetic_ct_image, synthetic_gtv, regions, False)
        _, md_path = write_recommendations(recs, tmp_path)
        text = Path(md_path).read_text(encoding="utf-8")
        assert text.startswith("# BioSCART")
