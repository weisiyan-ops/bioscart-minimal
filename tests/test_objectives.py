"""Tests for BioSCART objective template generation."""

from __future__ import annotations

import json

from bioscart_minimal.objectives import build_objective_template, write_objective_template
from bioscart_minimal.regions import build_five_gtv_regions


class TestObjectives:

    def test_template_is_research_only(self, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        template = build_objective_template(synthetic_gtv, regions)
        assert template["clinical_status"] == "research_only_not_for_treatment_decisions"
        assert not template["objective_policy"]["prescription_doses_included"]
        assert len(template["objectives"]) == 6

    def test_write_template_json(self, tmp_path, synthetic_gtv, synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx):
        regions = build_five_gtv_regions(synthetic_ct_volume, synthetic_gtv_mask, synthetic_spacing_zyx)
        template = build_objective_template(synthetic_gtv, regions)
        path = write_objective_template(template, tmp_path)
        data = json.loads(path.read_text(encoding="utf-8"))
        assert data["schema_version"] == "bioscart.objectives.v0.2"

