"""Tests for the SFRT/SCART SIB dose prescription logic."""

import json
import numpy as np
import pytest

from bioscart_minimal.sfrt_dose_logic import (
    PROTOCOLS,
    SCARTDoseLevel,
    SCARTProtocol,
    assign_doses_to_regions,
    build_sfrt_plan,
    check_immune_safety,
    generate_sib_objectives,
    scart_cervix_protocol,
    scart_esophagus_protocol,
    scart_lung_protocol,
    scart_single_fraction_protocol,
    write_sfrt_plan,
)
from bioscart_minimal.regions import RegionMask


# ── Fixtures ────────────────────────────────────────────────────────

@pytest.fixture
def mock_regions():
    """Five BioSCART regions with realistic roles."""
    mask = np.ones((3, 3, 3), dtype=bool)
    return {
        "BSCART_Core": RegionMask(
            "BSCART_Core", mask, "peak_candidate",
            "geometric_ct_only", 5.0, ["core"],
        ),
        "BSCART_CT_High_Q75": RegionMask(
            "BSCART_CT_High_Q75", mask, "viable_or_enhancing_surrogate_review",
            "ct_intensity_surrogate", 4.0, ["high"],
        ),
        "BSCART_Texture_High_Q75": RegionMask(
            "BSCART_Texture_High_Q75", mask, "heterogeneity_candidate",
            "ct_texture_surrogate", 3.0, ["texture"],
        ),
        "BSCART_Rim_5mm": RegionMask(
            "BSCART_Rim_5mm", mask, "border_or_valley_candidate",
            "geometric_ct_only", 8.0, ["rim"],
        ),
        "BSCART_CT_Low_Q25": RegionMask(
            "BSCART_CT_Low_Q25", mask, "necrosis_or_hypoxia_surrogate_review",
            "ct_intensity_surrogate", 3.0, ["low"],
        ),
    }


# ── Protocol Tests ──────────────────────────────────────────────────

class TestProtocols:
    def test_lung_protocol_defaults(self):
        p = scart_lung_protocol(5)
        assert p.protocol_id == "SCART-LUNG-v1"
        assert p.site == "lung"
        assert p.fractions == 5
        assert "peak" in p.levels
        assert "intermediate" in p.levels
        assert "valley" in p.levels

    def test_lung_protocol_doses(self):
        p = scart_lung_protocol(5)
        assert p.levels["peak"].total_dose_gy == 66.7
        assert p.levels["valley"].total_dose_gy == 20.0
        assert p.pvdr_target == 5.0

    def test_esophagus_protocol_lumen_safe(self):
        p = scart_esophagus_protocol(28)
        assert p.protocol_id == "SCART-ESOPH-v2"
        assert p.levels["peak"].total_dose_gy == 58.8
        assert p.levels["peak"].dose_per_fraction_gy == 2.1
        assert p.levels["intermediate"].total_dose_gy == 50.4
        assert p.levels["valley"].total_dose_gy == 44.8
        assert p.levels["valley"].role == "lumen_sparing_immune_sanctuary"

    def test_esophagus_peak_below_artdeco(self):
        """Peak dose must be below ARTDECO failure threshold of 61.6 Gy."""
        p = scart_esophagus_protocol(28)
        assert p.levels["peak"].total_dose_gy < 61.6

    def test_cervix_protocol_defaults(self):
        p = scart_cervix_protocol(28)
        assert p.protocol_id == "SCART-CERVIX-v1"
        assert p.site == "cervix"
        assert p.fractions == 28

    def test_cervix_protocol_doses(self):
        p = scart_cervix_protocol(28)
        assert p.levels["peak"].dose_per_fraction_gy == 2.3
        assert p.levels["peak"].total_dose_gy == 64.4
        assert p.levels["intermediate"].total_dose_gy == 50.4
        assert p.levels["valley"].total_dose_gy == 44.8

    def test_cervix_peak_higher_than_esophagus(self):
        """Cervix tolerates higher boost -- no thin lumen risk."""
        c = scart_cervix_protocol(28)
        e = scart_esophagus_protocol(28)
        assert c.levels["peak"].dose_per_fraction_gy > e.levels["peak"].dose_per_fraction_gy

    def test_cervix_valley_protects_hollow_organs(self):
        p = scart_cervix_protocol(28)
        assert "hollow_organ" in p.levels["valley"].role

    def test_cervix_in_registry(self):
        assert "cervix" in PROTOCOLS

    def test_single_fraction_protocol(self):
        p = scart_single_fraction_protocol()
        assert p.fractions == 1
        assert p.levels["peak"].total_dose_gy == 15.0
        assert p.levels["valley"].total_dose_gy == 3.0

    def test_protocols_registry(self):
        assert "lung" in PROTOCOLS
        assert "esophagus" in PROTOCOLS
        assert "single_fraction" in PROTOCOLS

    def test_all_protocols_have_three_levels(self):
        for name, factory in PROTOCOLS.items():
            p = factory()
            assert len(p.levels) == 3, f"{name} protocol missing levels"
            assert "peak" in p.levels
            assert "intermediate" in p.levels
            assert "valley" in p.levels

    def test_protocol_to_dict(self):
        p = scart_lung_protocol(5)
        d = p.to_dict()
        assert d["protocol_id"] == "SCART-LUNG-v1"
        assert "peak" in d["levels"]
        assert d["levels"]["peak"]["total_dose_gy"] == 66.7

    def test_pvdr_ordering(self):
        """Peak dose must always exceed valley dose."""
        for name, factory in PROTOCOLS.items():
            p = factory()
            assert p.levels["peak"].total_dose_gy > p.levels["valley"].total_dose_gy


# ── Region Assignment Tests ─────────────────────────────────────────

class TestRegionAssignment:
    def test_assigns_all_regions(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        assert len(assignments) == 5

    def test_core_gets_peak(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        core = [a for a in assignments if a.region_name == "BSCART_Core"][0]
        assert core.dose_level == "peak"

    def test_rim_gets_valley(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        rim = [a for a in assignments if a.region_name == "BSCART_Rim_5mm"][0]
        assert rim.dose_level == "valley"

    def test_texture_gets_intermediate(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        texture = [a for a in assignments if a.region_name == "BSCART_Texture_High_Q75"][0]
        assert texture.dose_level == "intermediate"

    def test_low_density_gets_valley(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        low = [a for a in assignments if a.region_name == "BSCART_CT_Low_Q25"][0]
        assert low.dose_level == "valley"

    def test_high_density_gets_peak(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        high = [a for a in assignments if a.region_name == "BSCART_CT_High_Q75"][0]
        assert high.dose_level == "peak"

    def test_esophagus_rim_is_lumen_safe(self, mock_regions):
        """Rim region must get valley dose in esophageal protocol (lumen protection)."""
        p = scart_esophagus_protocol(28)
        assignments = assign_doses_to_regions(mock_regions, p)
        rim = [a for a in assignments if a.region_name == "BSCART_Rim_5mm"][0]
        assert rim.dose_level == "valley"
        assert rim.total_dose_gy == 44.8

    def test_cervix_rim_protects_hollow_organs(self, mock_regions):
        """Rim region must get valley dose in cervix protocol (bladder/rectum protection)."""
        p = scart_cervix_protocol(28)
        assignments = assign_doses_to_regions(mock_regions, p)
        rim = [a for a in assignments if a.region_name == "BSCART_Rim_5mm"][0]
        assert rim.dose_level == "valley"
        assert rim.total_dose_gy == 44.8

    def test_cervix_core_gets_higher_peak_than_esophagus(self, mock_regions):
        c = scart_cervix_protocol(28)
        e = scart_esophagus_protocol(28)
        c_assign = assign_doses_to_regions(mock_regions, c)
        e_assign = assign_doses_to_regions(mock_regions, e)
        c_core = [a for a in c_assign if a.region_name == "BSCART_Core"][0]
        e_core = [a for a in e_assign if a.region_name == "BSCART_Core"][0]
        assert c_core.total_dose_gy > e_core.total_dose_gy


# ── SIB Objectives Tests ────────────────────────────────────────────

class TestSIBObjectives:
    def test_generates_objectives(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        objectives = generate_sib_objectives(assignments, p)
        assert len(objectives) > 0

    def test_has_lower_and_upper_per_region(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        objectives = generate_sib_objectives(assignments, p)
        region_objs = [o for o in objectives if o.structure == "BSCART_Core"]
        types = {o.objective_type for o in region_objs}
        assert "Lower" in types
        assert "Upper" in types

    def test_includes_oar_constraints(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        objectives = generate_sib_objectives(assignments, p)
        structures = {o.structure for o in objectives}
        assert "Cord" in structures
        assert "Total_Lung" in structures
        assert "Heart" in structures

    def test_includes_nto(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        objectives = generate_sib_objectives(assignments, p)
        nto = [o for o in objectives if o.structure == "NTO"]
        assert len(nto) == 1
        assert nto[0].objective_type == "Automatic"
        assert nto[0].priority == 100

    def test_peak_priority_highest(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        objectives = generate_sib_objectives(assignments, p)
        core_objs = [o for o in objectives if o.structure == "BSCART_Core"]
        for o in core_objs:
            assert o.priority == 100

    def test_custom_oar_constraints(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        custom = {"Cord": [{"type": "Upper", "vol": 0.0, "dose_cgy": 3000, "priority": 90}]}
        objectives = generate_sib_objectives(assignments, p, oar_constraints=custom)
        cord = [o for o in objectives if o.structure == "Cord"]
        assert len(cord) == 1
        assert cord[0].dose_cgy == 3000

    def test_to_dict(self, mock_regions):
        p = scart_lung_protocol(5)
        assignments = assign_doses_to_regions(mock_regions, p)
        objectives = generate_sib_objectives(assignments, p)
        d = objectives[0].to_dict()
        assert "structure" in d
        assert "dose_cgy" in d
        assert "priority" in d


# ── Immune Safety Tests ─────────────────────────────────────────────

class TestImmuneSafety:
    def test_no_ice3_returns_safe(self):
        p = scart_lung_protocol(5)
        result = check_immune_safety(None, p)
        assert result.safe is True
        assert result.ice3_available is False

    def test_safe_nadir(self):
        p = scart_lung_protocol(5)
        metrics = {"alc_predicted_nadir": 0.7, "mean_blood_dose_gy": 0.1}
        result = check_immune_safety(metrics, p)
        assert result.safe is True
        assert result.ice3_available is True
        assert "IMMUNE SAFE" in result.recommendation

    def test_unsafe_nadir(self):
        p = scart_lung_protocol(5)
        metrics = {"alc_predicted_nadir": 0.3, "mean_blood_dose_gy": 0.5}
        result = check_immune_safety(metrics, p)
        assert result.safe is False
        assert "IMMUNE RISK" in result.recommendation

    def test_missing_nadir_key(self):
        p = scart_lung_protocol(5)
        metrics = {"mean_blood_dose_gy": 0.2}
        result = check_immune_safety(metrics, p)
        assert result.ice3_available is True
        assert result.alc_nadir_predicted is None

    def test_to_dict(self):
        p = scart_lung_protocol(5)
        result = check_immune_safety(None, p)
        d = result.to_dict()
        assert "safe" in d
        assert "recommendation" in d


# ── Full Plan Tests ─────────────────────────────────────────────────

class TestBuildSFRTPlan:
    def test_builds_plan(self, mock_regions):
        p = scart_lung_protocol(5)
        plan = build_sfrt_plan(mock_regions, p)
        assert plan["schema_version"] == "bioscart.sfrt_plan.v0.1"
        assert "region_dose_assignments" in plan
        assert "eclipse_sib_objectives" in plan
        assert "dose_summary" in plan
        assert "immune_safety" in plan

    def test_pvdr_calculated(self, mock_regions):
        p = scart_lung_protocol(5)
        plan = build_sfrt_plan(mock_regions, p)
        assert plan["dose_summary"]["actual_pvdr"] is not None
        assert plan["dose_summary"]["actual_pvdr"] > 1.0

    def test_esophagus_pvdr_modest(self, mock_regions):
        p = scart_esophagus_protocol(28)
        plan = build_sfrt_plan(mock_regions, p)
        pvdr = plan["dose_summary"]["actual_pvdr"]
        assert pvdr < 2.0, "Esophageal PVDR should be modest for conventional fractionation"

    def test_cervix_plan(self, mock_regions):
        p = scart_cervix_protocol(28)
        plan = build_sfrt_plan(mock_regions, p)
        assert plan["protocol"]["protocol_id"] == "SCART-CERVIX-v1"
        assert plan["dose_summary"]["peak_total_gy"] == 64.4
        assert plan["dose_summary"]["valley_total_gy"] == 44.8
        assert plan["dose_summary"]["actual_pvdr"] is not None

    def test_write_sfrt_plan(self, mock_regions, tmp_path):
        p = scart_lung_protocol(5)
        plan = build_sfrt_plan(mock_regions, p)
        path = write_sfrt_plan(plan, tmp_path)
        assert path.exists()
        loaded = json.loads(path.read_text(encoding="utf-8"))
        assert loaded["protocol"]["protocol_id"] == "SCART-LUNG-v1"

    def test_plan_with_ice3_safe(self, mock_regions):
        p = scart_lung_protocol(5)
        ice3 = {"alc_predicted_nadir": 0.8}
        plan = build_sfrt_plan(mock_regions, p, ice3_metrics=ice3)
        assert plan["immune_safety"]["safe"] is True

    def test_plan_with_ice3_risk(self, mock_regions):
        p = scart_lung_protocol(5)
        ice3 = {"alc_predicted_nadir": 0.2}
        plan = build_sfrt_plan(mock_regions, p, ice3_metrics=ice3)
        assert plan["immune_safety"]["safe"] is False

    def test_plan_json_serializable(self, mock_regions):
        p = scart_lung_protocol(5)
        plan = build_sfrt_plan(mock_regions, p)
        serialized = json.dumps(plan)
        assert len(serialized) > 0

    def test_research_only_status(self, mock_regions):
        p = scart_lung_protocol(5)
        plan = build_sfrt_plan(mock_regions, p)
        assert plan["clinical_status"] == "research_only_not_for_treatment_decisions"
