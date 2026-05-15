"""Tests for bioscart_minimal.radiomics_runner."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest import mock

import numpy as np
import pytest
import SimpleITK as sitk

from bioscart_minimal.radiomics_runner import RadiomicsRunResult, run_pyradiomics
from bioscart_minimal.regions import RegionMask


def _make_dummy_regions(shape):
    mask = np.zeros(shape, dtype=bool)
    mask[10:20, 30:50, 30:50] = True
    return {
        "BSCART_Test": RegionMask(
            name="BSCART_Test",
            mask=mask,
            role="test",
            evidence_level="test",
            volume_cc=1.0,
            notes=[],
        )
    }


class TestRunPyradiomics:

    def test_graceful_fallback_when_not_installed(self, tmp_path, synthetic_ct_image, synthetic_gtv_mask):
        with mock.patch.dict(sys.modules, {"radiomics": None, "radiomics.featureextractor": None}):
            import importlib
            import bioscart_minimal.radiomics_runner as mod
            importlib.reload(mod)

            regions = _make_dummy_regions(synthetic_ct_image.array_hu.shape)
            result = mod.run_pyradiomics(
                image=synthetic_ct_image.sitk_image,
                gtv_mask=synthetic_gtv_mask,
                regions=regions,
                output_dir=tmp_path / "radiomics",
            )
            assert not result.available
            assert "not installed" in result.message.lower()

            importlib.reload(mod)

    def test_status_json_written_when_unavailable(self, tmp_path, synthetic_ct_image, synthetic_gtv_mask):
        with mock.patch.dict(sys.modules, {"radiomics": None, "radiomics.featureextractor": None}):
            import importlib
            import bioscart_minimal.radiomics_runner as mod
            importlib.reload(mod)

            regions = _make_dummy_regions(synthetic_ct_image.array_hu.shape)
            result = mod.run_pyradiomics(
                image=synthetic_ct_image.sitk_image,
                gtv_mask=synthetic_gtv_mask,
                regions=regions,
                output_dir=tmp_path / "radiomics",
            )
            status_path = Path(result.status_json)
            assert status_path.exists()
            status = json.loads(status_path.read_text(encoding="utf-8"))
            assert status["available"] is False

            importlib.reload(mod)

    def test_output_directory_created(self, tmp_path, synthetic_ct_image, synthetic_gtv_mask):
        out = tmp_path / "new_radiomics_dir"
        assert not out.exists()

        with mock.patch.dict(sys.modules, {"radiomics": None, "radiomics.featureextractor": None}):
            import importlib
            import bioscart_minimal.radiomics_runner as mod
            importlib.reload(mod)

            regions = _make_dummy_regions(synthetic_ct_image.array_hu.shape)
            mod.run_pyradiomics(
                image=synthetic_ct_image.sitk_image,
                gtv_mask=synthetic_gtv_mask,
                regions=regions,
                output_dir=out,
            )
            assert out.exists()

            importlib.reload(mod)
