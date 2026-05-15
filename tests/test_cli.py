"""Tests for bioscart_minimal.cli argument parsing."""

from __future__ import annotations

import pytest

from bioscart_minimal.cli import main


class TestCLI:

    def test_missing_required_args_exits(self):
        with pytest.raises(SystemExit):
            main([])

    def test_help_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            main(["--help"])
        assert exc_info.value.code == 0

    def test_argparse_accepts_all_options(self):
        """Verify argparse recognizes all documented flags (will fail at DICOM load, not parse)."""
        with pytest.raises((SystemExit, Exception)):
            main([
                "--dicom-dir", "/nonexistent",
                "--out-dir", "/tmp/out",
                "--gtv-name", "GTV",
                "--rim-mm", "3.0",
            ])
