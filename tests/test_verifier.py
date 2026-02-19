"""Tests for hla_pipeline.verifier."""

from pathlib import Path

import pytest

from hla_pipeline.verifier import ImputationVerifier, SubBatchStatus


@pytest.fixture()
def complete_sub_batch(tmp_path: Path) -> Path:
    prefix = tmp_path / "sub_001"
    (tmp_path / "sub_001.bed").write_bytes(b"\x00")
    (tmp_path / "sub_001.bim").write_text("6\tHLA_DRB1_0101\t0\t32000000\tA\tG\n")
    (tmp_path / "sub_001.fam").write_text("S001 S001 0 0 1 -9\n")
    (tmp_path / "sub_001.dosage").write_text("header\n")
    (tmp_path / "sub_001.bgl.r2").write_text("HLA_DRB1_0101 0.95\n")
    (tmp_path / "sub_001.bgl.log").write_text("Analysis finished successfully.\n")
    return prefix


class TestImputationVerifier:
    def test_complete_sub_batch(self, complete_sub_batch: Path) -> None:
        verifier = ImputationVerifier()
        status = verifier.verify_sub_batch(complete_sub_batch)
        assert status.is_complete
        assert status.hla_marker_count == 1
        assert status.beagle_completed

    def test_incomplete_sub_batch(self, tmp_path: Path) -> None:
        prefix = tmp_path / "sub_999"
        (tmp_path / "sub_999.fam").write_text("S001 S001 0 0 1 -9\n")
        verifier = ImputationVerifier()
        status = verifier.verify_sub_batch(prefix)
        assert not status.is_complete

    def test_verify_batch(self, complete_sub_batch: Path) -> None:
        verifier = ImputationVerifier()
        report = verifier.verify_batch(complete_sub_batch.parent)
        assert report.total_sub_batches >= 1
        text = verifier.format_report(report)
        assert "OK" in text or "INCOMPLETE" in text
