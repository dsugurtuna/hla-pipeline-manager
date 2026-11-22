"""Tests for hla_pipeline.verifier."""

from pathlib import Path

import pytest

from hla_pipeline.verifier import ImputationVerifier


def _write_outputs(prefix: Path, *, log: str = "Analysis finished.\n") -> None:
    Path(f"{prefix}.bed").write_bytes(b"\x6c\x1b\x01\x00")
    Path(f"{prefix}.bim").write_text("6\tHLA_DRB1_0101\t0\t32000000\tP\tA\n")
    Path(f"{prefix}.fam").write_text("S001 S001 0 0 1 -9\n")
    Path(f"{prefix}.dosage").write_text("HLA_DRB1_0101\tP\tA\t1.000\n")
    Path(f"{prefix}.bgl.r2").write_text("HLA_DRB1_0101 0.95\n")
    Path(f"{prefix}.bgl.log").write_text(log)


@pytest.fixture()
def complete_sub_batch(tmp_path: Path) -> Path:
    prefix = tmp_path / "sub_001"
    _write_outputs(prefix)
    return prefix


def test_complete_sub_batch(complete_sub_batch: Path) -> None:
    status = ImputationVerifier().verify_sub_batch(complete_sub_batch)
    assert status.is_complete
    assert status.hla_marker_count == 1
    assert status.hla_r2_count == 1
    assert status.beagle_completed


def test_incomplete_sub_batch(tmp_path: Path) -> None:
    (tmp_path / "sub_999.fam").write_text("S001 S001 0 0 1 -9\n")
    assert not ImputationVerifier().verify_sub_batch(tmp_path / "sub_999").is_complete


def test_empty_file_counts_as_missing(complete_sub_batch: Path) -> None:
    Path(f"{complete_sub_batch}.dosage").write_text("")
    status = ImputationVerifier().verify_sub_batch(complete_sub_batch)
    assert not status.has_dosage
    assert not status.is_complete


def test_r2_without_hla_rows_is_incomplete(complete_sub_batch: Path) -> None:
    Path(f"{complete_sub_batch}.bgl.r2").write_text("rs1 0.99\n")
    assert not ImputationVerifier().verify_sub_batch(complete_sub_batch).is_complete


def test_log_without_completion_marker(tmp_path: Path) -> None:
    prefix = tmp_path / "sub_002"
    _write_outputs(prefix, log="java.lang.OutOfMemoryError\n")
    status = ImputationVerifier().verify_sub_batch(prefix)
    assert not status.beagle_completed
    assert status.problems == ["no completion line in Beagle log"]


def test_dotted_prefix(tmp_path: Path) -> None:
    prefix = tmp_path / "b28.sub_001"
    _write_outputs(prefix)
    assert ImputationVerifier().verify_sub_batch(prefix).is_complete


def test_verify_batch_ignores_pre_imputation_filesets(
    complete_sub_batch: Path,
) -> None:
    # A pre-imputation PLINK fileset in the same directory has no .dosage.
    (complete_sub_batch.parent / "sub_001_input.fam").write_text("S1 S1 0 0 1 -9\n")
    verifier = ImputationVerifier()
    report = verifier.verify_batch(complete_sub_batch.parent, batch_id="b1")
    assert report.total_sub_batches == 1
    assert report.all_ok
    assert "sub_001: OK" in verifier.format_report(report)


def test_expected_sub_batches_detects_missing(complete_sub_batch: Path) -> None:
    verifier = ImputationVerifier()
    report = verifier.verify_batch(
        complete_sub_batch.parent, batch_id="b1", expected_sub_batches=3
    )
    assert report.missing_sub_batches("b1") == 2
    assert not report.all_ok
    assert "MISSING: 2" in verifier.format_report(report)
