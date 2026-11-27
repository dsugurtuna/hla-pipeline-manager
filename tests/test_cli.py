"""Tests for the command-line interface, run against examples/."""

from pathlib import Path

import pytest

import hla_pipeline
from hla_pipeline.__main__ import main

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def test_plan_prints_commands(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(
        [
            "plan",
            "--batch-id",
            "demo",
            "--bfile",
            str(EXAMPLES / "synthetic"),
            "--work-dir",
            str(tmp_path),
            "--reference",
            "HM_CEU_REF",
            "--sub-batch-size",
            "5",
        ]
    )
    out = capsys.readouterr().out
    assert code == 0
    assert "12 samples -> 3 sub-batches" in out
    assert out.count("SNP2HLA.csh") == 3


def test_verify_flags_incomplete_example(capsys: pytest.CaptureFixture[str]) -> None:
    code = main(["verify", str(EXAMPLES / "batch"), "--expected", "2"])
    out = capsys.readouterr().out
    assert code == 1
    assert "sub_batch_001_imputed: OK" in out
    assert "sub_batch_002_imputed: INCOMPLETE" in out
    assert "- missing or empty .bgl.log" in out


def test_report_snp2hla(capsys: pytest.CaptureFixture[str]) -> None:
    batch = EXAMPLES / "batch"
    code = main(
        [
            "report",
            str(batch / "sub_batch_001_imputed.dosage"),
            "--fam",
            str(batch / "sub_batch_001_imputed.fam"),
            "--allele",
            "HLA_DRB1_0401",
            "--r2",
            str(batch / "sub_batch_001_imputed.bgl.r2"),
        ]
    )
    captured = capsys.readouterr()
    assert code == 0
    assert "SYN03\t1.968\thomozygous  r2=0.96" in captured.out
    assert "2 carriers of 4 samples" in captured.err


def test_report_raw_to_csv(tmp_path: Path) -> None:
    out = tmp_path / "carriers.csv"
    main(
        [
            "report",
            str(EXAMPLES / "synthetic_hla.raw"),
            "--allele",
            "HLA_DRB1_1501",
            "--out",
            str(out),
        ]
    )
    rows = out.read_text().splitlines()
    assert rows[0] == "participant_id,allele,dosage,call,r2_score"
    assert len(rows) == 13


def test_plan_cli_does_not_create_work_dir(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    work = tmp_path / "nonexistent_work"
    code = main(
        [
            "plan",
            "--batch-id",
            "demo",
            "--bfile",
            str(EXAMPLES / "synthetic"),
            "--work-dir",
            str(work),
            "--reference",
            "HM_CEU_REF",
            "--sub-batch-size",
            "5",
        ]
    )
    capsys.readouterr()
    assert code == 0
    assert not work.exists()


def test_verify_empty_dir_returns_nonzero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(["verify", str(tmp_path)])
    out = capsys.readouterr().out
    assert code == 1
    assert "NO OUTPUTS" in out


def test_verify_missing_dir_returns_nonzero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(["verify", str(tmp_path / "missing")])
    captured = capsys.readouterr()
    assert code == 1
    assert "error:" in captured.err


def test_verify_mismatched_samples_size_returns_nonzero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(["verify", str(tmp_path), "--samples", "10"])
    assert code == 1
    assert "error:" in capsys.readouterr().err


def test_plan_invalid_batch_id_returns_nonzero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code = main(
        [
            "plan",
            "--batch-id",
            "../bad",
            "--bfile",
            str(EXAMPLES / "synthetic"),
            "--work-dir",
            str(tmp_path),
            "--reference",
            "HM_CEU_REF",
        ]
    )
    assert code == 1
    assert "error:" in capsys.readouterr().err


def test_package_version() -> None:
    assert hla_pipeline.__version__ == "0.1.0"
