"""Tests for the command-line interface, run against examples/."""

from pathlib import Path

import pytest

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
