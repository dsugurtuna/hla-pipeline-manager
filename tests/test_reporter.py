"""Tests for hla_pipeline.reporter and hla_pipeline.dosage."""

from pathlib import Path

import pytest

from hla_pipeline.dosage import read_plink_raw, read_r2, read_snp2hla_dosage
from hla_pipeline.reporter import ClinicalReporter


@pytest.fixture()
def dosage_file(tmp_path: Path) -> Path:
    p = tmp_path / "dosage.raw"
    p.write_text(
        "FID\tIID\tHLA_DRB1_0101\nS001\tS001\t1.8\nS002\tS002\t0.9\nS003\tS003\t0.1\n"
    )
    return p


def test_genotype_calls(dosage_file: Path) -> None:
    report = ClinicalReporter().generate_report(dosage_file, "HLA_DRB1_0101")
    assert report.total_participants == 3
    assert report.carrier_count == 2
    calls = {g.participant_id: g.call for g in report.genotypes}
    assert calls == {"S001": "homozygous", "S002": "heterozygous", "S003": "negative"}


def test_export_csv(dosage_file: Path, tmp_path: Path) -> None:
    report = ClinicalReporter().generate_report(dosage_file, "HLA_DRB1_0101")
    out = tmp_path / "report.csv"
    ClinicalReporter.export_csv(report, out)
    content = out.read_text()
    assert content.startswith("participant_id,allele,dosage,call,r2_score")
    assert "homozygous" in content


def test_plink_raw_space_delimited_with_counted_allele(tmp_path: Path) -> None:
    raw = tmp_path / "x.raw"
    raw.write_text(
        "FID IID PAT MAT SEX PHENOTYPE HLA_DRB1_0401_P HLA_DRB1_1501_A\n"
        "F1 S1 0 0 1 -9 2 2\n"
        "F2 S2 0 0 2 -9 1 0\n"
        "F3 S3 0 0 2 -9 NA 1\n"
    )
    assert read_plink_raw(raw, "HLA_DRB1_0401") == {"S1": 2.0, "S2": 1.0}
    # Counted allele is A (absent): present-allele dosage is 2 - value.
    assert read_plink_raw(raw, "HLA_DRB1_1501") == {"S1": 0.0, "S2": 2.0, "S3": 1.0}
    report = ClinicalReporter().generate_report(raw, "HLA_DRB1_1501_A")
    assert report.allele == "HLA_DRB1_1501"
    assert report.carrier_count == 2


def test_non_hla_marker_is_not_flipped(tmp_path: Path) -> None:
    raw = tmp_path / "x.raw"
    raw.write_text("FID IID PAT MAT SEX PHENOTYPE rs123_A\nF1 S1 0 0 1 -9 2\n")
    assert read_plink_raw(raw, "rs123") == {"S1": 2.0}


def test_snp2hla_dosage_reader(tmp_path: Path) -> None:
    dosage = tmp_path / "out.dosage"
    fam = tmp_path / "out.fam"
    fam.write_text("F1 S1 0 0 1 -9\nF2 S2 0 0 2 -9\n")
    dosage.write_text(
        "HLA_DRB1_0401\tP\tA\t1.950\t0.020\nHLA_DRB1_1501\tA\tP\t0.000\t1.000\n"
    )
    assert read_snp2hla_dosage(dosage, fam, "HLA_DRB1_0401") == {
        "S1": 1.95,
        "S2": 0.02,
    }
    # allele1 is A here, so the present-allele dosage is flipped.
    assert read_snp2hla_dosage(dosage, fam, "HLA_DRB1_1501") == {"S1": 2.0, "S2": 1.0}
    with pytest.raises(KeyError):
        read_snp2hla_dosage(dosage, fam, "HLA_B_0801")


def test_snp2hla_dosage_sample_mismatch(tmp_path: Path) -> None:
    dosage = tmp_path / "out.dosage"
    fam = tmp_path / "out.fam"
    fam.write_text("F1 S1 0 0 1 -9\n")
    dosage.write_text("HLA_DRB1_0401\tP\tA\t1.0\t0.0\n")
    with pytest.raises(ValueError, match="2 dosage values but 1 samples"):
        read_snp2hla_dosage(dosage, fam, "HLA_DRB1_0401")


def test_r2_annotation(tmp_path: Path, dosage_file: Path) -> None:
    r2 = tmp_path / "x.bgl.r2"
    r2.write_text("HLA_DRB1_0101\t0.93\nHLA_B_0801\tNaN\n")
    scores = read_r2(r2)
    assert scores == {"HLA_DRB1_0101": 0.93}
    report = ClinicalReporter().generate_report(dosage_file, "HLA_DRB1_0101", scores)
    assert {g.r2_score for g in report.genotypes} == {0.93}


def test_invalid_thresholds() -> None:
    with pytest.raises(ValueError):
        ClinicalReporter(homozygous_threshold=0.4, heterozygous_threshold=0.5)
