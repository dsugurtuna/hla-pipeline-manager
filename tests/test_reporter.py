"""Tests for hla_pipeline.reporter."""

from pathlib import Path

import pytest

from hla_pipeline.reporter import ClinicalReporter, ClinicalReport


@pytest.fixture()
def dosage_file(tmp_path: Path) -> Path:
    p = tmp_path / "dosage.raw"
    p.write_text(
        "FID\tIID\tHLA_DRB1_0101\n"
        "S001\tS001\t1.8\n"
        "S002\tS002\t0.9\n"
        "S003\tS003\t0.1\n"
    )
    return p


class TestClinicalReporter:
    def test_genotype_calls(self, dosage_file: Path) -> None:
        reporter = ClinicalReporter()
        report = reporter.generate_report(dosage_file, "HLA_DRB1_0101")
        assert report.total_participants == 3
        assert report.carrier_count == 2
        calls = {g.participant_id: g.call for g in report.genotypes}
        assert calls["S001"] == "homozygous"
        assert calls["S002"] == "heterozygous"
        assert calls["S003"] == "negative"

    def test_export_csv(self, dosage_file: Path, tmp_path: Path) -> None:
        reporter = ClinicalReporter()
        report = reporter.generate_report(dosage_file, "HLA_DRB1_0101")
        out = tmp_path / "report.csv"
        ClinicalReporter.export_csv(report, out)
        content = out.read_text()
        assert "participant_id" in content
        assert "homozygous" in content
