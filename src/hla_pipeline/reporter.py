"""HLA allele carrier report.

Turns per-sample allele dosages into hard calls for one HLA allele and writes
a CSV. Calls use fixed dosage thresholds on the expected count of the allele:

- dosage >= 1.5 -> homozygous
- dosage >  0.5 -> heterozygous
- otherwise     -> negative

These are the thresholds the original scripts used. They are a convenience
for listing carriers, not a clinical genotype: imputed dosages carry
uncertainty, which is why the marker's Beagle r2 is reported alongside.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from .dosage import read_plink_raw, read_snp2hla_dosage


@dataclass
class ParticipantGenotype:
    """Hard call for one participant."""

    participant_id: str
    allele: str
    dosage: float
    call: str  # "homozygous", "heterozygous" or "negative"
    r2_score: float | None = None


@dataclass
class ClinicalReport:
    """Hard calls for one HLA allele."""

    allele: str
    genotypes: list[ParticipantGenotype] = field(default_factory=list)

    @property
    def carrier_count(self) -> int:
        return sum(1 for g in self.genotypes if g.call != "negative")

    @property
    def total_participants(self) -> int:
        return len(self.genotypes)


class ClinicalReporter:
    """Make carrier reports from imputed HLA dosages."""

    def __init__(
        self,
        homozygous_threshold: float = 1.5,
        heterozygous_threshold: float = 0.5,
    ) -> None:
        if not 0 <= heterozygous_threshold < homozygous_threshold <= 2:
            raise ValueError("need 0 <= heterozygous < homozygous <= 2")
        self.homo_threshold = homozygous_threshold
        self.hetero_threshold = heterozygous_threshold

    def _call_genotype(self, dosage: float) -> str:
        if dosage >= self.homo_threshold:
            return "homozygous"
        if dosage > self.hetero_threshold:
            return "heterozygous"
        return "negative"

    def report_from_dosages(
        self,
        dosages: dict[str, float],
        allele: str,
        r2_scores: dict[str, float] | None = None,
    ) -> ClinicalReport:
        """Build a report from ``{participant_id: dosage}``."""
        r2 = r2_scores.get(allele) if r2_scores else None
        report = ClinicalReport(allele=allele)
        for pid, dosage in dosages.items():
            report.genotypes.append(
                ParticipantGenotype(
                    participant_id=pid,
                    allele=allele,
                    dosage=dosage,
                    call=self._call_genotype(dosage),
                    r2_score=r2,
                )
            )
        return report

    def generate_report(
        self,
        dosage_path: str | Path,
        allele_column: str,
        r2_scores: dict[str, float] | None = None,
    ) -> ClinicalReport:
        """Build a report from a PLINK ``--recode A`` (.raw) file.

        ``allele_column`` may be the bare marker (``HLA_DRB1_0401``) or PLINK's
        suffixed column name (``HLA_DRB1_0401_P``).
        """
        marker = allele_column
        if marker.endswith(("_P", "_A")) and marker.count("_") > 2:
            marker = marker[:-2]
        return self.report_from_dosages(
            read_plink_raw(dosage_path, marker), marker, r2_scores
        )

    def generate_report_snp2hla(
        self,
        dosage_path: str | Path,
        fam_path: str | Path,
        allele: str,
        r2_scores: dict[str, float] | None = None,
    ) -> ClinicalReport:
        """Build a report from a SNP2HLA ``.dosage`` file and its ``.fam``."""
        return self.report_from_dosages(
            read_snp2hla_dosage(dosage_path, fam_path, allele), allele, r2_scores
        )

    @staticmethod
    def export_csv(report: ClinicalReport, output_path: str | Path) -> None:
        """Write the report as CSV."""
        with open(output_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["participant_id", "allele", "dosage", "call", "r2_score"])
            for g in report.genotypes:
                writer.writerow(
                    [
                        g.participant_id,
                        g.allele,
                        f"{g.dosage:.4f}",
                        g.call,
                        "" if g.r2_score is None else f"{g.r2_score:.4f}",
                    ]
                )
