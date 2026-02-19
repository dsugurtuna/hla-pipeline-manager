"""Clinical HLA report generator.

Integrates dosage data, quality scores, and participant metadata to produce
deliverable genotype reports for clinical teams.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class ParticipantGenotype:
    """Genotype call for a single participant."""

    participant_id: str
    allele: str
    dosage: float
    call: str  # "homozygous", "heterozygous", "negative"
    r2_score: Optional[float] = None


@dataclass
class ClinicalReport:
    """A deliverable clinical HLA report."""

    allele: str
    genotypes: List[ParticipantGenotype] = field(default_factory=list)

    @property
    def carrier_count(self) -> int:
        return sum(1 for g in self.genotypes if g.call != "negative")

    @property
    def total_participants(self) -> int:
        return len(self.genotypes)


class ClinicalReporter:
    """Generate clinical genotype reports from imputation data.

    Reads dosage files and calls genotypes based on dosage thresholds:
    - dosage >= 1.5 → homozygous
    - dosage > 0.5 → heterozygous
    - dosage <= 0.5 → negative

    Parameters
    ----------
    homozygous_threshold : float
        Dosage at or above which a homozygous call is made.
    heterozygous_threshold : float
        Dosage above which a heterozygous call is made.
    """

    def __init__(
        self,
        homozygous_threshold: float = 1.5,
        heterozygous_threshold: float = 0.5,
    ) -> None:
        self.homo_threshold = homozygous_threshold
        self.hetero_threshold = heterozygous_threshold

    def _call_genotype(self, dosage: float) -> str:
        if dosage >= self.homo_threshold:
            return "homozygous"
        elif dosage > self.hetero_threshold:
            return "heterozygous"
        return "negative"

    def generate_report(
        self,
        dosage_path: str | Path,
        allele_column: str,
        sample_col: str = "IID",
        r2_scores: Optional[Dict[str, float]] = None,
    ) -> ClinicalReport:
        """Generate a clinical report from a dosage file.

        Parameters
        ----------
        dosage_path : path
            Tab-separated dosage file.
        allele_column : str
            Column name of the target allele.
        sample_col : str
            Column containing participant IDs.
        r2_scores : dict, optional
            Marker → R² scores for quality annotation.
        """
        report = ClinicalReport(allele=allele_column)
        r2 = r2_scores.get(allele_column) if r2_scores else None

        with open(dosage_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                pid = row.get(sample_col, "").strip()
                try:
                    dosage = float(row.get(allele_column, "0"))
                except (ValueError, TypeError):
                    continue
                call = self._call_genotype(dosage)
                report.genotypes.append(
                    ParticipantGenotype(
                        participant_id=pid,
                        allele=allele_column,
                        dosage=dosage,
                        call=call,
                        r2_score=r2,
                    )
                )
        return report

    @staticmethod
    def export_csv(report: ClinicalReport, output_path: str | Path) -> None:
        """Export the clinical report to a CSV file."""
        with open(output_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["participant_id", "allele", "dosage", "call", "r2_score"])
            for g in report.genotypes:
                writer.writerow([
                    g.participant_id, g.allele, f"{g.dosage:.4f}",
                    g.call, g.r2_score or "",
                ])
