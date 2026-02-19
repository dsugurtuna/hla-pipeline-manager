"""Imputation verification module.

Checks the completeness of HLA imputation output files across batches
and sub-batches: .bed/.bim/.fam/.dosage/.bgl.r2/.bgl.log presence,
HLA marker counts, and Beagle completion messages.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Set


@dataclass
class SubBatchStatus:
    """Status of a single sub-batch."""

    name: str
    has_bed: bool = False
    has_bim: bool = False
    has_fam: bool = False
    has_dosage: bool = False
    has_r2: bool = False
    has_log: bool = False
    hla_marker_count: int = 0
    beagle_completed: bool = False

    @property
    def is_complete(self) -> bool:
        return all([
            self.has_bed, self.has_bim, self.has_fam,
            self.has_dosage, self.has_r2, self.has_log,
            self.hla_marker_count > 0, self.beagle_completed,
        ])


@dataclass
class VerificationReport:
    """Verification report for one or more batches."""

    batch_statuses: Dict[str, List[SubBatchStatus]] = field(default_factory=dict)

    @property
    def total_sub_batches(self) -> int:
        return sum(len(v) for v in self.batch_statuses.values())

    @property
    def complete_sub_batches(self) -> int:
        return sum(
            1 for subs in self.batch_statuses.values()
            for s in subs if s.is_complete
        )

    @property
    def completeness_rate(self) -> float:
        if self.total_sub_batches == 0:
            return 0.0
        return self.complete_sub_batches / self.total_sub_batches


class ImputationVerifier:
    """Verify completeness of HLA imputation outputs.

    Parameters
    ----------
    expected_extensions : list of str
        File extensions to check for each sub-batch.
    """

    EXPECTED_EXTENSIONS = [".bed", ".bim", ".fam", ".dosage", ".bgl.r2", ".bgl.log"]

    def _count_hla_markers(self, bim_path: Path) -> int:
        """Count markers containing 'HLA_' in a .bim file."""
        count = 0
        if bim_path.exists():
            with open(bim_path) as fh:
                for line in fh:
                    if "HLA_" in line:
                        count += 1
        return count

    def _check_beagle_log(self, log_path: Path) -> bool:
        """Check if a Beagle log contains a completion marker."""
        if not log_path.exists():
            return False
        with open(log_path) as fh:
            for line in fh:
                if "finished" in line.lower() or "completed" in line.lower():
                    return True
        return False

    def verify_sub_batch(self, prefix: str | Path) -> SubBatchStatus:
        """Verify a single sub-batch by file prefix.

        Parameters
        ----------
        prefix : path
            Common prefix (e.g. ``/data/batch1/sub_001``). The method
            checks for ``prefix.bed``, ``prefix.bim``, etc.
        """
        prefix = Path(prefix)
        status = SubBatchStatus(name=prefix.name)
        status.has_bed = prefix.with_suffix(".bed").exists()
        status.has_bim = prefix.with_suffix(".bim").exists()
        status.has_fam = prefix.with_suffix(".fam").exists()
        status.has_dosage = prefix.with_suffix(".dosage").exists()
        status.has_r2 = Path(str(prefix) + ".bgl.r2").exists()
        status.has_log = Path(str(prefix) + ".bgl.log").exists()
        if status.has_bim:
            status.hla_marker_count = self._count_hla_markers(
                prefix.with_suffix(".bim")
            )
        if status.has_log:
            status.beagle_completed = self._check_beagle_log(
                Path(str(prefix) + ".bgl.log")
            )
        return status

    def verify_batch(
        self,
        batch_dir: str | Path,
        batch_id: str = "",
    ) -> VerificationReport:
        """Verify all sub-batches in a directory.

        Discovers sub-batches by looking for .fam files.
        """
        batch_dir = Path(batch_dir)
        bid = batch_id or batch_dir.name
        report = VerificationReport()
        sub_statuses: List[SubBatchStatus] = []

        for fam in sorted(batch_dir.glob("*.fam")):
            prefix = fam.with_suffix("")
            sub_statuses.append(self.verify_sub_batch(prefix))

        report.batch_statuses[bid] = sub_statuses
        return report

    @staticmethod
    def format_report(report: VerificationReport) -> str:
        """Format a human-readable verification report."""
        lines = [
            "HLA Imputation Verification Report",
            "=" * 45,
            f"Total sub-batches:    {report.total_sub_batches}",
            f"Complete:             {report.complete_sub_batches}",
            f"Completeness rate:    {report.completeness_rate:.1%}",
            "",
        ]
        for batch, subs in sorted(report.batch_statuses.items()):
            lines.append(f"Batch: {batch}")
            for s in subs:
                status = "OK" if s.is_complete else "INCOMPLETE"
                lines.append(f"  {s.name}: {status}  (HLA markers: {s.hla_marker_count})")
        return "\n".join(lines)
