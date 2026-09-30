"""Imputation output verification.

Checks each SNP2HLA sub-batch output (``<prefix>.bed/.bim/.fam/.dosage/
.bgl.r2/.bgl.log``) for presence, non-empty content, HLA allele markers in
the .bim and .bgl.r2 files, and a Beagle completion line in the log. These
are the same checks the original shell script made, in testable form.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

OUTPUT_SUFFIXES = (".bed", ".bim", ".fam", ".dosage", ".bgl.r2", ".bgl.log")

# Heuristic end-of-run markers. The original script accepted "finished" or
# "End time" near the end of the log; keep them configurable.
DEFAULT_COMPLETION_MARKERS = ("finished", "completed", "end time")


def _with_suffix(prefix: Path, suffix: str) -> Path:
    # Path.with_suffix would replace the text after the last dot, which breaks
    # prefixes such as "b28.sub_001". Append instead.
    return Path(str(prefix) + suffix)


def _non_empty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


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
    hla_r2_count: int = 0
    beagle_completed: bool = False

    @property
    def problems(self) -> list[str]:
        """Human-readable reasons this sub-batch is incomplete (empty if none)."""
        checks = [
            (self.has_bed, "missing or empty .bed"),
            (self.has_bim, "missing or empty .bim"),
            (self.has_fam, "missing or empty .fam"),
            (self.has_dosage, "missing or empty .dosage"),
            (self.has_r2, "missing or empty .bgl.r2"),
            (self.has_log, "missing or empty .bgl.log"),
            (self.hla_marker_count > 0, "no HLA markers in .bim"),
            (self.hla_r2_count > 0, "no HLA markers in .bgl.r2"),
            (self.beagle_completed, "no completion line in Beagle log"),
        ]
        return [message for ok, message in checks if not ok]

    @property
    def is_complete(self) -> bool:
        return not self.problems


@dataclass
class VerificationReport:
    """Verification report for one or more batches."""

    batch_statuses: dict[str, list[SubBatchStatus]] = field(default_factory=dict)
    expected_sub_batches: dict[str, int] = field(default_factory=dict)

    @property
    def total_sub_batches(self) -> int:
        return sum(len(v) for v in self.batch_statuses.values())

    @property
    def complete_sub_batches(self) -> int:
        return sum(
            1 for subs in self.batch_statuses.values() for s in subs if s.is_complete
        )

    @property
    def completeness_rate(self) -> float:
        if self.total_sub_batches == 0:
            return 0.0
        return self.complete_sub_batches / self.total_sub_batches

    def missing_sub_batches(self, batch_id: str) -> int:
        """Expected minus found sub-batches (0 when no expectation was given)."""
        expected = self.expected_sub_batches.get(batch_id)
        if expected is None:
            return 0
        return max(0, expected - len(self.batch_statuses.get(batch_id, [])))

    @property
    def all_ok(self) -> bool:
        """True when every sub-batch is complete and none are missing."""
        return all(
            s.is_complete for subs in self.batch_statuses.values() for s in subs
        ) and not any(self.missing_sub_batches(b) for b in self.batch_statuses)


class ImputationVerifier:
    """Verify completeness of SNP2HLA outputs."""

    def __init__(
        self, completion_markers: Sequence[str] = DEFAULT_COMPLETION_MARKERS
    ) -> None:
        self.completion_markers = tuple(m.lower() for m in completion_markers)

    @staticmethod
    def _count_hla_markers(bim_path: Path) -> int:
        """Count .bim rows whose variant ID (column 2) starts with ``HLA_``."""
        count = 0
        with open(bim_path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 2 and parts[1].startswith("HLA_"):
                    count += 1
        return count

    @staticmethod
    def _count_hla_r2(r2_path: Path) -> int:
        """Count .bgl.r2 rows for HLA allele markers."""
        with open(r2_path) as fh:
            return sum(1 for line in fh if line.startswith("HLA_"))

    def _check_beagle_log(self, log_path: Path, tail_lines: int = 20) -> bool:
        """True if a completion marker appears in the last lines of the log."""
        with open(log_path, errors="replace") as fh:
            tail = fh.readlines()[-tail_lines:]
        text = "".join(tail).lower()
        return any(marker in text for marker in self.completion_markers)

    def verify_sub_batch(self, prefix: str | Path) -> SubBatchStatus:
        """Verify one sub-batch given its output prefix (no extension)."""
        prefix = Path(prefix)
        status = SubBatchStatus(name=prefix.name)
        bim = _with_suffix(prefix, ".bim")
        r2 = _with_suffix(prefix, ".bgl.r2")
        log = _with_suffix(prefix, ".bgl.log")
        status.has_bed = _non_empty(_with_suffix(prefix, ".bed"))
        status.has_bim = _non_empty(bim)
        status.has_fam = _non_empty(_with_suffix(prefix, ".fam"))
        status.has_dosage = _non_empty(_with_suffix(prefix, ".dosage"))
        status.has_r2 = _non_empty(r2)
        status.has_log = _non_empty(log)
        if status.has_bim:
            status.hla_marker_count = self._count_hla_markers(bim)
        if status.has_r2:
            status.hla_r2_count = self._count_hla_r2(r2)
        if status.has_log:
            status.beagle_completed = self._check_beagle_log(log)
        return status

    def verify_batch(
        self,
        batch_dir: str | Path,
        batch_id: str = "",
        expected_sub_batches: int | None = None,
    ) -> VerificationReport:
        """Verify every sub-batch in a directory.

        Sub-batches are discovered from ``*.dosage`` files, which only SNP2HLA
        writes, so the pre-imputation PLINK filesets in the same directory are
        not mistaken for outputs. Pass ``expected_sub_batches`` (for example
        ``ceil(samples / sub_batch_size)``) to detect sub-batches that never
        produced output at all.
        """
        batch_dir = Path(batch_dir)
        bid = batch_id or batch_dir.name
        report = VerificationReport()
        report.batch_statuses[bid] = [
            self.verify_sub_batch(d.with_suffix(""))
            for d in sorted(batch_dir.glob("*.dosage"))
        ]
        if expected_sub_batches is not None:
            report.expected_sub_batches[bid] = expected_sub_batches
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
            missing = report.missing_sub_batches(batch)
            if missing:
                lines.append(f"  MISSING: {missing} expected sub-batch(es) not found")
            for s in subs:
                state = "OK" if s.is_complete else "INCOMPLETE"
                lines.append(
                    f"  {s.name}: {state}  (HLA markers in .bim: "
                    f"{s.hla_marker_count}, in .bgl.r2: {s.hla_r2_count})"
                )
                lines.extend(f"    - {problem}" for problem in s.problems)
        return "\n".join(lines)
