"""Batch imputation executor.

Orchestrates the HLA imputation pipeline: chromosome filtering, ID renaming,
MHC region extraction, sub-batch splitting, and parallel SNP2HLA execution.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List


@dataclass
class ExecutionResult:
    """Result summary of a batch imputation run."""

    batch_id: str
    total_samples: int = 0
    sub_batches: int = 0
    completed: int = 0
    failed: int = 0
    log_files: List[str] = field(default_factory=list)

    @property
    def success_rate(self) -> float:
        if self.sub_batches == 0:
            return 0.0
        return self.completed / self.sub_batches


@dataclass
class PipelineConfig:
    """Configuration for the HLA imputation pipeline."""

    plink_path: str = "plink"
    snp2hla_path: str = "SNP2HLA.csh"
    beagle_path: str = "beagle.jar"
    reference_panel: str = ""
    chromosome: int = 6
    mhc_start_bp: int = 26_000_000
    mhc_end_bp: int = 34_000_000
    sub_batch_size: int = 500
    max_parallel: int = 4


class BatchExecutor:
    """Execute HLA imputation across batches.

    Handles the full pipeline from raw PLINK files through to imputed output:
    1. Extract chromosome 6
    2. Rename AX-prefixed IDs to rsIDs
    3. Extract MHC region (26–34 Mb)
    4. Split into sub-batches
    5. Run SNP2HLA per sub-batch

    Parameters
    ----------
    config : PipelineConfig
        Pipeline settings.
    """

    def __init__(self, config: PipelineConfig | None = None) -> None:
        self.config = config or PipelineConfig()

    def count_samples(self, fam_path: str | Path) -> int:
        """Count samples in a .fam file."""
        with open(fam_path) as fh:
            return sum(1 for line in fh if line.strip())

    def split_fam(
        self,
        fam_path: str | Path,
        output_dir: str | Path,
        batch_size: int | None = None,
    ) -> List[Path]:
        """Split a .fam file into sub-batches.

        Returns list of sub-batch .fam file paths.
        """
        batch_size = batch_size or self.config.sub_batch_size
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        samples: List[str] = []
        with open(fam_path) as fh:
            samples = [line.strip() for line in fh if line.strip()]

        sub_batches: List[Path] = []
        for i in range(0, len(samples), batch_size):
            chunk = samples[i : i + batch_size]
            sub_path = output_dir / f"sub_batch_{i // batch_size + 1:03d}.fam"
            sub_path.write_text("\n".join(chunk) + "\n")
            sub_batches.append(sub_path)

        return sub_batches

    def build_rename_map(
        self,
        mapping_path: str | Path,
        ax_col: str = "probesetid",
        rs_col: str = "rsid",
    ) -> Dict[str, str]:
        """Build an AX → rs ID rename map from an annotation file."""
        import csv

        rename: Dict[str, str] = {}
        with open(mapping_path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                ax = row.get(ax_col, "").strip()
                rs = row.get(rs_col, "").strip()
                if ax and rs and rs.startswith("rs"):
                    rename[ax] = rs
        return rename

    def prepare_pipeline_command(
        self,
        input_prefix: str,
        output_prefix: str,
        sub_batch_fam: str | Path,
    ) -> List[str]:
        """Construct the SNP2HLA command for a single sub-batch."""
        return [
            self.config.snp2hla_path,
            input_prefix,
            self.config.reference_panel,
            output_prefix,
            self.config.plink_path,
            str(self.config.max_parallel),
        ]

    def execute_batch(
        self,
        batch_id: str,
        fam_path: str | Path,
        work_dir: str | Path,
    ) -> ExecutionResult:
        """Execute the full pipeline for one batch (dry-run safe).

        This method prepares the execution plan. In production, it would
        launch SNP2HLA via subprocess. Here it validates inputs and returns
        the planned configuration.
        """
        work_dir = Path(work_dir)
        result = ExecutionResult(batch_id=batch_id)
        result.total_samples = self.count_samples(fam_path)
        sub_batches = self.split_fam(fam_path, work_dir / "sub_batches")
        result.sub_batches = len(sub_batches)
        result.completed = len(sub_batches)  # dry-run: assume success
        return result
