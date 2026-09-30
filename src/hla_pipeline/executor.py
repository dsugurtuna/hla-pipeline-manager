"""Batch imputation planner and runner.

Builds the PLINK and SNP2HLA commands for one genotyping batch:

1. extract the extended MHC region on chromosome 6 (default 26-34 Mb),
2. optionally rename Affymetrix probe set IDs (AX-...) to rsIDs,
3. split the samples into sub-batches of a fixed size,
4. make one PLINK fileset per sub-batch, and
5. run SNP2HLA on each sub-batch.

Planning never touches the tools. Running is explicit (``dry_run=False``) and
goes through an injectable runner, so the orchestration logic is testable
without PLINK, Java or a reference panel installed.
"""

from __future__ import annotations

import csv
import subprocess
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

Runner = Callable[[list[str]], int]
"""Runs one command (argv list) and returns its exit code."""


def subprocess_runner(argv: list[str]) -> int:
    """Default runner: execute the command and return its exit code."""
    return subprocess.run(argv, check=False).returncode


@dataclass
class PipelineConfig:
    """Tool paths and parameters for the imputation pipeline."""

    plink_path: str = "plink"
    snp2hla_path: str = "SNP2HLA.csh"
    reference_panel: str = ""
    chromosome: int = 6
    mhc_start_bp: int = 26_000_000
    mhc_end_bp: int = 34_000_000
    sub_batch_size: int = 500
    # SNP2HLA.csh optional arguments 5 and 6. The script's own defaults are
    # 2000 MB of Java heap and a Beagle window of 1000 markers.
    java_max_memory_mb: int = 2000
    marker_window_size: int = 1000


@dataclass
class PlannedCommand:
    """One command in an execution plan."""

    step: str
    argv: list[str]
    sub_batch: str | None = None


@dataclass
class ExecutionPlan:
    """All commands needed to impute one batch, in order."""

    batch_id: str
    total_samples: int
    shared_steps: list[PlannedCommand] = field(default_factory=list)
    sub_batch_steps: dict[str, list[PlannedCommand]] = field(default_factory=dict)

    @property
    def commands(self) -> list[list[str]]:
        """Every argv in execution order."""
        out = [c.argv for c in self.shared_steps]
        for steps in self.sub_batch_steps.values():
            out.extend(c.argv for c in steps)
        return out


@dataclass
class ExecutionResult:
    """Outcome of planning or running one batch."""

    batch_id: str
    total_samples: int = 0
    sub_batches: int = 0
    completed: int = 0
    failed_sub_batches: list[str] = field(default_factory=list)
    commands: list[list[str]] = field(default_factory=list)
    dry_run: bool = True

    @property
    def failed(self) -> int:
        return len(self.failed_sub_batches)

    @property
    def success_rate(self) -> float:
        """Share of sub-batches whose commands all exited with status 0.

        Always 0.0 for a dry run, because nothing was executed.
        """
        if self.sub_batches == 0:
            return 0.0
        return self.completed / self.sub_batches


class BatchExecutor:
    """Plan and (optionally) run HLA imputation for one batch."""

    def __init__(self, config: PipelineConfig | None = None) -> None:
        self.config = config or PipelineConfig()

    def count_samples(self, fam_path: str | Path) -> int:
        """Count samples (non-blank lines) in a .fam file."""
        with open(fam_path) as fh:
            return sum(1 for line in fh if line.strip())

    def split_fam(
        self,
        fam_path: str | Path,
        output_dir: str | Path,
        batch_size: int | None = None,
    ) -> list[Path]:
        """Split a .fam file into PLINK ``--keep`` files of ``batch_size`` samples.

        Keep files use the ``.keep`` extension so that PLINK writing
        ``sub_batch_001.fam`` can never overwrite its own input.
        """
        batch_size = batch_size or self.config.sub_batch_size
        if batch_size < 1:
            raise ValueError("batch_size must be at least 1")
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        with open(fam_path) as fh:
            samples = [line.strip() for line in fh if line.strip()]

        keep_files: list[Path] = []
        for i in range(0, len(samples), batch_size):
            chunk = samples[i : i + batch_size]
            path = output_dir / f"sub_batch_{i // batch_size + 1:03d}.keep"
            path.write_text("\n".join(chunk) + "\n")
            keep_files.append(path)
        return keep_files

    def build_rename_map(
        self,
        mapping_path: str | Path,
        ax_col: str = "probesetid",
        rs_col: str = "rsid",
    ) -> dict[str, str]:
        """Build an AX -> rs ID map from an array annotation CSV.

        Rows without an ``rs`` identifier are skipped, so unmapped probes keep
        their AX- names.
        """
        rename: dict[str, str] = {}
        with open(mapping_path, newline="") as fh:
            for row in csv.DictReader(fh):
                ax = (row.get(ax_col) or "").strip()
                rs = (row.get(rs_col) or "").strip()
                if ax and rs.startswith("rs"):
                    rename[ax] = rs
        return rename

    @staticmethod
    def write_rename_file(mapping: Mapping[str, str], path: str | Path) -> Path:
        """Write ``old<TAB>new`` lines for ``plink --update-name <file> 2 1``."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            for old, new in sorted(mapping.items()):
                fh.write(f"{old}\t{new}\n")
        return path

    def snp2hla_command(self, input_prefix: str, output_prefix: str) -> list[str]:
        """SNP2HLA.csh DATA REFERENCE OUTPUT plink [java_max_memory_mb] [window]."""
        cfg = self.config
        return [
            cfg.snp2hla_path,
            input_prefix,
            cfg.reference_panel,
            output_prefix,
            cfg.plink_path,
            str(cfg.java_max_memory_mb),
            str(cfg.marker_window_size),
        ]

    def plan_batch(
        self,
        batch_id: str,
        bfile: str | Path,
        work_dir: str | Path,
        rename_file: str | Path | None = None,
    ) -> ExecutionPlan:
        """Build the full command plan for one batch.

        Parameters
        ----------
        batch_id : str
            Label used in output file names.
        bfile : path
            PLINK binary fileset prefix (``<bfile>.bed/.bim/.fam``).
        work_dir : path
            Directory for intermediate and output files.
        rename_file : path, optional
            Two-column ``old new`` file (see :meth:`write_rename_file`).
        """
        cfg = self.config
        if not cfg.reference_panel:
            raise ValueError("PipelineConfig.reference_panel must be set")
        work = Path(work_dir)
        bfile = str(bfile)
        fam = Path(bfile + ".fam")

        plan = ExecutionPlan(batch_id=batch_id, total_samples=self.count_samples(fam))

        mhc = str(work / f"{batch_id}_mhc")
        plan.shared_steps.append(
            PlannedCommand(
                step="extract_mhc",
                argv=[
                    cfg.plink_path,
                    "--bfile",
                    bfile,
                    "--chr",
                    str(cfg.chromosome),
                    "--from-bp",
                    str(cfg.mhc_start_bp),
                    "--to-bp",
                    str(cfg.mhc_end_bp),
                    "--make-bed",
                    "--out",
                    mhc,
                ],
            )
        )
        source = mhc
        if rename_file is not None:
            renamed = str(work / f"{batch_id}_mhc_renamed")
            plan.shared_steps.append(
                PlannedCommand(
                    step="rename_ids",
                    argv=[
                        cfg.plink_path,
                        "--bfile",
                        mhc,
                        "--update-name",
                        str(rename_file),
                        "2",
                        "1",
                        "--make-bed",
                        "--out",
                        renamed,
                    ],
                )
            )
            source = renamed

        # Region filtering removes variants, never samples, so the input .fam
        # defines the sub-batches before any tool has run.
        keep_files = self.split_fam(fam, work / "keep")
        sub_dir = work / "sub_batches"
        # PLINK does not create missing output directories.
        sub_dir.mkdir(parents=True, exist_ok=True)
        for keep in keep_files:
            name = keep.stem
            sub_prefix = str(sub_dir / name)
            plan.sub_batch_steps[name] = [
                PlannedCommand(
                    step="make_sub_batch",
                    sub_batch=name,
                    argv=[
                        cfg.plink_path,
                        "--bfile",
                        source,
                        "--keep",
                        str(keep),
                        "--make-bed",
                        "--out",
                        sub_prefix,
                    ],
                ),
                PlannedCommand(
                    step="snp2hla",
                    sub_batch=name,
                    argv=self.snp2hla_command(sub_prefix, sub_prefix + "_imputed"),
                ),
            ]
        return plan

    def run_plan(
        self, plan: ExecutionPlan, runner: Runner | None = None
    ) -> ExecutionResult:
        """Run a plan, stopping a sub-batch at its first failing command.

        If a shared step fails, no sub-batch is attempted and all are reported
        as failed.
        """
        run = runner or subprocess_runner
        result = ExecutionResult(
            batch_id=plan.batch_id,
            total_samples=plan.total_samples,
            sub_batches=len(plan.sub_batch_steps),
            commands=plan.commands,
            dry_run=False,
        )
        for cmd in plan.shared_steps:
            if run(cmd.argv) != 0:
                result.failed_sub_batches = list(plan.sub_batch_steps)
                return result
        for name, steps in plan.sub_batch_steps.items():
            if all(run(cmd.argv) == 0 for cmd in steps):
                result.completed += 1
            else:
                result.failed_sub_batches.append(name)
        return result

    def execute_batch(
        self,
        batch_id: str,
        bfile: str | Path,
        work_dir: str | Path,
        rename_file: str | Path | None = None,
        dry_run: bool = True,
        runner: Runner | None = None,
    ) -> ExecutionResult:
        """Plan a batch and, unless ``dry_run`` is True, run it.

        A dry run returns the planned commands with ``completed == 0``: it
        reports what would happen, not what did.
        """
        plan = self.plan_batch(batch_id, bfile, work_dir, rename_file)
        if dry_run:
            return ExecutionResult(
                batch_id=batch_id,
                total_samples=plan.total_samples,
                sub_batches=len(plan.sub_batch_steps),
                commands=plan.commands,
                dry_run=True,
            )
        return self.run_plan(plan, runner)
