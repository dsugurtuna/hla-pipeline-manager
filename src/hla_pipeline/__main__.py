"""Command-line interface: ``python -m hla_pipeline <command>``.

Commands only plan, read and report. Deployment (which overwrites release
files) is deliberately left to the Python API, where it has to be called on
purpose.
"""

from __future__ import annotations

import argparse
import math
import shlex
import sys
from collections.abc import Sequence

from .dosage import read_r2
from .executor import BatchExecutor, PipelineConfig
from .reporter import ClinicalReporter
from .verifier import ImputationVerifier


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="hla_pipeline")
    sub = parser.add_subparsers(dest="command", required=True)

    plan = sub.add_parser("plan", help="print the PLINK/SNP2HLA commands for a batch")
    plan.add_argument("--batch-id", required=True)
    plan.add_argument("--bfile", required=True, help="PLINK fileset prefix")
    plan.add_argument("--work-dir", required=True)
    plan.add_argument("--reference", required=True, help="SNP2HLA reference prefix")
    plan.add_argument("--rename-file", help="two-column old/new variant ID file")
    plan.add_argument("--sub-batch-size", type=int, default=500)

    verify = sub.add_parser("verify", help="check SNP2HLA outputs in a directory")
    verify.add_argument("batch_dir")
    verify.add_argument("--expected", type=int, help="expected number of sub-batches")
    verify.add_argument("--samples", type=int, help="samples in batch (with --size)")
    verify.add_argument("--size", type=int, help="sub-batch size (with --samples)")

    report = sub.add_parser("report", help="carrier report for one HLA allele")
    report.add_argument("dosage", help="PLINK .raw, or SNP2HLA .dosage with --fam")
    report.add_argument("--allele", required=True, help="e.g. HLA_DRB1_0401")
    report.add_argument("--fam", help=".fam from the same SNP2HLA run")
    report.add_argument("--r2", help="Beagle .bgl.r2 file for quality annotation")
    report.add_argument("--out", help="write CSV here instead of printing")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.command == "plan":
        config = PipelineConfig(
            reference_panel=args.reference, sub_batch_size=args.sub_batch_size
        )
        result = BatchExecutor(config).execute_batch(
            args.batch_id, args.bfile, args.work_dir, args.rename_file, dry_run=True
        )
        print(
            f"# {result.total_samples} samples -> {result.sub_batches} sub-batches"
            " (dry run: nothing executed)"
        )
        for cmd in result.commands:
            print(shlex.join(cmd))
        return 0

    if args.command == "verify":
        expected = args.expected
        if expected is None and args.samples and args.size:
            expected = math.ceil(args.samples / args.size)
        verifier = ImputationVerifier()
        vreport = verifier.verify_batch(args.batch_dir, expected_sub_batches=expected)
        print(verifier.format_report(vreport))
        return 0 if vreport.all_ok else 1

    reporter = ClinicalReporter()
    r2 = read_r2(args.r2) if args.r2 else None
    if args.fam:
        creport = reporter.generate_report_snp2hla(
            args.dosage, args.fam, args.allele, r2
        )
    else:
        creport = reporter.generate_report(args.dosage, args.allele, r2)
    if args.out:
        reporter.export_csv(creport, args.out)
    else:
        for g in creport.genotypes:
            r2_text = "" if g.r2_score is None else f"  r2={g.r2_score:.2f}"
            print(f"{g.participant_id}\t{g.dosage:.3f}\t{g.call}{r2_text}")
    print(
        f"# {creport.allele}: {creport.carrier_count} carriers of "
        f"{creport.total_participants} samples",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
