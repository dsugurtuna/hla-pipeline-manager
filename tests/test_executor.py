"""Tests for hla_pipeline.executor."""

from pathlib import Path

import pytest

from hla_pipeline.executor import BatchExecutor, PipelineConfig


@pytest.fixture()
def bfile(tmp_path: Path) -> Path:
    prefix = tmp_path / "input" / "b01"
    prefix.parent.mkdir()
    lines = [f"S{i:03d} S{i:03d} 0 0 1 -9" for i in range(1, 11)]
    Path(str(prefix) + ".fam").write_text("\n".join(lines) + "\n")
    return prefix


def _executor(size: int = 5) -> BatchExecutor:
    return BatchExecutor(PipelineConfig(reference_panel="REF", sub_batch_size=size))


def test_count_samples(bfile: Path) -> None:
    assert BatchExecutor().count_samples(str(bfile) + ".fam") == 10


def test_split_fam_writes_keep_files(bfile: Path, tmp_path: Path) -> None:
    keeps = _executor(3).split_fam(str(bfile) + ".fam", tmp_path / "keep")
    assert len(keeps) == 4  # 3 + 3 + 3 + 1
    assert all(k.suffix == ".keep" for k in keeps)
    assert len(keeps[0].read_text().splitlines()) == 3
    assert len(keeps[-1].read_text().splitlines()) == 1


def test_plan_extracts_mhc_then_sub_batches(bfile: Path, tmp_path: Path) -> None:
    plan = _executor().plan_batch("b01", bfile, tmp_path / "work")
    extract = plan.shared_steps[0].argv
    assert extract[:3] == ["plink", "--bfile", str(bfile)]
    assert extract[extract.index("--chr") + 1] == "6"
    assert extract[extract.index("--from-bp") + 1] == "26000000"
    assert extract[extract.index("--to-bp") + 1] == "34000000"
    assert list(plan.sub_batch_steps) == ["sub_batch_001", "sub_batch_002"]
    # PLINK will not create --out directories itself.
    assert (tmp_path / "work" / "sub_batches").is_dir()


def test_rename_step_uses_update_name_columns(bfile: Path, tmp_path: Path) -> None:
    rename = BatchExecutor.write_rename_file(
        {"AX-1": "rs1", "AX-2": "rs2"}, tmp_path / "rename.txt"
    )
    assert rename.read_text() == "AX-1\trs1\nAX-2\trs2\n"
    plan = _executor().plan_batch("b01", bfile, tmp_path / "work", rename)
    argv = plan.shared_steps[1].argv
    i = argv.index("--update-name")
    assert argv[i + 1 : i + 4] == [str(rename), "2", "1"]
    # Sub-batches are cut from the renamed fileset.
    make_sub = plan.sub_batch_steps["sub_batch_001"][0].argv
    assert make_sub[make_sub.index("--bfile") + 1].endswith("b01_mhc_renamed")


def test_snp2hla_argument_order() -> None:
    cfg = PipelineConfig(
        reference_panel="HM_CEU_REF", java_max_memory_mb=8000, marker_window_size=500
    )
    argv = BatchExecutor(cfg).snp2hla_command("in", "out")
    # SNP2HLA.csh DATA REFERENCE OUTPUT plink [java_max_memory_mb] [window]
    assert argv == ["SNP2HLA.csh", "in", "HM_CEU_REF", "out", "plink", "8000", "500"]


def test_plan_requires_reference(bfile: Path, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="reference_panel"):
        BatchExecutor().plan_batch("b01", bfile, tmp_path / "work")


def test_dry_run_completes_nothing(bfile: Path, tmp_path: Path) -> None:
    result = _executor().execute_batch("b01", bfile, tmp_path / "work")
    assert result.dry_run
    assert result.total_samples == 10
    assert result.sub_batches == 2
    assert result.completed == 0
    assert result.success_rate == 0.0
    assert len(result.commands) == 1 + 2 * 2


def test_run_all_succeed(bfile: Path, tmp_path: Path) -> None:
    calls: list[list[str]] = []

    def runner(argv: list[str]) -> int:
        calls.append(argv)
        return 0

    result = _executor().execute_batch(
        "b01", bfile, tmp_path / "work", dry_run=False, runner=runner
    )
    assert result.completed == 2
    assert result.success_rate == 1.0
    assert calls == result.commands


def test_run_reports_failed_sub_batch(bfile: Path, tmp_path: Path) -> None:
    def runner(argv: list[str]) -> int:
        is_snp2hla = argv[0] == "SNP2HLA.csh"
        return 1 if is_snp2hla and "sub_batch_002" in argv[1] else 0

    result = _executor().execute_batch(
        "b01", bfile, tmp_path / "work", dry_run=False, runner=runner
    )
    assert result.completed == 1
    assert result.failed_sub_batches == ["sub_batch_002"]
    assert result.success_rate == 0.5


def test_shared_step_failure_stops_everything(bfile: Path, tmp_path: Path) -> None:
    calls: list[list[str]] = []

    def runner(argv: list[str]) -> int:
        calls.append(argv)
        return 1

    result = _executor().execute_batch(
        "b01", bfile, tmp_path / "work", dry_run=False, runner=runner
    )
    assert len(calls) == 1
    assert result.completed == 0
    assert result.failed == 2


def test_build_rename_map_skips_unmapped(tmp_path: Path) -> None:
    annot = tmp_path / "annot.csv"
    annot.write_text("probesetid,rsid\nAX-1,rs10\nAX-2,---\nAX-3,rs30\n")
    assert BatchExecutor().build_rename_map(annot) == {"AX-1": "rs10", "AX-3": "rs30"}
