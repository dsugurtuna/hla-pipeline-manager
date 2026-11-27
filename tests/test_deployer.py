"""Tests for hla_pipeline.deployer."""

from pathlib import Path

import pytest

from hla_pipeline.deployer import MANIFEST_FILENAME, ResultDeployer, sha256
from hla_pipeline.verifier import OUTPUT_SUFFIXES


@pytest.fixture()
def source_dir(tmp_path: Path) -> Path:
    src = tmp_path / "source"
    src.mkdir()
    (src / "imputed.bed").write_bytes(b"\x6c\x1b\x01\x00")
    (src / "imputed.bim").write_text("6\tHLA_DRB1_0101\t0\t32000000\tP\tA\n")
    (src / "imputed.fam").write_text("S001 S001 0 0 1 -9\n")
    (src / "imputed.dosage").write_text("HLA_DRB1_0101\tP\tA\t1.000\n")
    (src / "imputed.bgl.r2").write_text("HLA_DRB1_0101\t0.95\n")
    (src / "imputed.bgl.log").write_text("Analysis finished.\n")
    return src


def test_default_extensions_match_verifier_outputs() -> None:
    assert ResultDeployer.DEFAULT_EXTENSIONS == OUTPUT_SUFFIXES


def test_deploy_fresh(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    report = ResultDeployer(target).deploy(source_dir)
    assert report.deployment_count == 6
    assert report.verified
    assert not report.rolled_back
    assert report.backup_dir == ""
    assert report.checksums["imputed.bed"] == sha256(target / "imputed.bed")
    assert (target / "imputed.dosage").is_file()
    assert (target / "imputed.bgl.r2").is_file()
    assert (target / "imputed.bgl.log").is_file()
    assert Path(report.manifest_path) == target / MANIFEST_FILENAME
    assert ResultDeployer.verify_manifest(target) == []


def test_deploy_with_backup_and_stale_cleanup(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    target.mkdir()
    (target / "old.bed").write_bytes(b"\xff")
    report = ResultDeployer(target, backup_root=tmp_path / "backups").deploy(source_dir)
    assert report.deployment_count == 6
    assert report.verified
    assert report.files_backed_up == ["old.bed"]
    assert (Path(report.backup_dir) / "old.bed").read_bytes() == b"\xff"
    assert (target / "imputed.bed").exists()
    assert not (target / "old.bed").exists()


def test_two_deploys_keep_separate_backups(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    deployer = ResultDeployer(target, backup_root=tmp_path / "backups")
    deployer.deploy(source_dir)
    first = deployer.deploy(source_dir)
    second = deployer.deploy(source_dir)
    assert first.backup_dir != second.backup_dir
    assert (Path(second.backup_dir) / MANIFEST_FILENAME).is_file()


def test_dry_run_writes_nothing(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    report = ResultDeployer(target).deploy(source_dir, dry_run=True)
    assert report.deployment_count == 6
    assert report.dry_run
    assert not report.verified
    assert not target.exists()


def test_corrupted_copy_is_not_verified(
    source_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def bad_copy(src: Path, dst: Path) -> None:
        Path(dst).write_bytes(b"truncated")

    monkeypatch.setattr("hla_pipeline.deployer.shutil.copy2", bad_copy)
    target = tmp_path / "production"
    report = ResultDeployer(target).deploy(source_dir)
    assert not report.verified
    assert report.rolled_back
    assert set(report.mismatched) == {f"imputed{ext}" for ext in OUTPUT_SUFFIXES}
    assert not (target / "imputed.bed").exists()
    assert not (target / MANIFEST_FILENAME).exists()


def test_corrupted_copy_rolls_back_existing_release(
    source_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "production"
    backups = tmp_path / "backups"
    deployer = ResultDeployer(target, backup_root=backups)
    first = deployer.deploy(source_dir)
    assert first.verified
    original_dosage = (target / "imputed.dosage").read_text()
    original_manifest = (target / MANIFEST_FILENAME).read_text()

    src_v2 = tmp_path / "source_v2"
    src_v2.mkdir()
    for path in source_dir.iterdir():
        (src_v2 / path.name).write_bytes(path.read_bytes())
    (src_v2 / "imputed.dosage").write_text("HLA_DRB1_0101\tP\tA\t1.950\n")

    def corrupt_atomic_copy(src: Path, dst: Path) -> None:
        dst.write_bytes(b"corrupted")

    monkeypatch.setattr("hla_pipeline.deployer._atomic_copy", corrupt_atomic_copy)
    second = deployer.deploy(src_v2)
    assert not second.verified
    assert second.rolled_back
    assert (target / "imputed.dosage").read_text() == original_dosage
    assert (target / MANIFEST_FILENAME).read_text() == original_manifest
    assert ResultDeployer.verify_manifest(target) == []


def test_deploy_rejects_incomplete_source_by_default(
    source_dir: Path, tmp_path: Path
) -> None:
    (source_dir / "imputed.bgl.log").unlink()
    target = tmp_path / "production"
    with pytest.raises(ValueError, match="failed imputation verification"):
        ResultDeployer(target).deploy(source_dir)
    assert not target.exists()


def test_deploy_rejects_missing_source_dir(tmp_path: Path) -> None:
    with pytest.raises(NotADirectoryError):
        ResultDeployer(tmp_path / "production").deploy(tmp_path / "missing")


def test_deploy_allows_unverified_when_disabled(tmp_path: Path) -> None:
    src = tmp_path / "plink_only"
    src.mkdir()
    (src / "cohort.bed").write_bytes(b"\x6c\x1b\x01\x00")
    (src / "cohort.bim").write_text("6\trs1\t0\t32000000\tA\tG\n")
    (src / "cohort.fam").write_text("S001 S001 0 0 1 -9\n")
    target = tmp_path / "production"
    report = ResultDeployer(
        target, extensions=[".bed", ".bim", ".fam"], require_verified=False
    ).deploy(src)
    assert report.verified
    assert report.deployment_count == 3
    assert ResultDeployer.verify_manifest(target) == []


def test_empty_extensions_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="extensions must not be empty"):
        ResultDeployer(tmp_path / "production", extensions=[])


def test_verify_manifest_detects_tampering(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    ResultDeployer(target).deploy(source_dir)
    (target / "imputed.dosage").write_text("tampered\n")
    (target / "imputed.bgl.log").unlink()
    assert ResultDeployer.verify_manifest(target) == [
        "imputed.bgl.log",
        "imputed.dosage",
    ]
