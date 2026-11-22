"""Tests for hla_pipeline.deployer."""

from pathlib import Path

import pytest

from hla_pipeline.deployer import ResultDeployer, sha256


@pytest.fixture()
def source_dir(tmp_path: Path) -> Path:
    src = tmp_path / "source"
    src.mkdir()
    (src / "imputed.bed").write_bytes(b"\x6c\x1b\x01\x00")
    (src / "imputed.bim").write_text("6\tHLA_DRB1_0101\t0\t32000000\tP\tA\n")
    (src / "imputed.fam").write_text("S001 S001 0 0 1 -9\n")
    return src


def test_deploy_fresh(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    report = ResultDeployer(target).deploy(source_dir)
    assert report.deployment_count == 3
    assert report.verified
    assert report.backup_dir == ""
    assert report.checksums["imputed.bed"] == sha256(target / "imputed.bed")


def test_deploy_with_backup(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    target.mkdir()
    (target / "old.bed").write_bytes(b"\xff")
    report = ResultDeployer(target, backup_root=tmp_path / "backups").deploy(source_dir)
    assert report.deployment_count == 3
    assert report.files_backed_up == ["old.bed"]
    assert (Path(report.backup_dir) / "old.bed").read_bytes() == b"\xff"
    assert (target / "imputed.bed").exists()


def test_two_deploys_keep_separate_backups(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    deployer = ResultDeployer(target, backup_root=tmp_path / "backups")
    deployer.deploy(source_dir)
    first = deployer.deploy(source_dir)
    second = deployer.deploy(source_dir)
    assert first.backup_dir != second.backup_dir


def test_dry_run_writes_nothing(source_dir: Path, tmp_path: Path) -> None:
    target = tmp_path / "production"
    report = ResultDeployer(target).deploy(source_dir, dry_run=True)
    assert report.deployment_count == 3
    assert report.dry_run
    assert not report.verified
    assert not target.exists()


def test_corrupted_copy_is_not_verified(
    source_dir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def bad_copy(src: Path, dst: Path) -> None:
        Path(dst).write_bytes(b"truncated")

    monkeypatch.setattr("hla_pipeline.deployer.shutil.copy2", bad_copy)
    report = ResultDeployer(tmp_path / "production").deploy(source_dir)
    assert not report.verified
    assert set(report.mismatched) == {"imputed.bed", "imputed.bim", "imputed.fam"}
