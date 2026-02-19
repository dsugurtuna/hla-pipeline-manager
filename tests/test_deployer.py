"""Tests for hla_pipeline.deployer."""

from pathlib import Path

import pytest

from hla_pipeline.deployer import ResultDeployer


@pytest.fixture()
def source_dir(tmp_path: Path) -> Path:
    src = tmp_path / "source"
    src.mkdir()
    (src / "imputed.bed").write_bytes(b"\x00")
    (src / "imputed.bim").write_text("6\tHLA_DRB1_0101\t0\t32000000\tA\tG\n")
    (src / "imputed.fam").write_text("S001 S001 0 0 1 -9\n")
    return src


class TestResultDeployer:
    def test_deploy_fresh(self, source_dir: Path, tmp_path: Path) -> None:
        target = tmp_path / "production"
        deployer = ResultDeployer(target)
        report = deployer.deploy(source_dir)
        assert report.deployment_count == 3
        assert report.verified
        assert (target / "imputed.bed").exists()

    def test_deploy_with_backup(self, source_dir: Path, tmp_path: Path) -> None:
        target = tmp_path / "production"
        target.mkdir()
        (target / "old.bed").write_bytes(b"\xFF")
        deployer = ResultDeployer(target, backup_root=tmp_path / "backups")
        report = deployer.deploy(source_dir)
        assert report.deployment_count == 3
        assert len(report.files_backed_up) == 1
        assert (target / "imputed.bed").exists()

    def test_dry_run(self, source_dir: Path, tmp_path: Path) -> None:
        target = tmp_path / "production"
        deployer = ResultDeployer(target)
        report = deployer.deploy(source_dir, dry_run=True)
        assert report.deployment_count == 3
        assert not report.verified
        assert not target.exists()
