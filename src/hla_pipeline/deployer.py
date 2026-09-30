"""Backup-first deployment of imputed results.

Copies result files into a release directory after taking a timestamped
backup of whatever is already there, then verifies each copy by SHA-256.
"""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path


def sha256(path: Path) -> str:
    """SHA-256 hex digest of a file, read in 1 MiB chunks."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass
class DeploymentReport:
    """Summary of a deployment."""

    source_dir: str
    target_dir: str
    backup_dir: str = ""
    files_deployed: list[str] = field(default_factory=list)
    files_backed_up: list[str] = field(default_factory=list)
    checksums: dict[str, str] = field(default_factory=dict)
    mismatched: list[str] = field(default_factory=list)
    verified: bool = False
    dry_run: bool = False

    @property
    def deployment_count(self) -> int:
        return len(self.files_deployed)


class ResultDeployer:
    """Deploy result files with a backup-first strategy.

    Parameters
    ----------
    target_dir : path
        Release directory to deploy into.
    backup_root : path, optional
        Root for backups. Defaults to ``<target_dir>/../backups``.
    extensions : list of str, optional
        File extensions to deploy. Defaults to the PLINK binary set.
    """

    DEFAULT_EXTENSIONS = (".bed", ".bim", ".fam")

    def __init__(
        self,
        target_dir: str | Path,
        backup_root: str | Path | None = None,
        extensions: list[str] | None = None,
    ) -> None:
        self.target_dir = Path(target_dir)
        self.backup_root = (
            Path(backup_root) if backup_root else self.target_dir.parent / "backups"
        )
        self.extensions = list(extensions or self.DEFAULT_EXTENSIONS)

    def _existing_files(self) -> list[Path]:
        if not self.target_dir.exists():
            return []
        found: list[Path] = []
        for ext in self.extensions:
            found.extend(sorted(self.target_dir.glob(f"*{ext}")))
        return found

    def _create_backup(self, files: list[Path]) -> Path:
        # Microseconds keep two deployments in the same second from sharing
        # (and overwriting) one backup directory.
        ts = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        backup_dir = self.backup_root / ts
        backup_dir.mkdir(parents=True, exist_ok=False)
        for f in files:
            shutil.copy2(f, backup_dir / f.name)
        return backup_dir

    def deploy(self, source_dir: str | Path, dry_run: bool = False) -> DeploymentReport:
        """Back up the target, copy new files in, and verify the copies.

        With ``dry_run=True`` nothing is written; the report lists what would
        be backed up and deployed.
        """
        source = Path(source_dir)
        report = DeploymentReport(
            source_dir=str(source), target_dir=str(self.target_dir), dry_run=dry_run
        )

        to_deploy: list[Path] = []
        for ext in self.extensions:
            to_deploy.extend(sorted(source.glob(f"*{ext}")))
        if not to_deploy:
            return report

        existing = self._existing_files()
        report.files_backed_up = [f.name for f in existing]
        report.files_deployed = [f.name for f in to_deploy]
        report.checksums = {f.name: sha256(f) for f in to_deploy}
        if dry_run:
            return report

        if existing:
            report.backup_dir = str(self._create_backup(existing))
        self.target_dir.mkdir(parents=True, exist_ok=True)
        for f in to_deploy:
            shutil.copy2(f, self.target_dir / f.name)

        report.mismatched = [
            name
            for name, expected in report.checksums.items()
            if not (self.target_dir / name).is_file()
            or sha256(self.target_dir / name) != expected
        ]
        report.verified = not report.mismatched
        return report
