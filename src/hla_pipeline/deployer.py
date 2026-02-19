"""Safe deployment module.

Deploys imputed HLA results to production with automatic timestamped
backup of existing files.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List


@dataclass
class DeploymentReport:
    """Summary of a deployment operation."""

    source_dir: str
    target_dir: str
    backup_dir: str = ""
    files_deployed: List[str] = field(default_factory=list)
    files_backed_up: List[str] = field(default_factory=list)
    verified: bool = False

    @property
    def deployment_count(self) -> int:
        return len(self.files_deployed)


class ResultDeployer:
    """Deploy imputed results with backup-first strategy.

    Creates a timestamped backup of existing production files before
    copying new results into the target directory.

    Parameters
    ----------
    target_dir : str or Path
        Production directory to deploy into.
    backup_root : str or Path, optional
        Root directory for backups. Defaults to ``target_dir/../backups``.
    extensions : list of str
        File extensions to deploy. Defaults to PLINK binary set.
    """

    DEFAULT_EXTENSIONS = [".bed", ".bim", ".fam"]

    def __init__(
        self,
        target_dir: str | Path,
        backup_root: str | Path | None = None,
        extensions: List[str] | None = None,
    ) -> None:
        self.target_dir = Path(target_dir)
        self.backup_root = (
            Path(backup_root) if backup_root else self.target_dir.parent / "backups"
        )
        self.extensions = extensions or self.DEFAULT_EXTENSIONS

    def _create_backup(self) -> Path:
        """Create a timestamped backup of current production files."""
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_dir = self.backup_root / ts
        backup_dir.mkdir(parents=True, exist_ok=True)

        backed_up: List[str] = []
        for ext in self.extensions:
            for f in self.target_dir.glob(f"*{ext}"):
                dest = backup_dir / f.name
                shutil.copy2(f, dest)
                backed_up.append(f.name)

        return backup_dir

    def deploy(
        self,
        source_dir: str | Path,
        dry_run: bool = False,
    ) -> DeploymentReport:
        """Deploy results from source to production.

        Parameters
        ----------
        source_dir : path
            Directory containing the new result files.
        dry_run : bool
            If True, report what would be done without copying.
        """
        source = Path(source_dir)
        report = DeploymentReport(
            source_dir=str(source),
            target_dir=str(self.target_dir),
        )

        # Identify files to deploy
        files_to_deploy: List[Path] = []
        for ext in self.extensions:
            files_to_deploy.extend(source.glob(f"*{ext}"))

        if not files_to_deploy:
            return report

        # Backup existing
        if not dry_run and self.target_dir.exists():
            backup_dir = self._create_backup()
            report.backup_dir = str(backup_dir)
            report.files_backed_up = [f.name for f in backup_dir.iterdir()]

        # Deploy
        self.target_dir.mkdir(parents=True, exist_ok=True)
        for f in files_to_deploy:
            if not dry_run:
                shutil.copy2(f, self.target_dir / f.name)
            report.files_deployed.append(f.name)

        # Verify
        if not dry_run:
            report.verified = all(
                (self.target_dir / fname).exists() for fname in report.files_deployed
            )

        return report
