"""Backup-first deployment of imputed results.

Copies verified SNP2HLA result files into a release directory after taking a
timestamped backup of any existing release files, checks every copy by
SHA-256, rolls back on mismatch, cleans stale release files, and writes a
``SHA256SUMS.tsv`` manifest.
"""

from __future__ import annotations

import hashlib
import os
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from .verifier import OUTPUT_SUFFIXES, ImputationVerifier

MANIFEST_FILENAME = "SHA256SUMS.tsv"


def sha256(path: Path) -> str:
    """SHA-256 hex digest of a file, read in 1 MiB chunks."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_copy(src: Path, dst: Path) -> None:
    tmp = dst.parent / f".{dst.name}.tmp.{os.getpid()}"
    try:
        shutil.copy2(src, tmp)
        os.replace(tmp, dst)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _atomic_write_text(dst: Path, content: str) -> None:
    tmp = dst.parent / f".{dst.name}.tmp.{os.getpid()}"
    try:
        tmp.write_text(content, encoding="utf-8")
        os.replace(tmp, dst)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _restore_file(src: Path, dst: Path) -> None:
    tmp = dst.parent / f".{dst.name}.restore.{os.getpid()}"
    try:
        shutil.copyfile(src, tmp)
        shutil.copystat(src, tmp)
        os.replace(tmp, dst)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


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
    manifest_path: str = ""
    rolled_back: bool = False

    @property
    def deployment_count(self) -> int:
        return len(self.files_deployed)


class ResultDeployer:
    """Deploy result files with a verification gate and backup-first rollback.

    Parameters
    ----------
    target_dir : path
        Release directory to deploy into.
    backup_root : path, optional
        Root for backups. Defaults to ``<target_dir>/../backups``.
    extensions : sequence of str, optional
        File suffixes to deploy. Defaults to
        :data:`hla_pipeline.verifier.OUTPUT_SUFFIXES`.
    require_verified : bool, default True
        Run :class:`~hla_pipeline.verifier.ImputationVerifier` on
        ``source_dir`` before deploying.
    """

    DEFAULT_EXTENSIONS = OUTPUT_SUFFIXES

    def __init__(
        self,
        target_dir: str | Path,
        backup_root: str | Path | None = None,
        extensions: Sequence[str] | None = None,
        require_verified: bool = True,
    ) -> None:
        self.target_dir = Path(target_dir)
        self.backup_root = (
            Path(backup_root) if backup_root else self.target_dir.parent / "backups"
        )
        exts = extensions if extensions is not None else self.DEFAULT_EXTENSIONS
        self.extensions = list(dict.fromkeys(exts))
        if not self.extensions:
            raise ValueError("extensions must not be empty")
        self.require_verified = require_verified

    def _collect_matching_files(self, directory: Path) -> list[Path]:
        if not directory.is_dir():
            return []
        seen: set[Path] = set()
        found: list[Path] = []
        for ext in self.extensions:
            for path in sorted(directory.glob(f"*{ext}")):
                if path.is_file() and path not in seen:
                    seen.add(path)
                    found.append(path)
        return found

    def _existing_files(self) -> list[Path]:
        return self._collect_matching_files(self.target_dir)

    def _create_backup(self, files: list[Path]) -> Path:
        # Microseconds keep two deployments in the same second from sharing
        # (and overwriting) one backup directory.
        ts = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        backup_dir = self.backup_root / ts
        backup_dir.mkdir(parents=True, exist_ok=False)
        for f in files:
            _restore_file(f, backup_dir / f.name)
        return backup_dir

    def _rollback(self, report: DeploymentReport) -> None:
        for name in report.files_deployed:
            (self.target_dir / name).unlink(missing_ok=True)
        if report.backup_dir:
            backup_path = Path(report.backup_dir)
            for backed_up in sorted(backup_path.iterdir()):
                if backed_up.is_file():
                    _restore_file(backed_up, self.target_dir / backed_up.name)
        report.rolled_back = True

    def _write_manifest(self, checksums: Mapping[str, str]) -> Path:
        manifest_path = self.target_dir / MANIFEST_FILENAME
        lines = ["sha256\tfile"]
        for name in sorted(checksums):
            lines.append(f"{checksums[name]}\t{name}")
        _atomic_write_text(manifest_path, "\n".join(lines) + "\n")
        return manifest_path

    @staticmethod
    def verify_manifest(target_dir: str | Path) -> list[str]:
        """Verify files in ``target_dir`` against ``SHA256SUMS.tsv``.

        Returns a list of missing or checksum-mismatched filenames (empty if
        all entries verify cleanly).
        """
        target = Path(target_dir)
        manifest_path = target / MANIFEST_FILENAME
        if not manifest_path.is_file():
            raise FileNotFoundError(f"Manifest not found: {manifest_path}")
        lines = [
            ln.strip()
            for ln in manifest_path.read_text(encoding="utf-8").splitlines()
            if ln.strip()
        ]
        if not lines or lines[0] != "sha256\tfile" or len(lines) < 2:
            raise ValueError(f"Malformed or empty manifest: {manifest_path}")
        mismatched: list[str] = []
        for entry in lines[1:]:
            parts = entry.split("\t", 1)
            if len(parts) != 2:
                raise ValueError(f"Malformed manifest line: {entry!r}")
            expected_hash, name = parts
            candidate = target / name
            if not candidate.is_file() or sha256(candidate) != expected_hash:
                mismatched.append(name)
        return mismatched

    def deploy(
        self,
        source_dir: str | Path,
        dry_run: bool = False,
        require_verified: bool | None = None,
    ) -> DeploymentReport:
        """Verify, back up target, copy new files atomically, and verify SHA-256.

        With ``dry_run=True`` nothing is written; the report lists what would
        be backed up and deployed. On checksum mismatch after copying, deployed
        files are rolled back from the backup directory.
        """
        source = Path(source_dir)
        if not source.is_dir():
            raise NotADirectoryError(f"Source directory not found: {source}")

        verify_gate = (
            self.require_verified if require_verified is None else require_verified
        )
        if verify_gate:
            vreport = ImputationVerifier().verify_batch(source)
            if not vreport.all_ok:
                raise ValueError(
                    "Source directory failed imputation verification:\n"
                    f"{ImputationVerifier.format_report(vreport)}"
                )

        report = DeploymentReport(
            source_dir=str(source), target_dir=str(self.target_dir), dry_run=dry_run
        )

        to_deploy = self._collect_matching_files(source)
        if not to_deploy:
            return report

        existing = self._existing_files()
        existing_manifest = self.target_dir / MANIFEST_FILENAME
        backup_candidates = list(existing)
        if existing_manifest.is_file():
            backup_candidates.append(existing_manifest)

        report.files_backed_up = [f.name for f in existing]
        report.files_deployed = [f.name for f in to_deploy]
        report.checksums = {f.name: sha256(f) for f in to_deploy}
        if dry_run:
            return report

        if backup_candidates:
            report.backup_dir = str(self._create_backup(backup_candidates))
        self.target_dir.mkdir(parents=True, exist_ok=True)

        for f in to_deploy:
            _atomic_copy(f, self.target_dir / f.name)

        report.mismatched = [
            name
            for name, expected in report.checksums.items()
            if not (self.target_dir / name).is_file()
            or sha256(self.target_dir / name) != expected
        ]
        if report.mismatched:
            self._rollback(report)
            report.verified = False
            return report

        deployed_names = set(report.files_deployed)
        for old_file in existing:
            if old_file.name not in deployed_names:
                old_file.unlink(missing_ok=True)

        manifest_path = self._write_manifest(report.checksums)
        report.manifest_path = str(manifest_path)
        report.verified = True
        return report
