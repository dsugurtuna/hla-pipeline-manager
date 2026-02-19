"""Tests for hla_pipeline.executor."""

from pathlib import Path

import pytest

from hla_pipeline.executor import BatchExecutor, PipelineConfig


@pytest.fixture()
def fam_file(tmp_path: Path) -> Path:
    p = tmp_path / "test.fam"
    lines = [f"S{i:03d} S{i:03d} 0 0 1 -9" for i in range(1, 11)]
    p.write_text("\n".join(lines) + "\n")
    return p


class TestBatchExecutor:
    def test_count_samples(self, fam_file: Path) -> None:
        executor = BatchExecutor()
        assert executor.count_samples(fam_file) == 10

    def test_split_fam(self, fam_file: Path, tmp_path: Path) -> None:
        executor = BatchExecutor(PipelineConfig(sub_batch_size=3))
        subs = executor.split_fam(fam_file, tmp_path / "splits")
        assert len(subs) == 4  # 10 / 3 = 3 full + 1 remainder
        # First sub-batch has 3 lines
        assert len(subs[0].read_text().strip().split("\n")) == 3

    def test_execute_batch(self, fam_file: Path, tmp_path: Path) -> None:
        config = PipelineConfig(sub_batch_size=5)
        executor = BatchExecutor(config)
        result = executor.execute_batch("B001", fam_file, tmp_path / "work")
        assert result.batch_id == "B001"
        assert result.total_samples == 10
        assert result.sub_batches == 2
        assert result.success_rate == 1.0
