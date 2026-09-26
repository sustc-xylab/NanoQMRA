from pathlib import Path

import pytest

from nanoqmra.workflow import prepare_fasta


def test_prepare_fasta_hardlinks_other_suffix(tmp_path: Path):
    src = tmp_path / "reads.fa"
    src.write_text(">r1\nACGT\n")
    dest_dir = tmp_path / "out"
    dest = prepare_fasta(src, dest_dir, "reads")
    assert dest.name == "reads.fasta"
    assert dest.read_text() == src.read_text()
    assert dest.stat().st_ino == src.stat().st_ino


def test_prepare_fasta_rejects_unknown_suffix(tmp_path: Path):
    src = tmp_path / "reads.fastq"
    src.write_text("@r1\nACGT\n+\n!!!!\n")
    with pytest.raises(ValueError):
        prepare_fasta(src, tmp_path / "out", "reads")
