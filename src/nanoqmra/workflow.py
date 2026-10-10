from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from nanoqmra.pipeline import run_score


FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".fna", ".ffn"}


def count_fasta_records(path: Path) -> int:
    n = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                n += 1
    return n


def place_fasta(src: Path, dest: Path) -> None:
    """Make dest available without duplicating a large FASTA when possible.

    A hard link is another name for the same file. It does not follow the
    source the way a symbolic link does, and removing dest does not remove
    the original. If the two paths are on different filesystems, fall back
    to a copy. Plascad looks up ``<prefix>.fasta`` beside its other outputs,
    so dest should use that suffix even when src is ``.fa`` or ``.fna``.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        if dest.resolve() == src.resolve() and not dest.is_symlink():
            return
        dest.unlink()
    try:
        os.link(src, dest)
    except OSError:
        shutil.copy2(src, dest)


def prepare_fasta(src: Path, dest_dir: Path, prefix: str) -> Path:
    suffix = src.suffix.lower()
    if suffix not in FASTA_SUFFIXES:
        accepted = ", ".join(sorted(FASTA_SUFFIXES))
        raise ValueError(f"unsupported FASTA suffix {suffix!r}; accepted: {accepted}")
    dest = dest_dir / f"{prefix}.fasta"
    place_fasta(src, dest)
    return dest


def find_argpore_taxa(argpore_cwd: Path) -> Path:
    dirs = sorted(
        [
            p
            for p in argpore_cwd.iterdir()
            if p.is_dir()
            and (
                "_ARGpore2_" in p.name
                or "_ARGpore_" in p.name
            )
        ]
    )
    if not dirs:
        raise FileNotFoundError(f"ARGpore2 output directory not found in {argpore_cwd}")
    out_dir = dirs[-1]
    for pattern in ("*.arg.w.taxa.tab", "*.w.taxa.tab"):
        hits = sorted(out_dir.glob(pattern))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"ARG taxa table not found in {out_dir}")


def find_plascad_sum(plascad_cwd: Path, prefix: str) -> Path:
    hits = sorted(plascad_cwd.glob(f"{prefix}*_plasmids_classification_sum.txt"))
    if not hits:
        hits = sorted(plascad_cwd.glob("*_plasmids_classification_sum.txt"))
    if not hits:
        raise FileNotFoundError(f"Plascad summary not found in {plascad_cwd}")
    path = hits[-1]
    if path.stat().st_size == 0:
        raise RuntimeError(f"Plascad summary is empty: {path}")
    return path


def _conda_exe() -> str:
    found = shutil.which("conda")
    if found:
        return found
    root = _conda_root()
    if root:
        for cand in (root / "condabin" / "conda", root / "bin" / "conda"):
            if cand.is_file():
                return str(cand)
    raise RuntimeError("conda not found on PATH; activate an env or add condabin")


def _conda_run_env() -> dict:
    """Drop PERL5LIB so conda env perl (Bio::SeqIO) is not shadowed."""
    env = os.environ.copy()
    env.pop("PERL5LIB", None)
    return env


def _conda_root() -> Path | None:
    which = shutil.which("conda")
    if which:
        exe = Path(which).resolve()
        if exe.parent.name in {"bin", "condabin"}:
            return exe.parent.parent
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        pp = Path(prefix)
        if pp.parent.name == "envs":
            return pp.parent.parent
        return pp
    return None


def _conda_run_path(env_name: str, *extra_bins: str) -> str:
    """Put the target env's bin ahead of base.

    ``conda run`` on this machine leaves ``miniforge3/bin`` first, so
    ``perl`` misses Bio::SeqIO installed in the env.
    """
    parts: list[str] = []
    bindir = _conda_env_bin(env_name)
    if bindir is not None:
        parts.append(str(bindir))
    parts.extend(b for b in extra_bins if b)
    parts.append("$PATH")
    return ":".join(parts)


def _conda_env_bin(env_name: str) -> Path | None:
    root = _conda_root()
    if root is None:
        return None
    bindir = root / "envs" / env_name / "bin"
    return bindir if bindir.is_dir() else None


def run_argpore(
    argpore_home: Path,
    fasta: Path,
    threads: int,
    workdir: Path,
    log: Path,
    env_name: str,
    r_env: str = "",
) -> None:
    script = argpore_home / "argpore.sh"
    if not script.is_file():
        raise FileNotFoundError(f"argpore.sh not found: {argpore_home}")
    workdir.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    r_bin = _conda_env_bin(r_env) if r_env else None
    r_prefix = str(r_bin) if r_bin and (r_bin / "Rscript").is_file() else ""
    inner = (
        "set -eu; unset PERL5LIB; "
        + f'export PATH="{_conda_run_path(env_name, r_prefix)}"; '
        + f'export ARGPORE_HOME="{argpore_home}"; '
        + f'exec argpore2 -f "{fasta}" -t "{threads}" -o "{workdir}"'
    )
    with log.open("w") as handle:
        subprocess.run(
            [
                _conda_exe(),
                "run",
                "-n",
                env_name,
                "--no-capture-output",
                "bash",
                "-c",
                inner,
            ],
            cwd=workdir,
            env=_conda_run_env(),
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )


def run_plascad(fasta: Path, env_name: str, workdir: Path, log: Path) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    local_fa = workdir / fasta.name
    if (
        not local_fa.exists()
        or local_fa.is_symlink()
        or local_fa.resolve() != fasta.resolve()
    ):
        place_fasta(fasta, local_fa)
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as handle:
        subprocess.run(
            [
                _conda_exe(),
                "run",
                "-n",
                env_name,
                "--no-capture-output",
                "bash",
                "-c",
                f'export PATH="{_conda_run_path(env_name)}"; exec Plascad -i "{local_fa.name}"',
            ],
            cwd=workdir,
            env=_conda_run_env(),
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )


def run_workflow(
    fasta: Path,
    threads: int,
    outdir: Path,
    argpore_home: Path,
    plascad_env: str,
    argpore_env: str,
    argpore_r_env: str = "",
) -> dict:
    fasta = fasta.resolve()
    outdir = outdir.resolve()
    argpore_home = argpore_home.resolve()
    if not fasta.is_file():
        raise FileNotFoundError(f"input fasta not found: {fasta}")

    prefix = fasta.stem
    logs = outdir / "logs"
    results = outdir / "results"
    arg_dir = outdir / "intermediate" / "argpore2"
    plas_dir = outdir / "intermediate" / "plascad"
    for path in (logs, results, arg_dir, plas_dir):
        path.mkdir(parents=True, exist_ok=True)

    work_fa = prepare_fasta(fasta, outdir / "input", prefix)

    print(f"[1/4] ARGpore2  {argpore_home}  env={argpore_env}")
    run_argpore(
        argpore_home,
        work_fa,
        threads,
        arg_dir,
        logs / f"{prefix}_ARGpore.log",
        env_name=argpore_env,
        r_env=argpore_r_env,
    )
    arg_taxa = find_argpore_taxa(arg_dir)
    print(f"[INFO] ARG taxa: {arg_taxa}")

    print(f"[2/4] Plascad env={plascad_env}")
    run_plascad(work_fa, plascad_env, plas_dir, logs / f"{prefix}_Plascad.log")
    plasmid_sum = find_plascad_sum(plas_dir, prefix)
    n_lines = sum(1 for _ in plasmid_sum.open())
    if n_lines < 1:
        raise RuntimeError("Plascad produced no classifications")
    print(f"[INFO] Plascad summary: {plasmid_sum} ({n_lines} rows)")

    print("[3/4] Scoring")
    score_out = run_score(
        arg_taxa=arg_taxa,
        plasmid_sum=plasmid_sum,
        outdir=results,
        prefix=prefix,
        total_reads=count_fasta_records(fasta),
    )

    print("[4/4] Done")
    print(f"[OK] mother table : {score_out['mother_table']}")
    print(f"[OK] risk list    : {score_out['risk_list']}")
    print(f"[OK] risk_score   : {score_out['risk_score']}")
    return score_out
