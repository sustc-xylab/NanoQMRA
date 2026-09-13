from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from nanoqmra import __version__
from nanoqmra.pipeline import run_score
from nanoqmra.resources import arg_csv, rules_dir
from nanoqmra.workflow import run_workflow


def _which_conda() -> str | None:
    found = shutil.which("conda")
    if found:
        return found
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        pp = Path(prefix)
        root = pp.parent.parent if pp.parent.name == "envs" else Path(prefix)
        for cand in (root / "condabin" / "conda", root / "bin" / "conda"):
            if cand.is_file():
                return str(cand)
    return None


def _add_score_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--arg-taxa", required=True, help="ARGpore2 *.arg.w.taxa.tab")
    parser.add_argument("--plasmid-sum", required=True, help="Plascad *_plasmids_classification_sum.txt")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("--prefix", default="", help="Output file prefix (default: from --arg-taxa)")
    parser.add_argument("--rules", default="", help="Rules directory")
    parser.add_argument("--arg-csv", default="", help="ARG rank CSV")
    parser.add_argument("--total-reads", type=int, default=0, help="Optional FASTA record count for sample risk_score")


def cmd_score(args: argparse.Namespace) -> int:
    arg_taxa = Path(args.arg_taxa)
    plasmid_sum = Path(args.plasmid_sum)
    prefix = args.prefix or arg_taxa.name.split(".")[0]
    out = run_score(
        arg_taxa=arg_taxa,
        plasmid_sum=plasmid_sum,
        outdir=Path(args.outdir),
        prefix=prefix,
        total_reads=args.total_reads or None,
        rules=Path(args.rules) if args.rules else None,
        arg_csv_path=Path(args.arg_csv) if args.arg_csv else None,
    )
    print(f"[OK] mother table : {out['mother_table']}")
    print(f"[OK] risk list    : {out['risk_list']}")
    print(f"[OK] rows         : {out['n_rows']}")
    if out["risk_score"] is not None:
        print(f"[OK] risk_score   : {out['risk_score']}")
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    ok = True
    print(f"NanoQMRA {__version__}")
    print(f"python     : {sys.executable}")

    try:
        import pandas as pd
        print(f"pandas     : {pd.__version__}  OK")
    except Exception as exc:
        print(f"pandas     : MISSING ({exc})")
        ok = False

    try:
        print(f"rules      : {rules_dir()}  OK")
        print(f"ARG.csv    : {arg_csv()}  OK")
    except FileNotFoundError as exc:
        print(f"data       : MISSING ({exc})")
        ok = False

    conda = _which_conda()
    print(f"conda      : {conda or 'MISSING'}")

    env_name = args.plascad_env
    if conda:
        proc = subprocess.run(
            [conda, "run", "-n", env_name, "which", "Plascad"],
            capture_output=True,
            text=True,
        )
        found = proc.stdout.strip() if proc.returncode == 0 else ""
        print(f"Plascad    : {found or 'MISSING'}  (env={env_name})")
    else:
        print("Plascad    : SKIP (no conda)")

    argpore_env = args.argpore_env
    if conda:
        perl_ok = subprocess.run(
            [conda, "run", "-n", argpore_env, "perl", "-MBio::SeqIO", "-e", "1"],
            capture_output=True,
            text=True,
        ).returncode == 0
        par_ok = subprocess.run(
            [conda, "run", "-n", argpore_env, "which", "parallel"],
            capture_output=True,
            text=True,
        ).returncode == 0
        r_env = args.argpore_r_env
        r_ok = subprocess.run(
            [
                conda,
                "run",
                "-n",
                r_env,
                "Rscript",
                "-e",
                'stopifnot(all(sapply(c("plyr","data.table","dplyr","foreach","doParallel"), requireNamespace, quietly=TRUE)))',
            ],
            capture_output=True,
            text=True,
        ).returncode == 0
        print(f"ARGpore2 env: {argpore_env}")
        print(f"  Bio::SeqIO: {'OK' if perl_ok else 'MISSING'}")
        print(f"  parallel  : {'OK' if par_ok else 'MISSING'}")
        print(f"  R env     : {r_env}  packages {'OK' if r_ok else 'MISSING'}")
    else:
        print("ARGpore2 env: SKIP (no conda)")

    home = Path(args.argpore_home) if args.argpore_home else None
    if not home:
        print("ARGpore2   : NOT SET  (--argpore-home or ARGPORE_HOME)")
    else:
        script = home / "argpore.sh"
        print(f"ARGpore2   : {home}")
        print(f"  argpore.sh: {'OK' if script.is_file() else 'MISSING'}")
        db = home / "database"
        markers_ok = db.is_dir() and (
            (db / "markers.fasta").is_file() and (db / "markers.fasta").stat().st_size > 1000
            or any(db.glob("markers.lastindex*"))
        )
        sarg_ok = (db / "SARG_20211207_14210_filter.ffn").is_file() if db.is_dir() else False
        plsdb_ok = (db / "PLSDB_2020_03_04.fna.lastindex.prj").is_file() if db.is_dir() else False
        cen_ok = bool(list(db.glob("p_compressed+h+v*.cf"))) if db.is_dir() else False
        print(f"  markers   : {'OK' if markers_ok else 'MISSING (run setup)'}")
        print(f"  SARG      : {'OK' if sarg_ok else 'MISSING (run setup)'}")
        print(f"  PLSDB     : {'OK' if plsdb_ok else 'MISSING (run setup)'}")
        print(f"  centrifuge: {'OK' if cen_ok else 'MISSING (run setup)'}")

    print("score/check do not need ARGpore2 databases.")
    return 0 if ok else 1


def cmd_run(args: argparse.Namespace) -> int:
    home = args.argpore_home
    if not home:
        print("[ERROR] --argpore-home or ARGPORE_HOME is required for run", file=sys.stderr)
        return 2
    run_workflow(
        fasta=Path(args.fasta),
        threads=args.threads,
        outdir=Path(args.outdir),
        argpore_home=Path(home),
        plascad_env=args.plascad_env,
        argpore_env=args.argpore_env,
        argpore_r_env=args.argpore_r_env,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nanoqmra",
        description="NanoQMRA: ARGpore2 + Plascad + per-read risk scoring",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_score = sub.add_parser("score", help="Score existing ARGpore2 + Plascad tables")
    _add_score_args(p_score)
    p_score.set_defaults(func=cmd_score)

    p_check = sub.add_parser("check", help="Check scoring data and optional run-time tools")
    p_check.add_argument("--argpore-home", default=os.environ.get("ARGPORE_HOME", ""), help="ARGpore2 install directory")
    p_check.add_argument("--argpore-env", default=os.environ.get("ARGPORE_ENV", "argpore2"), help="Conda env for ARGpore2 (perl/parallel)")
    p_check.add_argument("--argpore-r-env", default=os.environ.get("ARGPORE_R_ENV", "argpore2-r"), help="Conda env for ARGpore2 R packages")
    p_check.add_argument("--plascad-env", default="Plascad", help="Conda environment that provides Plascad")
    p_check.set_defaults(func=cmd_check)

    p_run = sub.add_parser("run", help="Full workflow: ARGpore2 -> Plascad -> score")
    p_run.add_argument("-f", "--fasta", required=True, help="Input nanopore FASTA (.fa or .fasta)")
    p_run.add_argument("-t", "--threads", type=int, required=True, help="Threads for ARGpore2")
    p_run.add_argument("-o", "--outdir", required=True, help="Output directory")
    p_run.add_argument("--argpore-home", default=os.environ.get("ARGPORE_HOME", ""), help="ARGpore2 install directory")
    p_run.add_argument("--argpore-env", default=os.environ.get("ARGPORE_ENV", "argpore2"), help="Conda env for ARGpore2 (default: argpore2)")
    p_run.add_argument("--argpore-r-env", default=os.environ.get("ARGPORE_R_ENV", "argpore2-r"), help="Conda env for ARGpore2 R packages (default: argpore2-r)")
    p_run.add_argument("--plascad-env", default="Plascad", help="Conda environment that provides Plascad")
    p_run.set_defaults(func=cmd_run)
    return parser


def main(argv=None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        code = args.func(args)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    raise SystemExit(code)
