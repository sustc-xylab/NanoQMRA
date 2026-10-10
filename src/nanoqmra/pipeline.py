from __future__ import annotations

from pathlib import Path

import pandas as pd

from nanoqmra.resources import arg_csv, rules_dir
from nanoqmra.risk_list import append_sample_risk_score, build_risk_score_list
from nanoqmra.score import build_mother_table, score_mother_table


def run_score(
    arg_taxa: Path,
    plasmid_sum: Path,
    outdir: Path,
    prefix: str = "sample",
    total_reads: int | None = None,
    rules: Path | None = None,
    arg_csv_path: Path | None = None,
) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rules = Path(rules) if rules else rules_dir()
    arg_csv_path = Path(arg_csv_path) if arg_csv_path else arg_csv()

    mother = build_mother_table(arg_taxa, plasmid_sum, rules)
    scored = score_mother_table(mother, arg_csv_path)
    mother_out = outdir / f"{prefix}_mother_table.scored.tsv"
    scored.to_csv(mother_out, sep="\t", index=False)

    risk = build_risk_score_list(scored)
    risk_out = outdir / f"{prefix}_risk_score_list.tsv"
    risk.to_csv(risk_out, sep="\t", index=False)

    sample_score = None
    if total_reads:
        if total_reads <= 0:
            raise ValueError("total_reads must be > 0")
        sample_score = append_sample_risk_score(risk_out, total_reads)

    return {
        "mother_table": mother_out,
        "risk_list": risk_out,
        "n_rows": int(len(scored)),
        "risk_score": sample_score,
    }
