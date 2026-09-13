from pathlib import Path

from nanoqmra.pipeline import run_score
from nanoqmra.risk_list import build_risk_score_list
from nanoqmra.score import ARG_RANK_SCORE, build_mother_table, match_arg_rank, score_mother_table
from nanoqmra.resources import arg_csv, rules_dir

FIXTURES = Path(__file__).parent / "fixtures"


def test_kpc_conjugative_score_is_500():
    mother = build_mother_table(
        FIXTURES / "arg.w.taxa.tab",
        FIXTURES / "plasmids_classification_sum.txt",
        rules_dir(),
    )
    scored = score_mother_table(mother, arg_csv())
    row = scored.set_index("query").loc["read1"]
    assert row["WHO_priority"] == "Critical Priority"
    assert int(row["WHO_priority_score"]) == 10
    assert int(row["ARG_rank_score"]) == 5
    assert int(row["mobile_factor"]) == 2
    assert int(row["conjugative_plasmid_factor"]) == 5
    assert int(row["risk_score"]) == 500


def test_plasmid_display_keeps_host():
    mother = build_mother_table(
        FIXTURES / "arg.w.taxa.tab",
        FIXTURES / "plasmids_classification_sum.txt",
        rules_dir(),
    )
    scored = score_mother_table(mother, arg_csv())
    risk = build_risk_score_list(scored)
    row = risk[risk["risk_score"] == 500].iloc[0]
    assert row["species"] == "plasmid"
    assert row["host"] == "Klebsiella pneumoniae"


def test_tetm_duplicate_rank_uses_higher_score():
    from nanoqmra.score import load_arg_rank_map

    mapping = load_arg_rank_map(arg_csv())
    rank, method = match_arg_rank("tetM", mapping)
    assert rank == "Rank I"
    assert method == "exact"
    assert ARG_RANK_SCORE[rank] == 5


def test_missing_query_column_raises(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("subtype\ttype\nKPC-2\tbeta-lactams\n")
    plas = FIXTURES / "plasmids_classification_sum.txt"
    try:
        build_mother_table(bad, plas, rules_dir())
        raised = False
    except ValueError as exc:
        raised = "query" in str(exc)
    assert raised


def test_score_cli_pipeline(tmp_path):
    out = run_score(
        arg_taxa=FIXTURES / "arg.w.taxa.tab",
        plasmid_sum=FIXTURES / "plasmids_classification_sum.txt",
        outdir=tmp_path,
        prefix="demo",
        total_reads=10,
    )
    assert out["n_rows"] == 5
    assert out["risk_score"] is not None
    text = Path(out["risk_list"]).read_text()
    assert "risk_score:" in text
