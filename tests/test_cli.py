from nanoqmra.cli import build_parser, cmd_check


def test_parser_has_score_check_run():
    parser = build_parser()
    args = parser.parse_args(["check"])
    assert args.cmd == "check"
    args = parser.parse_args(
        ["score", "--arg-taxa", "a.tsv", "--plasmid-sum", "b.txt", "-o", "out"]
    )
    assert args.cmd == "score"


def test_check_scoring_deps_ok(capsys):
    parser = build_parser()
    args = parser.parse_args(["check"])
    code = cmd_check(args)
    captured = capsys.readouterr().out
    assert code == 0
    assert "pandas" in captured
    assert "rules" in captured
