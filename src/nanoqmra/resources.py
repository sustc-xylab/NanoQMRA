from pathlib import Path


def package_data_dir() -> Path:
    return Path(__file__).resolve().parent / "data"


def rules_dir() -> Path:
    bundled = package_data_dir() / "rules"
    if (bundled / "who_rules.tsv").is_file():
        return bundled
    raise FileNotFoundError("Cannot find packaged rules/who_rules.tsv")


def arg_csv() -> Path:
    bundled = package_data_dir() / "reference" / "ARG.csv"
    if bundled.is_file():
        return bundled
    raise FileNotFoundError("Cannot find packaged ARG.csv")
