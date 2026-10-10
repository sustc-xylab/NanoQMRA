"""NanoQMRA: nanopore workflow for species-resolved microbial risk scoring."""

try:
    from importlib.metadata import version as _pkg_version

    __version__ = _pkg_version("nanoqmra")
except Exception:
    __version__ = "1.2.0"
