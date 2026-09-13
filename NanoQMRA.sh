#!/usr/bin/env bash
set -euo pipefail

show_help() {
cat << EOF
Usage: ${0##*/} -f reads.fasta -t threads -o outdir [-p Plascad_env] [--argpore-home DIR] [--argpore-env NAME]

version 1.2.0

This is a compatibility wrapper around:
  nanoqmra run -f ... -t ... -o ... --plascad-env ... --argpore-home ...

Prefer installing the package:
  pip install -e .
  nanoqmra check
  nanoqmra score --arg-taxa ... --plasmid-sum ... -o outdir
  nanoqmra run -f sample.fasta -t 8 -o outdir --argpore-home /path/to/ARGpore2
EOF
}

INPUT_FA=""
THREADS=""
OUTDIR=""
PLASCAD_ENV="Plascad"
ARGPORE_HOME="${ARGPORE_HOME:-}"
ARGPORE_ENV="${ARGPORE_ENV:-argpore2}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      show_help
      exit 0
      ;;
    -f)
      INPUT_FA="$2"; shift 2
      ;;
    -t)
      THREADS="$2"; shift 2
      ;;
    -o)
      OUTDIR="$2"; shift 2
      ;;
    -p)
      PLASCAD_ENV="$2"; shift 2
      ;;
    --argpore-home)
      ARGPORE_HOME="$2"; shift 2
      ;;
    --argpore-env)
      ARGPORE_ENV="$2"; shift 2
      ;;
    *)
      echo "[ERROR] Invalid option: $1" >&2
      show_help >&2
      exit 1
      ;;
  esac
done

if [[ -z "${INPUT_FA}" || -z "${THREADS}" ]]; then
  echo "[ERROR] -f and -t are required" >&2
  exit 1
fi

if [[ -z "${OUTDIR}" ]]; then
  prefix="$(basename "${INPUT_FA}")"
  prefix="${prefix%.*}"
  OUTDIR="./nanoqmra_out/${prefix}"
  echo "[INFO] -o not set, using ${OUTDIR}"
fi

if [[ -z "${ARGPORE_HOME}" ]]; then
  echo "[ERROR] set --argpore-home or ARGPORE_HOME" >&2
  exit 1
fi

if command -v nanoqmra >/dev/null 2>&1; then
  exec nanoqmra run -f "${INPUT_FA}" -t "${THREADS}" -o "${OUTDIR}" \
    --argpore-home "${ARGPORE_HOME}" --argpore-env "${ARGPORE_ENV}" --plascad-env "${PLASCAD_ENV}"
fi

SCRIPT=$(realpath "$0")
DIR=$(dirname "$SCRIPT")
export PYTHONPATH="${DIR}/src${PYTHONPATH:+:$PYTHONPATH}"
exec python3 -m nanoqmra run -f "${INPUT_FA}" -t "${THREADS}" -o "${OUTDIR}" \
  --argpore-home "${ARGPORE_HOME}" --argpore-env "${ARGPORE_ENV}" --plascad-env "${PLASCAD_ENV}"
