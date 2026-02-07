#!/bin/bash
set -euo pipefail

# One-step runner: py2disc -> disc2r (inside Docker)
#
# Usage:
#   ./run_snaptosign.sh --rna_h5ad <rna_h5ad> --atac_h5ad <atac_h5ad> --disc_out_dir <dir> --out_rds <out_rds> [--fragment_files <frag1> <frag2> ...]
#
# Backwards-compatible positional usage:
#   ./run_snaptosign.sh <rna_h5ad> <atac_h5ad> <disc_out_dir> <out_rds> [fragment_glob...]
#
# Example:
#   ./run_snaptosign.sh --rna_h5ad data/pbmc10k_multiome_rna.h5ad --atac_h5ad data/pbmc10k_multiome_atac.h5ad --disc_out_dir out --out_rds data/pbmc10k_multiome.rds

usage() {
  cat >&2 <<'EOF'
One-step runner: py2disc -> disc2r (inside Docker)

Required:
  --rna_h5ad <path>         RNA AnnData .h5ad (relative to current directory)
  --atac_h5ad <path>        ATAC AnnData .h5ad (relative to current directory)
  --disc_out_dir <dir>      Output directory for exported "disc" structure
  --out_rds <path>          Output .rds file path

Optional:
  --fragment_files <...>    One or more fragment .tsv.gz files (or globs expanded by your shell)
  --image <name:tag>        Docker image to use (overrides $SNAPTOSIGN_IMAGE)
  -h, --help                Show this help

Backwards-compatible positional form:
  run_snaptosign.sh <rna_h5ad> <atac_h5ad> <disc_out_dir> <out_rds> [fragment_glob...]
EOF
}

to_container_path() {
  local host_path="$1"

  # We only mount the current working directory into the container at /work.
  # Require all paths to be within $PWD (relative paths are fine).
  if [[ "$host_path" = /* ]]; then
    if [[ "$host_path" == "$PWD/"* ]]; then
      echo "/work/${host_path#"$PWD/"}"
      return 0
    fi
    echo "Error: absolute path '$host_path' is not under current directory '$PWD'." >&2
    echo "       Run this script from the directory you want to mount, and pass relative paths." >&2
    return 2
  fi

  echo "/work/${host_path}"
}

RNA_H5AD=""
ATAC_H5AD=""
DISC_OUT_DIR=""
OUT_RDS=""
IMAGE_OVERRIDE=""
FRAGMENT_FILES=()

if [[ "$#" -eq 0 ]]; then
  usage
  exit 2
fi

if [[ "$1" == --* ]]; then
  while [[ "$#" -gt 0 ]]; do
    case "$1" in
      --rna_h5ad)
        RNA_H5AD="${2:-}"
        shift 2
        ;;
      --atac_h5ad)
        ATAC_H5AD="${2:-}"
        shift 2
        ;;
      --disc_out_dir|--out_dir)
        DISC_OUT_DIR="${2:-}"
        shift 2
        ;;
      --out_rds|--out_file)
        OUT_RDS="${2:-}"
        shift 2
        ;;
      --image)
        IMAGE_OVERRIDE="${2:-}"
        shift 2
        ;;
      --fragment_files|--fragments)
        shift
        while [[ "$#" -gt 0 && "$1" != --* ]]; do
          FRAGMENT_FILES+=("$1")
          shift
        done
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      --)
        shift
        break
        ;;
      *)
        echo "Error: unknown option '$1'" >&2
        usage
        exit 2
        ;;
    esac
  done

  # Treat any remaining args as fragment files for convenience.
  if [[ "$#" -gt 0 ]]; then
    FRAGMENT_FILES+=("$@")
  fi
else
  # Positional (legacy) mode.
  if [[ "$#" -lt 4 ]]; then
    usage
    exit 2
  fi
  RNA_H5AD="$1"
  ATAC_H5AD="$2"
  DISC_OUT_DIR="$3"
  OUT_RDS="$4"
  shift 4
  if [[ "$#" -gt 0 ]]; then
    FRAGMENT_FILES+=("$@")
  fi
fi

if [[ -z "$RNA_H5AD" || -z "$ATAC_H5AD" || -z "$DISC_OUT_DIR" || -z "$OUT_RDS" ]]; then
  echo "Error: missing required arguments." >&2
  usage
  exit 2
fi

IMAGE="${IMAGE_OVERRIDE:-${SNAPTOSIGN_IMAGE:-ghcr.io/lukasadam/snaptosign:latest}}"

mkdir -p "$DISC_OUT_DIR"

PY_ARGS=(
  --rna_h5ad "$(to_container_path "$RNA_H5AD")"
  --atac_h5ad "$(to_container_path "$ATAC_H5AD")"
  --out_dir "$(to_container_path "$DISC_OUT_DIR")"
)

if [[ "${#FRAGMENT_FILES[@]}" -gt 0 ]]; then
  FRAG_ARGS=()
  for frag in "${FRAGMENT_FILES[@]}"; do
    FRAG_ARGS+=("$(to_container_path "$frag")")
  done
  PY_ARGS+=(--fragment_files)
  PY_ARGS+=("${FRAG_ARGS[@]}")
fi

echo "[1/2] py2disc -> ${DISC_OUT_DIR}" >&2
docker run --rm -v "$PWD":/work -w /work "$IMAGE" py2disc "${PY_ARGS[@]}"

echo "[2/2] disc2r -> ${OUT_RDS}" >&2
docker run --rm -v "$PWD":/work -w /work "$IMAGE" disc2r --disc_dir "$(to_container_path "$DISC_OUT_DIR")" --out_file "$(to_container_path "$OUT_RDS")"
