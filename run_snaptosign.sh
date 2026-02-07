#!/bin/bash
set -euo pipefail

# One-step runner: py2disc -> disc2r (inside Docker)
#
# Usage:
#   ./run_snaptosign.sh <rna_h5ad> <atac_h5ad> <disc_out_dir> <out_rds> [fragment_glob...]
#
# Example:
#   ./run_snaptosign.sh data/pbmc10k_multiome_rna.h5ad data/pbmc10k_multiome_atac.h5ad out data/pbmc10k_multiome.rds

if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <rna_h5ad> <atac_h5ad> <disc_out_dir> <out_rds> [fragment_glob...]" >&2
  echo "Example: $0 data/rna.h5ad data/atac.h5ad out data/object.rds data/fragments/*.tsv.gz" >&2
  exit 2
fi

RNA_H5AD="$1"
ATAC_H5AD="$2"
DISC_OUT_DIR="$3"
OUT_RDS="$4"
shift 4

IMAGE="${SNAPTOSIGN_IMAGE:-snaptosign:local}"

mkdir -p "$DISC_OUT_DIR"

PY_ARGS=(
  --rna_h5ad "/work/${RNA_H5AD}"
  --atac_h5ad "/work/${ATAC_H5AD}"
  --out_dir "/work/${DISC_OUT_DIR}"
)

if [ "$#" -gt 0 ]; then
  FRAG_ARGS=()
  for frag in "$@"; do
    FRAG_ARGS+=("/work/${frag}")
  done
  PY_ARGS+=(--fragment_files)
  PY_ARGS+=("${FRAG_ARGS[@]}")
fi

echo "[1/2] py2disc -> ${DISC_OUT_DIR}" >&2
docker run --rm -v "$PWD":/work -w /work "$IMAGE" py2disc "${PY_ARGS[@]}"

echo "[2/2] disc2r -> ${OUT_RDS}" >&2
docker run --rm -v "$PWD":/work -w /work "$IMAGE" disc2r --disc_dir "/work/${DISC_OUT_DIR}" --out_file "/work/${OUT_RDS}"
