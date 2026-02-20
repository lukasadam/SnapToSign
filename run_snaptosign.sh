#!/bin/bash
set -euo pipefail

# One-step runner: py2disc -> disc2r (inside Apptainer/Singularity)
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
One-step runner: py2disc -> disc2r (inside Apptainer/Singularity)

Required:
  --rna_h5ad <path>         RNA AnnData .h5ad (relative or absolute)
  --atac_h5ad <path>        ATAC AnnData .h5ad (relative or absolute)
  --disc_out_dir <dir>      Output directory for exported "disc" structure
  --out_rds <path>          Output .rds file path (relative or absolute)

Optional:
  --fragment_files <...>    One or more fragment .tsv.gz files (or globs expanded by your shell)
  --image <ref>             Apptainer image ref (local .sif or docker://...; overrides $SNAPTOSIGN_IMAGE)
  -h, --help                Show this help

Backwards-compatible positional form:
  run_snaptosign.sh <rna_h5ad> <atac_h5ad> <disc_out_dir> <out_rds> [fragment_glob...]
EOF
}

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN=python3
  else
    PYTHON_BIN=python
  fi
fi

abspath() {
  "$PYTHON_BIN" -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$1"
}

APPTAINER_BIN="${APPTAINER_BIN:-${SINGULARITY_BIN:-apptainer}}"

DOCKER_MOUNTS=()
MOUNT_HOST_DIRS=()
MOUNT_CONTAINER_DIRS=()

ensure_mount_dir() {
  local host_dir="$1"
  local container_dir=""

  for i in "${!MOUNT_HOST_DIRS[@]}"; do
    if [[ "${MOUNT_HOST_DIRS[$i]}" == "$host_dir" ]]; then
      echo "${MOUNT_CONTAINER_DIRS[$i]}"
      return 0
    fi
  done

  local idx="${#MOUNT_HOST_DIRS[@]}"
  container_dir="/mnt/${idx}"
  MOUNT_HOST_DIRS+=("$host_dir")
  MOUNT_CONTAINER_DIRS+=("$container_dir")
  DOCKER_MOUNTS+=("--bind" "${host_dir}:${container_dir}:rw")
  echo "$container_dir"
}

normalize_image_ref() {
  local ref="$1"

  # If it is a local file path (e.g. ./snaptosign.sif or /path/snaptosign.sif), keep as-is.
  if [[ "$ref" == /* || "$ref" == ./* || "$ref" == ../* ]]; then
    echo "$ref"
    return 0
  fi
  if [[ -f "$ref" ]]; then
    echo "$ref"
    return 0
  fi
  if [[ "$ref" == *.sif ]]; then
    echo "$ref"
    return 0
  fi

  # If it already has a transport (docker://, oras://, etc), keep as-is.
  if [[ "$ref" == *"://"* ]]; then
    echo "$ref"
    return 0
  fi

  # Otherwise assume it's a docker registry reference.
  echo "docker://${ref}"
}

to_container_path() {
  local host_path="$1"

  # Always resolve to an absolute host path so we can bind/mount reliably.
  local abs_host_path
  abs_host_path="$(abspath "$host_path")"

  # We only mount the current working directory into the container at /work.
  # Require all paths to be within $PWD (relative paths are fine).
  if [[ "$abs_host_path" = /* ]]; then
    if [[ "$abs_host_path" == "$PWD/"* ]]; then
      echo "/work/${abs_host_path#"$PWD/"}"
      return 0
    fi

    # Absolute path outside $PWD: mount its parent directory and translate to /mnt/N/<basename>.
    local host_dir
    host_dir="$(dirname "$abs_host_path")"
    local container_dir
    container_dir="$(ensure_mount_dir "$host_dir")"
    echo "${container_dir}/$(basename "$abs_host_path")"
    return 0
  fi

  # (Should be unreachable because abspath always returns absolute.)
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

DEFAULT_IMAGE="docker://ghcr.io/lukasadam/snaptosign:latest"
if [[ -z "${SNAPTOSIGN_IMAGE:-}" && -z "$IMAGE_OVERRIDE" ]]; then
  if [[ -f "./snaptosign.sif" ]]; then
    DEFAULT_IMAGE="./snaptosign.sif"
  fi
fi

IMAGE_RAW="${IMAGE_OVERRIDE:-${SNAPTOSIGN_IMAGE:-$DEFAULT_IMAGE}}"
IMAGE="$(normalize_image_ref "$IMAGE_RAW")"

mkdir -p "$DISC_OUT_DIR"
mkdir -p "$(dirname "$(abspath "$OUT_RDS")")"

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
"$APPTAINER_BIN" exec --bind "$PWD":/work --pwd /work "${DOCKER_MOUNTS[@]}" "$IMAGE" py2disc "${PY_ARGS[@]}"

echo "[2/2] disc2r -> ${OUT_RDS}" >&2
DISC_DIR_CONTAINER="$(to_container_path "$DISC_OUT_DIR")"
OUT_RDS_CONTAINER="$(to_container_path "$OUT_RDS")"
"$APPTAINER_BIN" exec --bind "$PWD":/work --pwd /work "${DOCKER_MOUNTS[@]}" "$IMAGE" disc2r \
  --disc_dir "$DISC_DIR_CONTAINER" \
  --out_file "$OUT_RDS_CONTAINER"