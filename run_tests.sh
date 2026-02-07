#!/bin/bash
set -euo pipefail

IMAGE="${SNAPTOSIGN_IMAGE:-ghcr.io/lukasadam/snaptosign:latest}"

echo "[1/2] Python tests (unittest)" >&2
docker run --rm -v "$PWD":/work -w /work "$IMAGE" \
  python -m unittest discover -s tests -p 'test*.py'

echo "[2/2] R tests (minimal)" >&2
docker run --rm -v "$PWD":/work -w /work "$IMAGE" \
  Rscript r/tests/run_tests.R
