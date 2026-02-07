#!/bin/bash
set -euo pipefail

if [[ -n "${SNAPTOSIGN_IMAGE:-}" ]]; then
  IMAGE="$SNAPTOSIGN_IMAGE"
elif docker image inspect snaptosign:local >/dev/null 2>&1; then
  IMAGE="snaptosign:local"
else
  IMAGE="ghcr.io/lukasadam/snaptosign:latest"
fi

echo "Using image: $IMAGE" >&2

echo "[1/2] Python tests (unittest)" >&2
docker run --rm -v "$PWD":/work -w /work "$IMAGE" \
  python -m unittest discover -s tests -p 'test*.py'

echo "[2/2] R tests (minimal)" >&2
docker run --rm -v "$PWD":/work -w /work "$IMAGE" \
  Rscript r/tests/run_tests.R
