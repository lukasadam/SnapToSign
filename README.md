# SnapToSign

[![tests](https://github.com/lukasadam/SnapToSign/actions/workflows/tests.yml/badge.svg)](https://github.com/lukasadam/SnapToSign/actions/workflows/tests.yml)
[![ruff](https://github.com/lukasadam/SnapToSign/actions/workflows/ruff.yml/badge.svg)](https://github.com/lukasadam/SnapToSign/actions/workflows/ruff.yml)
[![air](https://github.com/lukasadam/SnapToSign/actions/workflows/air.yml/badge.svg)](https://github.com/lukasadam/SnapToSign/actions/workflows/air.yml)
[![docker](https://github.com/lukasadam/SnapToSign/actions/workflows/docker.yml/badge.svg)](https://github.com/lukasadam/SnapToSign/actions/workflows/docker.yml)
[![codecov](https://codecov.io/gh/lukasadam/SnapToSign/branch/main/graph/badge.svg)](https://codecov.io/gh/lukasadam/SnapToSign)

Convert paired RNA+ATAC AnnData (.h5ad) to a Signac/Seurat-ready format:

1) `py2disc` exports a Signac-compatible directory structure
2) `disc2r` loads that directory and writes a Seurat object (`.rds`)

## Quickstart (Apptainer/Singularity)

This repo publishes a container image to GHCR. On many HPC clusters you cannot run Docker/Podman directly, but you can run the same image via Apptainer (or Singularity).

### Option A: pull a local `.sif` (recommended)

```bash
apptainer pull snaptosign.sif docker://ghcr.io/lukasadam/snaptosign:latest
```

If `snaptosign.sif` is in your current directory, `run_snaptosign.sh` will use it by default.

(If your cluster uses `singularity`, the commands are the same: replace `apptainer` with `singularity`.)

### Option B: run directly from GHCR

```bash
SNAPTOSIGN_IMAGE=docker://ghcr.io/lukasadam/snaptosign:latest \
	./run_snaptosign.sh data/rna.h5ad data/atac.h5ad out data/object.rds
```

Run everything in one step:

```bash
chmod +x run_snaptosign.sh
SNAPTOSIGN_IMAGE=$PWD/snaptosign.sif \
	./run_snaptosign.sh data/rna.h5ad data/atac.h5ad out data/object.rds
```

Absolute paths work too (useful on clusters):

```bash
SNAPTOSIGN_IMAGE=/path/to/snaptosign.sif \
	./run_snaptosign.sh \
		/abs/path/rna.h5ad \
		/abs/path/atac.h5ad \
		/abs/path/conversion_out \
		/abs/path/object.rds
```

Equivalent (explicit flags):

```bash
./run_snaptosign.sh --rna_h5ad data/rna.h5ad --atac_h5ad data/atac.h5ad --disc_out_dir out --out_rds data/object.rds
```

Optional fragments (globs work):

```bash
./run_snaptosign.sh \
	data/rna.h5ad \
	data/atac.h5ad \
	out \
	data/object.rds \
	data/fragments/*.tsv.gz
```

Or with flags:

```bash
./run_snaptosign.sh \
  --rna_h5ad data/rna.h5ad \
  --atac_h5ad data/atac.h5ad \
  --disc_out_dir out \
  --out_rds data/object.rds \
  --fragment_files data/fragments/*.tsv.gz
```

If you tagged the image differently:

```bash
SNAPTOSIGN_IMAGE=docker://ghcr.io/lukasadam/snaptosign:sha-<commit> ./run_snaptosign.sh \
	data/rna.h5ad \
	data/atac.h5ad \
	out \
	data/object.rds
```

### Build locally with Docker (optional, for development)

```bash
docker build -t snaptosign:local -f docker/snaptosign.Dockerfile .
echo "If you want to run a locally-built Docker image, use Docker directly (or push it to a registry and use the Apptainer workflow above)."
```


## CLI usage

Export directory structure:

```bash
py2disc --rna_h5ad rna.h5ad --atac_h5ad atac.h5ad --out_dir disc_out
```

Import to `.rds`:

```bash
disc2r --disc_dir disc_out --out_file object.rds
```

## Local (no Docker)

Python (installs `py2disc` and `snaptosign`):

```bash
uv pip install -e .
```

R: use `r/disc2r.R` (requires Seurat + Signac + optparse + GenomicRanges).