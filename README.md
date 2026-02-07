# SnapToSign

Convert paired RNA+ATAC AnnData (.h5ad) to a Signac/Seurat-ready format:

1) `py2disc` exports a Signac-compatible directory structure
2) `disc2r` loads that directory and writes a Seurat object (`.rds`)

## Quickstart (Docker)

Pull the latest image from GHCR:

```bash
docker pull ghcr.io/lukasadam/snaptosign:latest
```

Run everything in one step (mounted working directory):

```bash
chmod +x run_snaptosign.sh
./run_snaptosign.sh data/rna.h5ad data/atac.h5ad out data/object.rds
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

If you tagged the image differently:

```bash
SNAPTOSIGN_IMAGE=ghcr.io/lukasadam/snaptosign:sha-<commit> ./run_snaptosign.sh \
	data/rna.h5ad \
	data/atac.h5ad \
	out \
	data/object.rds
```

### Build locally (optional)

```bash
docker build -t snaptosign:local .
SNAPTOSIGN_IMAGE=snaptosign:local ./run_snaptosign.sh data/rna.h5ad data/atac.h5ad out data/object.rds
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