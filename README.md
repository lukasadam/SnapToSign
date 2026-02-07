# SnapToSign

Convert paired RNA+ATAC AnnData (.h5ad) to a Signac/Seurat-ready format:

1) `py2disc` exports a Signac-compatible directory structure
2) `disc2r` loads that directory and writes a Seurat object (`.rds`)

## Quickstart (Docker)

Build the image:

```bash
docker build -t snaptosign:local .
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
SNAPTOSIGN_IMAGE=snaptosign:ci ./run_snaptosign.sh \
	data/rna.h5ad \
	data/atac.h5ad \
	out \
	data/object.rds
```

## Apptainer/Singularity

Build a `.sif`:

```bash
apptainer build snaptosign.sif Apptainer.def
```

Run (no bind flags needed if your cluster auto-binds `$PWD`; otherwise add `-B $PWD`):

```bash
apptainer exec snaptosign.sif py2disc --help
apptainer exec snaptosign.sif disc2r --help
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