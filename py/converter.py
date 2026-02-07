"""Converter functions to export AnnData objects into Seurat/Signac-compatible directory structure."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.io import mmwrite


def _ensure_csr(matrix: sp.spmatrix) -> sp.csr_matrix:
    """Ensure the input sparse matrix is in CSR format.

    Parameters
    ----------
    matrix : sp.spmatrix
        Input sparse matrix.

    Returns
    -------
    sp.csr_matrix
        Sparse matrix in CSR format.
    """
    if sp.isspmatrix_csr(matrix):
        return matrix
    if sp.issparse(matrix):
        return matrix.tocsr()
    raise ValueError("Input matrix must be a sparse matrix.")


def _gzip_file(src: Path, dst: Path) -> None:
    """Gzip the file at src and save it to dst, then remove the original file.

    Parameters
    ----------
    src : Path
        Path to the source file.
    dst : Path
        Path to the destination file.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(src, "rb") as f_in:
        with gzip.open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    src.unlink()


def _write_tsv_gz(df_or_series, path_gz: Path, header: bool = False) -> None:
    """Write a DataFrame or Series to a gzipped TSV file.

    Parameters
    ----------
    df_or_series : pd.DataFrame or pd.Series
        The DataFrame or Series to write.
    path_gz : Path
        The path to the output gzipped TSV file.
    header : bool, optional
        Whether to write the column names (default is False).
    """
    path_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp = path_gz.with_suffix("")  # remove .gz for temp
    if isinstance(df_or_series, pd.Series):
        df_or_series.to_csv(tmp, sep="\t", index=False, header=header)
    else:
        df_or_series.to_csv(tmp, sep="\t", index=False, header=header)
    _gzip_file(tmp, path_gz)


def _write_csv_gz(df: pd.DataFrame, path_gz: Path) -> None:
    """Write a DataFrame to a gzipped CSV file.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to write.
    path_gz : Path
        The path to the output gzipped CSV file.
    """
    path_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp = path_gz.with_suffix("")
    df.to_csv(tmp)
    _gzip_file(tmp, path_gz)


def _write_mtx_gz(x_csr: sp.csr_matrix, path_gz: Path) -> None:
    """Write a CSR sparse matrix to a gzipped Matrix Market (.mtx) file.

    Parameters
    ----------
    x_csr : sp.csr_matrix
        The input CSR sparse matrix.
    path_gz : Path
        The path to the output gzipped Matrix Market file.
    """
    path_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp = path_gz.with_suffix("")  # .mtx temp

    x_out = x_csr.T.astype(np.float64)  # ensure float64 for 10x compatibility
    mmwrite(tmp, x_out)
    _gzip_file(tmp, path_gz)


def _export_obsm(adata, key: str, filename: str, out_dir: Path) -> None:
    """Export an obsm matrix to a gzipped CSV file.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object containing the obsm matrix.
    key : str
        The key in adata.obsm to export.
    filename : str
        The filename for the output gzipped CSV file (without path).
    out_dir : Path
        The output directory path.
    """
    if key in adata.obsm_keys():
        arr = np.asarray(adata.obsm[key])
        df = pd.DataFrame(arr, index=adata.obs_names.astype(str))
        _write_csv_gz(df, out_dir / "reductions" / filename)


def _peak_to_bed(features: pd.Index) -> pd.DataFrame:
    """Convert peak features to BED format DataFrame.

    Parameters
    ----------
    features : pd.Index
        The peak features to convert.

    Returns
    -------
    pd.DataFrame
        A DataFrame in BED format.
    """
    # robust split
    chromosome = features.to_series().str.split(":", n=1, expand=True)[0].astype(str)
    rest = features.to_series().str.split(":", n=1, expand=True)[1].astype(str)
    start = rest.str.split("-", n=1, expand=True)[0].astype(int)
    end = rest.str.split("-", n=1, expand=True)[1].astype(int)
    bed = pd.DataFrame(
        {"chr": chromosome.values, "start": start.values, "end": end.values}
    )
    bed["name"] = features.astype(str).values
    return bed


def export_h5ad_to_signac_dir(
    rna: ad.AnnData,
    atac: ad.AnnData,
    out_dir: Path,
    fragment_files: Optional[Dict[str, Path]] = None,
    rna_layer: Optional[str] = None,
    atac_layer: Optional[str] = None,
    copy_fragments: bool = True,
) -> None:
    """
    Exports RNA and ATAC AnnData objects into a fixed directory structure
    usable by Seurat/Signac in R.

    Parameters
    ----------
    rna : ad.AnnData
        AnnData object containing RNA data.
    atac : ad.AnnData
        AnnData object containing ATAC data.
    out_dir : Path
        Output directory path.
    fragment_files : Optional[Dict[str, Path]], optional
        Dictionary mapping sample IDs to fragment file paths (default is None).
    rna_layer : Optional[str], optional
        Layer name in rna to use for RNA counts (default is None, uses .X).
    atac_layer : Optional[str], optional
        Layer name in atac to use for ATAC counts (default is None, uses .X).
    copy_fragments : bool, optional
        Whether to copy fragment files instead of symlinking (default is True).
    """
    out_dir = Path(out_dir)
    (out_dir / "rna").mkdir(parents=True, exist_ok=True)
    (out_dir / "atac").mkdir(parents=True, exist_ok=True)
    (out_dir / "meta").mkdir(parents=True, exist_ok=True)
    (out_dir / "reductions").mkdir(parents=True, exist_ok=True)
    (out_dir / "fragments").mkdir(parents=True, exist_ok=True)

    # Align cells between modalities
    common = rna.obs_names.intersection(atac.obs_names)
    if len(common) == 0:
        raise ValueError("No overlapping cells between RNA and ATAC AnnData objects.")

    # Subset to common cells
    rna = rna[common].copy()
    atac = atac[common].copy()

    # Extract matrices
    x_rna = rna.layers[rna_layer] if (rna_layer is not None) else rna.X
    x_atac = atac.layers[atac_layer] if (atac_layer is not None) else atac.X
    x_rna = _ensure_csr(x_rna)
    x_atac = _ensure_csr(x_atac)

    # Export barcodes (cell IDs)
    barcodes = pd.Series(common.astype(str), name="barcode")
    _write_tsv_gz(barcodes, out_dir / "rna" / "barcodes.tsv.gz", header=False)
    _write_tsv_gz(barcodes, out_dir / "atac" / "barcodes.tsv.gz", header=False)

    # Export features
    # 10x features.tsv has 3 columns: id, name, type
    # We'll try to use var fields if present.
    var_rna = rna.var.copy()
    gene_id = (
        var_rna["gene_id"].astype(str)
        if "gene_id" in var_rna.columns
        else var_rna.index.astype(str)
    )
    gene_name = (
        var_rna["gene_name"].astype(str)
        if "gene_name" in var_rna.columns
        else var_rna.index.astype(str)
    )
    rna_features = pd.DataFrame(
        {"id": gene_id.values, "name": gene_name.values, "type": "Gene Expression"}
    )
    _write_tsv_gz(rna_features, out_dir / "rna" / "features.tsv.gz", header=False)

    # Parse ATAC peaks from var_names, which are expected to be in "chr:start-end" format
    peaks = atac.var_names.astype(str)
    atac_features = pd.DataFrame({"id": peaks, "name": peaks, "type": "Peaks"})
    _write_tsv_gz(atac_features, out_dir / "atac" / "features.tsv.gz", header=False)

    # Export ATAC peaks in BED format (chr, start, end, name)
    bed = _peak_to_bed(atac.var_names)
    _write_tsv_gz(bed, out_dir / "atac" / "peaks.bed.gz", header=False)

    # Export count matrices in Matrix Market format (cells x features, transposed for 10x compatibility)
    _write_mtx_gz(x_rna, out_dir / "rna" / "matrix.mtx.gz")
    _write_mtx_gz(x_atac, out_dir / "atac" / "matrix.mtx.gz")

    # Export metadata (obs)
    # Save once (shared barcodes index)
    meta = rna.obs.copy()
    # Ensure rownames in R match barcodes
    meta.index = common.astype(str)
    _write_csv_gz(meta, out_dir / "meta" / "metadata.csv.gz")

    # Export dimensionality reductions (obsm)
    _export_obsm(rna, "X_pca", "rna_pca.csv.gz", out_dir=out_dir)
    _export_obsm(rna, "X_umap", "rna_umap.csv.gz", out_dir=out_dir)
    _export_obsm(atac, "X_spectral", "atac_spectral.csv.gz", out_dir=out_dir)
    _export_obsm(atac, "X_umap", "atac_umap.csv.gz", out_dir=out_dir)

    # If fragment_files is provided, copy/symlink them into the fragments directory
    if fragment_files is not None:
        # Copy/symlink all provided fragments into out_dir/fragments
        for sample_id, fpath in fragment_files.items():
            fpath = Path(fpath)
            if not fpath.exists():
                continue
            dest = out_dir / "fragments" / f"{sample_id}.tsv.gz"
            if dest.exists():
                continue
            if copy_fragments:
                shutil.copy2(fpath, dest)
            else:
                dest.symlink_to(fpath)

            # Also copy the index file if present
            fpath_idx = fpath.with_suffix(".tsv.gz.tbi")
            if fpath_idx.exists():
                dest_idx = dest.with_suffix(".tbi")
                if copy_fragments:
                    shutil.copy2(fpath_idx, dest_idx)
                else:
                    dest_idx.symlink_to(fpath_idx)

    print(f"Export complete: {out_dir}")
