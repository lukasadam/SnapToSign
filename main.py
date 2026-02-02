import argparse
import scanpy as sc
from pathlib import Path
from typing import Dict, Optional
from converter import export_h5ad_to_signac_dir

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Convert h5ad files to Signac directory structure.")
    parser.add_argument("--rna_h5ad", type=str, required=True, help="Path to RNA h5ad file.")
    parser.add_argument("--atac_h5ad", type=str, required=True, help="Path to ATAC h5ad file.")
    parser.add_argument("--fragment_files", type=str, nargs='+', required=True, help="List of fragment file paths.")
    parser.add_argument("--out_dir", type=str, required=True, help="Output directory for Signac files.")
    parser.add_argument("--rna_layer", type=str, default=None, help="RNA layer to use.")
    parser.add_argument("--atac_layer", type=str, default=None, help="ATAC layer to use.")
    parser.add_argument("--copy_fragments", action='store_true', help="Whether to copy fragment files instead of symlinking.")
    args = parser.parse_args()
    
    # Load h5ad files
    adata_rna = sc.read_h5ad(args.rna_h5ad)
    adata_atac = sc.read_h5ad(args.atac_h5ad)
    fragment_files = [Path(f) for f in args.fragment_files]
    out_dir = Path(args.out_dir)
    rna_layer = args.rna_layer
    atac_layer = args.atac_layer
    copy_fragments = args.copy_fragments

    # Convert h5ad to Signac directory structure
    export_h5ad_to_signac_dir(
        adata_rna=adata_rna,
        adata_atac=adata_atac,
        fragment_files=fragment_files,
        out_dir=out_dir,
        rna_layer=rna_layer,
        atac_layer=atac_layer,
        copy_fragments=copy_fragments
    )

if __name__ == "__main__":
    main()
