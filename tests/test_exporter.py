import gzip
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

from py.converter import export_h5ad_to_signac_dir


def _read_tsv_gz(path: Path, sep: str = "\t") -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        return pd.read_csv(f, sep=sep, header=None)


def _read_mtx_gz(path: Path):
    with gzip.open(path, "rb") as f:
        return mmread(f)


class TestExporter(unittest.TestCase):
    def test_export_writes_expected_structure(self) -> None:
        cells = ["c1", "c2", "c3"]

        rna_x = sp.csr_matrix(np.array([[1, 0], [0, 2], [3, 4]], dtype=np.int64))
        rna = ad.AnnData(
            X=rna_x,
            obs=pd.DataFrame(index=cells),
            var=pd.DataFrame(index=["g1", "g2"]),
        )

        atac_x = sp.csr_matrix(np.array([[0, 1], [1, 0], [2, 3]], dtype=np.int64))
        peaks = ["chr1:1-10", "chr2:5-15"]
        atac = ad.AnnData(
            X=atac_x,
            obs=pd.DataFrame(index=cells),
            var=pd.DataFrame(index=peaks),
        )

        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "disc"
            export_h5ad_to_signac_dir(rna=rna, atac=atac, out_dir=out_dir)

            expected = [
                out_dir / "rna" / "matrix.mtx.gz",
                out_dir / "rna" / "barcodes.tsv.gz",
                out_dir / "rna" / "features.tsv.gz",
                out_dir / "atac" / "matrix.mtx.gz",
                out_dir / "atac" / "barcodes.tsv.gz",
                out_dir / "atac" / "features.tsv.gz",
                out_dir / "atac" / "peaks.bed.gz",
                out_dir / "meta" / "metadata.csv.gz",
            ]
            for path in expected:
                self.assertTrue(path.exists(), f"Missing expected file: {path}")

            barcodes = _read_tsv_gz(out_dir / "rna" / "barcodes.tsv.gz")
            self.assertEqual(barcodes[0].tolist(), cells)

            rna_features = _read_tsv_gz(out_dir / "rna" / "features.tsv.gz")
            self.assertEqual(rna_features.shape[1], 3)
            self.assertEqual(rna_features[0].tolist(), ["g1", "g2"])

            atac_features = _read_tsv_gz(out_dir / "atac" / "features.tsv.gz")
            self.assertEqual(atac_features[0].tolist(), peaks)

            bed = _read_tsv_gz(out_dir / "atac" / "peaks.bed.gz")
            self.assertEqual(bed.shape[1], 4)
            self.assertEqual(bed.iloc[0, 0], "chr1")

            rna_mtx = _read_mtx_gz(out_dir / "rna" / "matrix.mtx.gz")
            self.assertEqual(rna_mtx.shape, (2, 3))  # features x cells (10x convention)

    def test_export_accepts_fragment_list(self) -> None:
        cells = ["c1", "c2"]

        rna = ad.AnnData(
            X=sp.csr_matrix(np.array([[1], [2]], dtype=np.int64)),
            obs=pd.DataFrame(index=cells),
            var=pd.DataFrame(index=["g1"]),
        )

        atac = ad.AnnData(
            X=sp.csr_matrix(np.array([[1], [0]], dtype=np.int64)),
            obs=pd.DataFrame(index=cells),
            var=pd.DataFrame(index=["chr1:1-2"]),
        )

        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            frag = tmp_path / "fragA.tsv.gz"
            frag.write_bytes(b"")
            (tmp_path / "fragA.tsv.gz.tbi").write_bytes(b"")

            out_dir = tmp_path / "disc"
            export_h5ad_to_signac_dir(
                rna=rna,
                atac=atac,
                out_dir=out_dir,
                fragment_files=[frag],
            )

            self.assertTrue((out_dir / "fragments" / "fragA.tsv.gz").exists())
            self.assertTrue((out_dir / "fragments" / "fragA.tsv.tbi").exists())
