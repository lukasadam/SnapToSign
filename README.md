# SnapToSign (h5ad → Signac Converter)

A small command-line tool to convert paired RNA and ATAC .h5ad files into a Signac-compatible directory structure, including fragment files, for downstream analysis in R / Seurat–Signac.

## ✨ What this does
•	Loads RNA and ATAC AnnData (.h5ad) objects
•	Exports them into a directory layout expected by Signac
•	Links or copies ATAC fragment files
•	Supports optional layers for RNA and ATAC assays

## 🔧 Installation (using uv)

This project uses uv for fast, reproducible Python dependency management.

1️⃣ Install uv

If you don’t have uv yet:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2️⃣ Clone the repository

```bash
git clone https://github.com/lukasadam/SnapToSign
cd SnapToSign
```

3️⃣ Create a virtual environment

```bash
uv create venv
```

4️⃣ Install dependencies for the project

```bash
uv pip install -e .
```

## 🚀 Usage

Run the script from the command line:

```bash
py2disc --rna_h5ad path/to/rna.h5ad --atac_h5ad path/to/atac.h5ad --fragment_files path/to/fragments/* --out_dir path/to/output
```