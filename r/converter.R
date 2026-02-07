#' Load SnapToSign Directory Structure into Seurat/Signac Object
#'
#' This script loads RNA and ATAC data exported from Python AnnData objects
#' into a Seurat object with Signac for multimodal scRNA-seq and scATAC-seq analysis.

library(Seurat)
library(Signac)
library(Matrix)


#' Read 10x-style Matrix Market files
#'
#' @param dir_path Path to directory containing matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
#' @return A sparse matrix
#' @keywords internal
.read_10x_mtx <- function(dir_path) {
  mtx_path <- file.path(dir_path, "matrix.mtx.gz")
  barcodes_path <- file.path(dir_path, "barcodes.tsv.gz")
  features_path <- file.path(dir_path, "features.tsv.gz")
  
  # Read matrix
  mat <- Matrix::readMM(gzfile(mtx_path))
  
  # Read barcodes (cell names)
  barcodes <- read.table(gzfile(barcodes_path), header = FALSE, stringsAsFactors = FALSE)
  colnames(mat) <- barcodes$V1
  
  # Read features (gene/peak names)
  features <- read.table(gzfile(features_path), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  rownames(mat) <- features$V1
  
  return(mat)
}


#' Read metadata CSV
#'
#' @param meta_dir Path to meta directory
#' @return A data frame with metadata
#' @keywords internal
.read_metadata <- function(meta_dir) {
  meta_path <- file.path(meta_dir, "metadata.csv.gz")
  if (!file.exists(meta_path)) {
    return(NULL)
  }
  meta <- read.csv(gzfile(meta_path), row.names = 1, stringsAsFactors = FALSE)
  return(meta)
}


#' Read dimensionality reduction
#'
#' @param reductions_dir Path to reductions directory
#' @param filename Name of the reduction file (e.g., "rna_pca.csv.gz")
#' @return A matrix with reduction coordinates
#' @keywords internal
.read_reduction <- function(reductions_dir, filename) {
  red_path <- file.path(reductions_dir, filename)
  if (!file.exists(red_path)) {
    return(NULL)
  }
  red <- read.csv(gzfile(red_path), row.names = 1, stringsAsFactors = FALSE)
  return(as.matrix(red))
}


#' Import SnapToSign directory structure into Seurat/Signac object
#'
#' @param input_dir Path to the directory containing exported data
#' @param fragment_dir Optional path to fragment files directory (default: input_dir/fragments)
#' @return A Seurat object with RNA and ATAC assays
#' @export
import_snaptosign_dir <- function(input_dir, fragment_dir = NULL) {
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  
  if (is.null(fragment_dir)) {
    fragment_dir <- file.path(input_dir, "fragments")
  }
  
  message("Loading RNA data...")
  rna_dir <- file.path(input_dir, "rna")
  rna_mat <- .read_10x_mtx(rna_dir)
  
  message("Loading ATAC data...")
  atac_dir <- file.path(input_dir, "atac")
  atac_mat <- .read_10x_mtx(atac_dir)
  
  # Read peaks BED file for ATAC
  peaks_path <- file.path(atac_dir, "peaks.bed.gz")
  peaks <- read.table(gzfile(peaks_path), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  colnames(peaks) <- c("chr", "start", "end", "name")
  
  # Create GRanges object for peaks
  peak_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    peaks,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end"
  )
  names(peak_ranges) <- peaks$name
  
  message("Creating Seurat object with RNA assay...")
  # Create Seurat object with RNA data
  seurat_obj <- CreateSeuratObject(
    counts = rna_mat,
    assay = "RNA",
    project = "SnapToSign"
  )
  
  message("Adding ATAC assay...")
  # Create ChromatinAssay for ATAC data
  # Find fragment files
  fragment_files <- NULL
  if (dir.exists(fragment_dir)) {
    fragment_files <- list.files(fragment_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
    if (length(fragment_files) > 0) {
      message(sprintf("Found %d fragment file(s)", length(fragment_files)))
    } else {
      fragment_files <- NULL
    }
  }
  
  atac_assay <- CreateChromatinAssay(
    counts = atac_mat,
    ranges = peak_ranges,
    fragments = fragment_files
  )
  
  seurat_obj[["ATAC"]] <- atac_assay
  
  # Load metadata
  message("Loading metadata...")
  meta_dir <- file.path(input_dir, "meta")
  metadata <- .read_metadata(meta_dir)
  if (!is.null(metadata)) {
    # Add metadata to Seurat object
    # Ensure metadata rows match Seurat object cells
    shared_cells <- intersect(colnames(seurat_obj), rownames(metadata))
    if (length(shared_cells) > 0) {
      seurat_obj <- AddMetaData(seurat_obj, metadata[shared_cells, , drop = FALSE])
    }
  }
  
  # Load dimensionality reductions
  message("Loading dimensionality reductions...")
  reductions_dir <- file.path(input_dir, "reductions")
  
  # RNA PCA
  rna_pca <- .read_reduction(reductions_dir, "rna_pca.csv.gz")
  if (!is.null(rna_pca)) {
    shared_cells <- intersect(colnames(seurat_obj), rownames(rna_pca))
    if (length(shared_cells) > 0) {
      rna_pca <- rna_pca[shared_cells, , drop = FALSE]
      colnames(rna_pca) <- paste0("PC_", seq_len(ncol(rna_pca)))
      seurat_obj[["pca"]] <- CreateDimReducObject(
        embeddings = rna_pca,
        key = "PC_",
        assay = "RNA"
      )
    }
  }
  
  # RNA UMAP
  rna_umap <- .read_reduction(reductions_dir, "rna_umap.csv.gz")
  if (!is.null(rna_umap)) {
    shared_cells <- intersect(colnames(seurat_obj), rownames(rna_umap))
    if (length(shared_cells) > 0) {
      rna_umap <- rna_umap[shared_cells, , drop = FALSE]
      colnames(rna_umap) <- paste0("UMAP_", seq_len(ncol(rna_umap)))
      seurat_obj[["umap"]] <- CreateDimReducObject(
        embeddings = rna_umap,
        key = "UMAP_",
        assay = "RNA"
      )
    }
  }
  
  # ATAC Spectral (LSI)
  atac_spectral <- .read_reduction(reductions_dir, "atac_spectral.csv.gz")
  if (!is.null(atac_spectral)) {
    shared_cells <- intersect(colnames(seurat_obj), rownames(atac_spectral))
    if (length(shared_cells) > 0) {
      atac_spectral <- atac_spectral[shared_cells, , drop = FALSE]
      colnames(atac_spectral) <- paste0("LSI_", seq_len(ncol(atac_spectral)))
      seurat_obj[["lsi"]] <- CreateDimReducObject(
        embeddings = atac_spectral,
        key = "LSI_",
        assay = "ATAC"
      )
    }
  }
  
  # ATAC UMAP
  atac_umap <- .read_reduction(reductions_dir, "atac_umap.csv.gz")
  if (!is.null(atac_umap)) {
    shared_cells <- intersect(colnames(seurat_obj), rownames(atac_umap))
    if (length(shared_cells) > 0) {
      atac_umap <- atac_umap[shared_cells, , drop = FALSE]
      colnames(atac_umap) <- paste0("atacUMAP_", seq_len(ncol(atac_umap)))
      seurat_obj[["umap.atac"]] <- CreateDimReducObject(
        embeddings = atac_umap,
        key = "atacUMAP_",
        assay = "ATAC"
      )
    }
  }
  
  message("Seurat/Signac object created successfully!")
  return(seurat_obj)
}


#' Save Seurat object to RDS file
#'
#' @param seurat_obj Seurat object to save
#' @param output_path Path to output RDS file
#' @export
save_seurat_rds <- function(seurat_obj, output_path) {
  message(sprintf("Saving Seurat object to: %s", output_path))
  saveRDS(seurat_obj, file = output_path)
  message("Save complete!")
}
