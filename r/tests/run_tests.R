#!/usr/bin/env Rscript

# Minimal R tests (not a package): exit non-zero on failure.

gzip_file <- function(src, dst) {
    in_con <- file(src, open = "rb")
    on.exit(close(in_con), add = TRUE)
    out_con <- gzfile(dst, open = "wb")
    on.exit(close(out_con), add = TRUE)

    repeat {
        buf <- readBin(in_con, what = "raw", n = 1024 * 1024)
        if (length(buf) == 0) {
            break
        }
        writeBin(buf, out_con)
    }

    unlink(src)
}

write_tsv_gz <- function(df, path) {
    tmp <- sub("\\.gz$", "", path)
    write.table(
        df,
        file = tmp,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
    gzip_file(tmp, path)
}

write_csv_gz <- function(df, path) {
    tmp <- sub("\\.gz$", "", path)
    write.csv(df, file = tmp)
    gzip_file(tmp, path)
}

ca <- commandArgs(trailingOnly = FALSE)
file_arg <- ca[grepl("^--file=", ca)]
script_path <- if (length(file_arg) > 0) {
    sub("^--file=", "", file_arg[1])
} else {
    NA_character_
}
script_dir <- if (!is.na(script_path) && nzchar(script_path)) {
    dirname(normalizePath(script_path, mustWork = TRUE))
} else {
    getwd()
}

source(file.path(script_dir, "..", "converter.R"))

# Build a minimal exported structure
base <- tempfile("snaptosign_disc_")
dir.create(base)
dir.create(file.path(base, "rna"), recursive = TRUE)
dir.create(file.path(base, "atac"), recursive = TRUE)
dir.create(file.path(base, "meta"), recursive = TRUE)
dir.create(file.path(base, "reductions"), recursive = TRUE)
dir.create(file.path(base, "fragments"), recursive = TRUE)

barcodes <- data.frame(V1 = c("c1", "c2"))
write_tsv_gz(barcodes, file.path(base, "rna", "barcodes.tsv.gz"))
write_tsv_gz(barcodes, file.path(base, "atac", "barcodes.tsv.gz"))

rna_features <- data.frame(
    V1 = c("g1", "g2"),
    V2 = c("g1", "g2"),
    V3 = c("Gene Expression", "Gene Expression")
)
write_tsv_gz(rna_features, file.path(base, "rna", "features.tsv.gz"))

peaks <- c("chr1:1-2", "chr2:5-6")
atac_features <- data.frame(V1 = peaks, V2 = peaks, V3 = c("Peaks", "Peaks"))
write_tsv_gz(atac_features, file.path(base, "atac", "features.tsv.gz"))

# 10x convention: features x cells
suppressPackageStartupMessages({
    library(Matrix)
})

rna_mat <- Matrix(c(1, 0, 2, 3), nrow = 2, sparse = TRUE)
atac_mat <- Matrix(c(0, 1, 1, 0), nrow = 2, sparse = TRUE)

rna_tmp <- file.path(base, "rna", "matrix.mtx")
Matrix::writeMM(rna_mat, rna_tmp)
gzip_file(rna_tmp, file.path(base, "rna", "matrix.mtx.gz"))

atac_tmp <- file.path(base, "atac", "matrix.mtx")
Matrix::writeMM(atac_mat, atac_tmp)
gzip_file(atac_tmp, file.path(base, "atac", "matrix.mtx.gz"))

bed <- data.frame(
    chr = c("chr1", "chr2"),
    start = c(1, 5),
    end = c(2, 6),
    name = peaks
)
write_tsv_gz(bed, file.path(base, "atac", "peaks.bed.gz"))

meta <- data.frame(celltype = c("a", "b"), row.names = c("c1", "c2"))
write_csv_gz(meta, file.path(base, "meta", "metadata.csv.gz"))

obj <- import_snaptosign_dir(base)

stopifnot(inherits(obj, "Seurat"))
stopifnot(all(c("RNA", "ATAC") %in% names(obj@assays)))
stopifnot(ncol(obj) == 2)

cat("OK: R minimal tests passed\n")
