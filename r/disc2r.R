#!/usr/bin/env Rscript

# Command-line interface for loading SnapToSign directory structure into Seurat/Signac
#
# Usage example:
#   ./r/disc2r.R --input_dir exported_data --output_rds seurat_object.rds --genome hg38

suppressPackageStartupMessages({
  library(optparse)
})

.get_script_dir <- function() {
  # Works for Rscript execution (commandArgs contains --file=...)
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- ca[grepl("^--file=", ca)]
  if (length(file_arg) == 0) {
    return(getwd())
  }
  script_path <- sub("^--file=", "", file_arg[1])
  return(dirname(normalizePath(script_path, mustWork = TRUE)))
}

script_dir <- .get_script_dir()


parse_cli_args <- function() {
  option_list <- list(
    make_option(
      c("--disc_dir"),
      type = "character",
      default = NULL,
      dest = "disc_dir",
      help = "Directory containing exported data structure (output of py2disc)",
      metavar = "PATH"
    ),
    make_option(
      c("--out_file"),
      type = "character",
      default = NULL,
      dest = "out_file",
      help = "Output .rds file path for the Seurat object",
      metavar = "FILE"
    ),
    # Backwards-compatible aliases
    make_option(
      c("--input_dir"),
      type = "character",
      default = NULL,
      dest = "disc_dir",
      help = "(deprecated) use --disc_dir",
      metavar = "PATH"
    ),
    make_option(
      c("--output_rds"),
      type = "character",
      default = NULL,
      dest = "out_file",
      help = "(deprecated) use --out_file",
      metavar = "FILE"
    ),
    make_option(
      c("--fragment_dir"),
      type = "character",
      default = NULL,
      help = "Directory containing fragment files (default: input_dir/fragments)",
      metavar = "PATH"
    )
  )

  opt_parser <- OptionParser(
    usage = "Usage: disc2r [options]",
    option_list = option_list,
    description = "\nLoad SnapToSign directory structure and create Seurat/Signac object"
  )

  opt <- parse_args(opt_parser)

  if (is.null(opt$disc_dir)) {
    print_help(opt_parser)
    stop("--disc_dir is required", call. = FALSE)
  }
  if (is.null(opt$out_file)) {
    print_help(opt_parser)
    stop("--out_file is required", call. = FALSE)
  }

  opt
}


main <- function() {
  opt <- parse_cli_args()

  # Only load the heavy Seurat/Signac converter after args are parsed.
  # This keeps `disc2r.R --help` usable even if those packages aren't installed.
  source(file.path(script_dir, "converter.R"))

  message("disc2r (SnapToSign)")
  message("=================")
  message(sprintf("disc_dir: %s", opt$disc_dir))
  message(sprintf("out_file: %s", opt$out_file))
  if (!is.null(opt$fragment_dir)) {
    message(sprintf("Fragment directory: %s", opt$fragment_dir))
  }
  message("")

  seurat_obj <- import_snaptosign_dir(
    input_dir = opt$disc_dir,
    fragment_dir = opt$fragment_dir
  )

  save_seurat_rds(seurat_obj, opt$out_file)

  message("\n======================")
  message("Conversion complete!")
  message(sprintf("Seurat object saved to: %s", opt$out_file))
  message(sprintf("Number of cells: %d", ncol(seurat_obj)))
  message(sprintf(
    "Assays: %s",
    paste(names(seurat_obj@assays), collapse = ", ")
  ))
  message(sprintf(
    "Reductions: %s",
    paste(names(seurat_obj@reductions), collapse = ", ")
  ))
}


if (!interactive()) {
  main()
}
