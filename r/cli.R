#!/usr/bin/env Rscript

# Compatibility wrapper (kept for existing docs/usages).
# The preferred CLI entrypoint is: r/disc2r.R

ca <- commandArgs(trailingOnly = FALSE)
file_arg <- ca[grepl("^--file=", ca)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()

source(file.path(script_dir, "disc2r.R"))
