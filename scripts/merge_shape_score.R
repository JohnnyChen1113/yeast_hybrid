#!/usr/bin/env Rscript

# Suppress all warnings and messages during package loading
suppressMessages({
  suppressWarnings({
    library(tidyverse)
    library(optparse)
  })
})

# Version number
VERSION <- "0.2.2"

# Define command line options
option_list <- list(
  make_option(c("-a", "--assigned"), type="character", default=NULL,
              help="Path to assignedClusters file", metavar="character"),
  make_option(c("-s", "--shape"), type="character", default=NULL,
              help="Path to shape file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Path to output file [default: assignedClusters prefix with _merged.txt]"),
  make_option(c("-d", "--dominate"), action="store_true", default=FALSE,
              help="Output only the dominant cluster for each gene based on maximum tags value"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE,
              help="Print version number and exit")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, print_help_and_exit=TRUE)  # Automatically handle -h/--help

# Print version if requested
if (opt$version) {
  cat("merge_shape_score.R version", VERSION, "\n")
  quit(status=0)
}

# If no arguments provided (excluding help/version), print help and exit
if (length(commandArgs(trailingOnly=TRUE)) == 0) {
  print_help(opt_parser)
  quit(status=0)
}

# Check if required inputs are provided
if (is.null(opt$assigned) || is.null(opt$shape)) {
  stop("Both --assigned and --shape arguments are required")
}

# Set default output file name if not provided
if (is.null(opt$output)) {
  prefix <- sub("\\.txt$", "", basename(opt$assigned))
  opt$output <- paste0(prefix, "_merged.txt")
}

# Read the input files with show_col_types = FALSE
table1 <- read_tsv(opt$assigned, show_col_types = FALSE)
table2 <- read_tsv(opt$shape, show_col_types = FALSE)

# Merge tables based on cluster column
merged_table <- table1 %>%
  left_join(select(table2, cluster, shape.score), by = "cluster")

# If dominate option is used, filter for dominant cluster per gene
if (opt$dominate) {
  merged_table <- merged_table %>%
    group_by(gene) %>%
    slice(which.max(tags)) %>%
    ungroup()
}

# Write the output
write_tsv(merged_table, opt$output)
