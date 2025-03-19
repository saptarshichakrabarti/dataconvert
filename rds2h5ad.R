#!/usr/bin/env Rscript
# rds2h5ad.R: Convert an RDS file (Seurat or monocle3 cell_data_set) to an AnnData (H5AD) file.
# Usage: Rscript master_converter.R <input_file.RDS>

# --- Step 1: Ensure BiocManager is installed and install Bioconductor version 3.20 ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(version = "3.20", ask = FALSE)

# --- Step 2: Install Bioconductor dependencies ---
bioc_dependencies <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr')
for(pkg in bioc_dependencies) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# --- Step 3: Install monocle3 from GitHub using devtools ---
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("monocle3", quietly = TRUE)) {
  devtools::install_github('cole-trapnell-lab/monocle3')
}

# --- Step 4: Install additional CRAN packages ---
cran_packages <- c("shiny", "ggplot2", "dplyr", "Matrix", "irlba")
for(pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# --- Step 5: Install conversion packages: SeuratDisk and zellkonverter ---
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  devtools::install_github("mojaveazure/seurat-disk")
}
if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  BiocManager::install("zellkonverter", ask = FALSE)
}

# --- Load Required Libraries ---
suppressPackageStartupMessages({
  library(monocle3)
  library(Matrix)
  library(irlba)
  library(methods)  # For handling S4 objects
})

# --- Command-line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript master_converter.R <input_file.RDS>")
}
input_file <- args[1]

# Generate output filename by replacing .RDS with .h5ad
output_file <- sub("\\.RDS$", ".h5ad", sub("\\.rds$", ".h5ad", input_file))

# --- Read the RDS File ---
object <- readRDS(input_file)

# --- Conversion Based on Object Type ---
if (inherits(object, "Seurat")) {
  suppressPackageStartupMessages({
    library(SeuratDisk)
  })
  cat("Detected a Seurat object. Converting to intermediary H5Seurat format...\n")

  temp_h5seurat <- tempfile(fileext = ".h5Seurat")
  SaveH5Seurat(object, filename = temp_h5seurat, overwrite = TRUE)

  cat("Converting H5Seurat to H5AD...\n")
  Convert(temp_h5seurat, dest = "h5ad", overwrite = TRUE)

  temp_h5ad <- sub(".h5Seurat$", ".h5ad", temp_h5seurat)
  file.rename(temp_h5ad, output_file)
  cat("Conversion complete! Output saved as:", output_file, "\n")

} else if (inherits(object, "cell_data_set")) {
  suppressPackageStartupMessages({
    library(zellkonverter)
  })
  cat("Detected a monocle3 cell_data_set object. Converting using zellkonverter...\n")

  writeH5AD(object, file = output_file)
  cat("Conversion complete! Output saved as:", output_file, "\n")

} else {
  stop("Unsupported object type. This script supports Seurat and monocle3 cell_data_set objects.")
}