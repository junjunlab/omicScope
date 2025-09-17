# omicScope

<!-- badges: start -->

[![R-CMD-check](https://github.com/junjunlab/omicScope/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/junjunlab/omicScope/actions/workflows/R-CMD-check.yaml) [![Codecov test coverage](https://codecov.io/gh/junjunlab/omicScope/branch/main/graph/badge.svg)](https://codecov.io/gh/junjunlab/omicScope?branch=main) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

## Overview

**omicScope** is a comprehensive R package for one-stop RNA-seq data analysis and visualization.
Built with the modern S4 object-oriented framework, omicScope streamlines the entire RNA-seq analysis
workflow from mapped bam files or raw counts to biological insights.

## Key Features

-   📦 **Data Normalization**: Multiple normalization methods (CPM, TPM, RPKM/FPKM)
-   📦 **Dimensionality Reduction**: PCA, t-SNE, UMAP analysis with interactive visualizations
-   📦 **Differential Expression**: Integration with DESeq2, edgeR, and limma
-   📦 **Functional Enrichment**: GO, KEGG, and custom gene set enrichment analysis
-   📦 **Activity Inference**: Pathway and transcription factor activity estimation using decoupleR
-   📦 **Rich Visualizations**: Publication-ready plots including heatmaps, volcano plots, and correlation matrices
-   📦 **S4 Object System**: Organized data management and reproducible analysis workflows

## Installation

You can install the development version of omicScope from GitHub:

``` r
# Install from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("junjunlab/omicScope")

# Or using pak (recommended)
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}
pak::pak("junjunlab/omicScope")
```

## Quick Start

Here's a basic workflow demonstrating omicScope's capabilities:

``` r
library(omicScope)

# Load example data or create omicscope object from count matrix
data("example_counts")  # Example count matrix
data("example_metadata") # Sample metadata

# Create omicscope object
os <- create_omicscope(counts = example_counts, 
                       colData = example_metadata,
                       organism = "mouse")

# 1. Data normalization
os <- get_normalized_data(os, method = "cpm")

# 2. Dimensionality reduction
os <- run_pca(os)
os <- run_umap(os)

# Visualize PCA
pca_plot(os, color_by = "treatment")

# 3. Differential expression analysis
os <- run_differential_expression(os, 
                                  method = "deseq2",
                                  design = "~ treatment")

# Volcano plot
volcano_plot(os, comparison = "treat_vs_control")

# 4. Functional enrichment analysis
os <- run_enrichment(os, 
                     method = "gsea",
                     gene_sets = "GO_BP")

# 5. Activity inference
# Pathway activity
os <- infer_activity(os, 
                     input_type = "diff_data",
                     infer_type = "pathway",
                     organism = "mouse")

activity_plot(os)

# Transcription factor activity  
os <- infer_activity(os, 
                     input_type = "diff_data",
                     infer_type = "tf",
                     statistics = "ulm",
                     organism = "mouse")

activity_plot(os, target_tf = c("Pou5f1", "Sox2"))
```
