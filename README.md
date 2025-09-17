<!-- badges: start -->

# omicScope <img src="man/figures/logo.png" align="right" width="130"/>

[![Codecov test coverage](https://codecov.io/gh/junjunlab/omicScope/branch/main/graph/badge.svg)](https://codecov.io/gh/junjunlab/omicScope?branch=main) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![GitHub issues](https://img.shields.io/github/issues/junjunlab/omicScope)](https://github.com/junjunlab/omicScope/issues)

<!-- badges: end -->


## Overview

**omicScope** is a comprehensive R package for one-stop RNA-seq data analysis and visualization. Built on the **SummarizedExperiment** data structure and modern S4 object-oriented framework, omicScope provides a unified interface for the entire RNA-seq analysis workflow - from mapped BAM files or raw count matrices to meaningful biological insights. By leveraging the well-established Bioconductor infrastructure, omicScope ensures seamless integration with the broader ecosystem of genomic analysis tools.

## Key Features

-   📦 **Data Normalization**: Multiple normalization methods (CPM, TPM, RPKM/FPKM)
-   📦 **Dimensionality Reduction**: PCA, t-SNE, UMAP analysis with interactive visualizations
-   📦 **Differential Expression**: Integration with DESeq2, edgeR, and limma
-   📦 **Functional Enrichment**: GO, KEGG, and custom gene set enrichment analysis
-   📦 **Activity Inference**: Pathway and transcription factor activity estimation using decoupleR
-   📦 **Rich Visualizations**: Publication-ready plots including heatmaps, volcano plots, and correlation matrices
-   📦 **S4 Object System**: Organized data management and reproducible analysis workflows

> ⚠️ **Development Status**: omicScope is currently under active development. We welcome feedback, suggestions, and contributions from the community. If you have ideas for new features or improvements, please feel free to submit pull requests or open issues!

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

## Contributing

omicScope is an open-source project and we actively encourage community contributions! Whether you're interested in:

-   💡 **Bug fixes** - Help us identify and resolve issues.
-   💡 **New features** - Propose and implement new analytical capabilities.
-   💡 **Documentation** - Improve examples, tutorials, or function documentation.
-   💡 **Testing** - Add unit tests or test with your own datasets.
-   💡 **Ideas** - Share suggestions for improvements or new directions.

## Quick Start

Here's a basic workflow demonstrating omicScope's capabilities:

``` r
library(omicScope)

# Load example data or create omicscope object from count matrix
bams <- c("../test-bam/0a.sorted.bam","../test-bam/0b.sorted.bam",
          "../test-bam/4a.sorted.bam","../test-bam/4b.sorted.bam",
          "../test-bam/10a.sorted.bam","../test-bam/10b.sorted.bam")

mta <- data.frame(sample = bams,
                  sample_name = c("day0-rep1","day0-rep2",
                                  "day4-rep1","day4-rep2",
                                  "day10-rep1","day10-rep2"),
                  group = rep(c("day0","day4","day10"),each = 2))


data("counts")

# Create omicscope object
os <- omicscope(gtfAnno = "Mus_musculus.GRCm38.102.gtf.gz",
                counts = counts,
                metadata = mta)

# 1. Data normalization
os <- get_normalized_data(os, method = "tpm")

# 2. Dimensionality reduction
os <- run_pca(os)


# Visualize PCA
pca_plot(os)

# 3. Differential expression analysis
os <- run_differential_expression(os, 
                                  method = "deseq2",
                                  selectedSample = c("day0-rep1","day0-rep2",
                                                     "day10-rep1","day10-rep2"),
                                  deseq2Contrast = c('group', 'day10', 'day0')
                                  )

# Volcano plot
volcano_plot(os, comparison = "treat_vs_control")

# 4. Functional enrichment analysis
library(org.Mm.eg.db)

os <- run_enrichment(os, 
                     enrich_type = "go",
                     OrgDb = org.Mm.eg.db)

# 5. Activity inference
# Pathway activity
os <- infer_activity(os, 
                     input_type = "counts",
                     infer_type = "pathway",
                     organism = "mouse",
                     use_local_netdata = TRUE)

activity_plot(os)

# Transcription factor activity  
os <- infer_activity(os, 
                     input_type = "counts",
                     infer_type = "tf",
                     statistics = "ulm",
                     organism = "mouse",
                     use_local_netdata = TRUE)

activity_plot(os, target_tf = c("Pou5f1", "Sox2"))
```
