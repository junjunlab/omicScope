
#' Create an omicscope object for RNA-seq analysis
#'
#' This function creates an omicscope object from various input types including
#' count matrices, BAM files, and GTF annotations. It serves as the main constructor
#' for initializing omicscope objects for comprehensive RNA-seq analysis workflows.
#'
#' @param gtfAnno Character string specifying the path to GTF/GFF annotation file,
#'   or NULL. If provided, gene annotations will be parsed and stored in the object
#' @param bamFile Character vector of BAM file paths for read counting.
#'   Used when starting analysis from aligned reads rather than count matrices
#' @param counts A count matrix, data.frame, or character string specifying
#'   the path to a count file. Accepts various formats:
#'   \itemize{
#'     \item Matrix or data.frame with genes as rows and samples as columns
#'     \item Path to a tab-delimited count file with row names as gene IDs
#'     \item featureCounts output file (specify \code{featureCountOutput = TRUE})
#'   }
#'   Default: \code{NULL}.
#' @param featureCountOutput Logical value indicating whether the count file is
#'   output from featureCounts. If TRUE, the first 7 columns (annotation columns)
#'   will be skipped. Default is FALSE
#' @param metadata Data.frame containing sample metadata. Must include columns
#'   'sample' and 'group'. Optional 'sample_name' column for custom sample names.
#'   If NULL, basic metadata will be inferred from file names
#'
#' @return An omicscope object containing:
#'   \itemize{
#'     \item Expression data in assays slot (if counts provided)
#'     \item Gene annotations in rowData and gtfAnno slots
#'     \item Sample metadata in colData slot
#'     \item Empty slots for downstream analysis results
#'   }
#'
#' @details
#' The function supports three main usage scenarios:
#'
#' \strong{1. From count matrix:}
#' \itemize{
#'   \item Provide \code{counts} file path and optional \code{metadata}
#'   \item Count matrix is loaded and converted to sparse matrix format
#'   \item Gene annotations are matched to count matrix if GTF is provided
#' }
#'
#' \strong{2. From BAM files:}
#' \itemize{
#'   \item Provide \code{bamFile} paths and \code{metadata}
#'   \item Creates object structure for subsequent read counting with \code{count_data()}
#'   \item Requires GTF annotation for gene definitions
#' }
#'
#' \strong{3. Empty object:}
#' \itemize{
#'   \item Call with no arguments to create minimal object
#'   \item Useful for step-by-step object construction
#' }
#'
#' \strong{File Format Requirements:}
#' \itemize{
#'   \item \strong{Count files}: Tab-separated, genes as rows, samples as columns,
#'     first column as gene IDs
#'   \item \strong{featureCounts output}: Standard format with annotation columns 1-7
#'   \item \strong{GTF files}: Standard GTF/GFF format with gene_id, gene_name,
#'     and gene_biotype attributes
#'   \item \strong{BAM files}: Aligned reads in standard BAM format
#' }
#'
#' \strong{Metadata Requirements:}
#' \itemize{
#'   \item \code{sample}: Must match count matrix column names or BAM file paths
#'   \item \code{group}: Experimental groups for differential analysis
#'   \item \code{sample_name}: Optional custom names for samples
#' }
#'
#' The function automatically:
#' \itemize{
#'   \item Converts count matrices to sparse format for memory efficiency
#'   \item Matches gene annotations to expression data
#'   \item Validates sample-metadata correspondence
#'   \item Sets up proper row and column names
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Create from featureCounts output
#' metadata <- data.frame(
#'   sample = c("sample1.bam", "sample2.bam", "sample3.bam", "sample4.bam"),
#'   group = rep(c("control", "treatment"), each = 2),
#'   sample_name = c("ctrl1", "ctrl2", "treat1", "treat2")
#' )
#'
#' os <- omicscope(
#'   gtfAnno = "path/to/annotation.gtf",
#'   counts = "path/to/featureCounts_output.txt",
#'   featureCountOutput = TRUE,
#'   metadata = metadata
#' )
#'
#' # Example 2: Create from custom count matrix
#' os <- omicscope(
#'   gtfAnno = "annotation.gtf",
#'   counts = "expression_matrix.txt",
#'   featureCountOutput = FALSE,
#'   metadata = metadata
#' )
#'
#' # Example 3: Create from BAM files for counting
#' bam_files <- c("sample1.sorted.bam", "sample2.sorted.bam",
#'                "sample3.sorted.bam", "sample4.sorted.bam")
#'
#' metadata <- data.frame(
#'   sample = bam_files,
#'   group = rep(c("day0", "day10"), each = 2),
#'   sample_name = c("day0-rep1", "day0-rep2", "day10-rep1", "day10-rep2")
#' )
#'
#' os <- omicscope(
#'   gtfAnno = "annotation.gtf",
#'   bamFile = bam_files,
#'   metadata = metadata
#' )
#'
#' # Then perform read counting
#' os <- count_data(os)
#'
#' # Example 4: Empty object for step-by-step construction
#' os <- omicscope()
#'
#' # Basic inspection
#' print(os)
#' dim(os)
#' colData(os)
#' rowData(os)
#' }
#'
#'
#'
#' @importFrom rtracklayer import.gff
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr inner_join distinct
#' @importFrom methods as
#' @importFrom utils read.table
#'
#' @export
omicscope <- function(gtfAnno = NULL,
                      bamFile = NULL,
                      counts = NULL,
                      featureCountOutput = FALSE,
                      metadata = NULL){
    # ==========================================================================
    # load gtf
    if(!is.null(gtfAnno)){
        if(!exists("rowData.rda")){
            gtf <- rtracklayer::import.gff(gtfAnno, format = "gtf")

            gtf.tmp <- gtf |> data.frame(check.names = FALSE)

            # feature
            if("gene_id" %in% colnames(gtf.tmp)){
                if(!("gene_name" %in% colnames(gtf.tmp))){
                    gtf.tmp$gene_name <- gtf.tmp$gene_id
                }

                if(!("gene_biotype" %in% colnames(gtf.tmp))){
                    if("gene_type" %in% colnames(gtf.tmp)){
                        gtf.tmp$gene_biotype <- gtf.tmp$gene_type
                    }else{
                        gtf.tmp$gene_biotype <- "gene"
                    }
                }

                rowData <- gtf.tmp[,c("gene_id","gene_name","gene_biotype")] |>
                    dplyr::distinct()
            }else{
                stop("gene_id column not in your gtf file!")
            }

            rownames(rowData) <- rowData$gene_id

            save(rowData, file = "rowData.rda")
        }else{
            load("rowData.rda")
        }

    }else{
        gtf <- GenomicRanges::GRanges(
            seqnames = character(0),
            ranges = IRanges(),
            strand = character(0),
            gene_id = character(0),
            score = numeric(0)
        )

        rowData <- NULL
    }

    # ==========================================================================
    # load counts
    if(!is.null(counts) && is.null(bamFile)){
        # check counts
        if(inherits(counts,c("data.frame","matrix"))){
            ct <- counts
        }else if(is.character(counts)){
            ct <- read.table(file = counts,skip = "#",
                             header = TRUE,
                             row.names = 1,check.names = FALSE)
        }

        # check featurecounts output?
        if(featureCountOutput == TRUE){
            ct <- ct[,c(8:ncol(ct))]
        }

        if(is.null(gtfAnno)){
            rowData <- data.frame(gene_id = rownames(ct),
                                  gene_name = rownames(ct),
                                  gene_biotype = "gene")

            rownames(rowData) <- rowData$gene_id
        }

        rowData <- subset(rowData, gene_id %in% rownames(ct))

        # keep same order
        rowData <- rowData[match(rownames(ct),rowData$gene_id),]

        # identical(rownames(ct),rowData$gene_id)

        # ======================================================================
        # sample info
        cd <- colnames(ct)
        mt <- data.frame(sample = cd,row.names = cd)

        # check metadata
        if(!is.null(metadata)){
            if(all(c("sample") %in% colnames(metadata))){
                mt <- mt |>
                    dplyr::inner_join(y = metadata, by = "sample")

                rownames(mt) <- mt$sample

                # check sample name
                if(!("sample_name" %in% colnames(mt))){
                    mt$sample_name <- mt$sample
                }

                # order
                mt <- mt[match(colnames(ct),rownames(mt)),]

                rownames(mt) <- mt$sample_name
                colnames(ct) <- mt$sample_name
            }else{
                stop("sample column must be included in metadata!")
            }
        }

        # sce
        sce <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = methods::as(as.matrix(ct), "dgCMatrix")),
            rowData = S4Vectors::DataFrame(rowData,
                                           row.names = rownames(rowData)),
            colData = S4Vectors::DataFrame(mt,
                                           row.names = rownames(mt))
        )

    }else if(is.null(counts) && !is.null(bamFile)){
        # ======================================================================
        # sample info
        mt <- data.frame(sample = bamFile,row.names = bamFile)

        # check metadata
        if(!is.null(metadata)){
            if(all(c("sample") %in% colnames(metadata))){
                mt <- mt |>
                    dplyr::inner_join(y = metadata, by = "sample")

                rownames(mt) <- mt$sample

                # check sample name
                if(!("sample_name" %in% colnames(mt))){
                    mt$sample_name <- mt$sample
                }
            }else{
                stop("sample column must be included in metadata!")
            }
        }

        # sce
        sce <- SummarizedExperiment::SummarizedExperiment(
            rowData = S4Vectors::DataFrame(rowData,
                                           row.names = rownames(rowData)),
            colData = S4Vectors::DataFrame(mt,
                                           row.names = rownames(mt))
        )
    }else{
        sce <- SummarizedExperiment::SummarizedExperiment()
    }

    # return
    return(.omicscope(sce,
                      gtfAnno = gtf,
                      gtfPath = gtfAnno))
}






#' Create omicscope Object from UCSC Xena Data
#'
#' @description
#' This function constructs an \code{omicscope} object by integrating multiple types of
#' genomic data downloaded from UCSC Xena, including gene expression counts, phenotype
#' data, survival data, and optionally GTEx normal tissue data. The function performs
#' data preprocessing, sample filtering, and creates a comprehensive SummarizedExperiment-based
#' object suitable for downstream omics analysis.
#'
#' @param gtf_anno Character string specifying the path to GTF annotation file
#'   (e.g., "gencode.v36.annotation.gtf.gz"). If NULL, an empty GRanges object
#'   will be used. The function will cache processed rowData as "rowData.rda"
#'   to speed up subsequent runs.
#' @param counts_data Character string specifying the path to gene expression counts
#'   file from TCGA (e.g., "TCGA-LUAD.star_counts.tsv.gz"). This parameter is
#'   required and contains log2-transformed expression data that will be converted
#'   back to raw counts using the formula: counts = 2^log2_counts - 1.
#' @param gtex_counts_data Character string specifying the path to GTEx expression
#'   data file (e.g., "gene_reads_v10_lung.gct.gz"). If provided, GTEx normal
#'   samples will be integrated with TCGA data. Default is NULL.
#' @param pheno_data Character string specifying the path to clinical/phenotype
#'   data file (e.g., "TCGA-LUAD.clinical.tsv.gz"). This parameter is required
#'   and contains sample metadata including clinical information.
#' @param survival_data Character string specifying the path to survival data file
#'   (e.g., "TCGA-LUAD.survival.tsv.gz"). If NULL, only phenotype data will be used
#'   for sample metadata. Default is NULL.
#'
#' @return An \code{omicscope} object containing:
#' \itemize{
#'   \item \strong{assays}: Raw gene expression counts as a sparse matrix (dgCMatrix)
#'   \item \strong{rowData}: Gene annotation including gene_id, gene_name, and gene_type
#'   \item \strong{colData}: Sample metadata with group assignment, batch information,
#'     and survival data (if provided)
#'   \item \strong{gtfAnno}: GenomicRanges object with GTF annotations
#'   \item \strong{gtfPath}: Path to the original GTF file
#' }
#'
#' @details
#' The function performs comprehensive data processing in several stages:
#'
#' \subsection{GTF Annotation Processing}{
#'   \itemize{
#'     \item Imports GTF file and extracts gene_id, gene_name, and gene_type
#'     \item Caches processed annotation as "rowData.rda" for efficiency
#'     \item Filters annotation to match genes present in expression data
#'   }
#' }
#'
#' \subsection{Expression Data Processing}{
#'   \itemize{
#'     \item Reads TCGA expression data (assumes log2-transformed input)
#'     \item Converts back to raw counts: \code{counts = 2^log2_counts - 1}
#'     \item Ensures integer count values for downstream analysis
#'     \item Filters genes to match available annotation
#'   }
#' }
#'
#' \subsection{Sample Classification}{
#'   \itemize{
#'     \item Extracts sample type from TCGA barcode (4th component, positions 1-2)
#'     \item Classifies samples: id2 >= 10 → "normal", id2 < 10 → "tumor"
#'     \item Assigns batch label "TCGA" to all TCGA samples
#'   }
#' }
#'
#' \subsection{GTEx Integration (Optional)}{
#'   \itemize{
#'     \item Reads GTEx expression data in GCT format
#'     \item Identifies common genes between TCGA and GTEx datasets
#'     \item Combines expression matrices and metadata
#'     \item Assigns batch label "GTEX" and group label "normal" to GTEx samples
#'   }
#' }
#'
#' @section File Format Requirements:
#' \describe{
#'   \item{GTF file}{Standard GTF format with gene_id, gene_name, and gene_type attributes}
#'   \item{Counts data}{Tab-separated file with genes as rows, samples as columns,
#'     first column containing gene IDs}
#'   \item{Phenotype data}{Tab-separated file with 'sample' column matching count
#'     data column names}
#'   \item{Survival data}{Tab-separated file with 'sample' column for joining
#'     with phenotype data}
#'   \item{GTEx data}{GCT format file (skip first 2 lines) with 'Name' column
#'     for gene IDs}
#' }
#'
#' @section TCGA Sample Type Codes:
#' TCGA sample barcodes follow the pattern: TCGA-XX-XXXX-XXX-XXX-XXXX-XX
#' The 4th component determines sample type:
#' \itemize{
#'   \item 01-09: Primary solid tumors
#'   \item 10-19: Normal tissues
#'   \item 20-29: Control analytes
#' }
#'
#' @note
#' \itemize{
#'   \item The function assumes TCGA expression data is log2-transformed
#'   \item GTEx data file path is currently hardcoded as "gene_reads_v10_lung.gct.gz"
#'   \item Row data caching improves performance for repeated runs with same GTF
#'   \item All data alignment and filtering is performed automatically
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with TCGA data only
#' obj <- ucscZenaToObj(
#'   gtf_anno = "gencode.v36.annotation.gtf.gz",
#'   counts_data = "TCGA-LUAD.star_counts.tsv.gz",
#'   pheno_data = "TCGA-LUAD.clinical.tsv.gz"
#' )
#'
#' # Include survival data
#' obj <- ucscZenaToObj(
#'   gtf_anno = "gencode.v36.annotation.gtf.gz",
#'   counts_data = "TCGA-LUAD.star_counts.tsv.gz",
#'   pheno_data = "TCGA-LUAD.clinical.tsv.gz",
#'   survival_data = "TCGA-LUAD.survival.tsv.gz"
#' )
#'
#' # Full integration with GTEx data
#' obj <- ucscZenaToObj(
#'   gtf_anno = "gencode.v36.annotation.gtf.gz",
#'   counts_data = "TCGA-LUAD.star_counts.tsv.gz",
#'   pheno_data = "TCGA-LUAD.clinical.tsv.gz",
#'   survival_data = "TCGA-LUAD.survival.tsv.gz",
#'   gtex_counts_data = "gene_reads_v10_lung.gct.gz"
#' )
#'
#' # Check the created object
#' dim(obj)
#' table(obj$group)
#' table(obj$batch)
#' }
#'
#'
#'
#'
#' @importFrom rtracklayer import.gff
#' @importFrom data.table fread
#' @importFrom dplyr select mutate across everything pull bind_rows left_join filter
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom methods as
#' @importFrom GenomicRanges GRanges
#' @importFrom utils read.delim
#'
#' @export
ucscZenaToObj <- function(gtf_anno = NULL,
                          counts_data = NULL,
                          gtex_counts_data = NULL,
                          pheno_data = NULL,
                          survival_data = NULL){
    # ==========================================================================
    # load gtf
    if(!is.null(gtf_anno)){
        if(!exists("rowData.rda")){
            gtf <- rtracklayer::import.gff(gtf_anno, format = "gtf")

            gtf.tmp <- gtf |> data.frame(check.names = FALSE)

            # feature
            if("gene_id" %in% colnames(gtf.tmp)){
                if(!("gene_name" %in% colnames(gtf.tmp))){
                    gtf.tmp$gene_name <- gtf.tmp$gene_id
                }

                if(!("gene_biotype" %in% colnames(gtf.tmp))){
                    if("gene_type" %in% colnames(gtf.tmp)){
                        gtf.tmp$gene_biotype <- gtf.tmp$gene_type
                    }else{
                        gtf.tmp$gene_biotype <- "gene"
                    }
                }

                rowData <- gtf.tmp[,c("gene_id","gene_name","gene_biotype")] |>
                    dplyr::distinct()
            }else{
                stop("gene_id column not in your gtf file!")
            }

            rownames(rowData) <- rowData$gene_id

            save(rowData, file = "rowData.rda")
        }else{
            load("rowData.rda")
        }

    }else{
        gtf <- GenomicRanges::GRanges()

        rowData <- NULL
    }
    # ==========================================================================
    # load counts
    if(is.null(counts_data)){
        stop("Please supply counts data for samples!")
    }

    ct.nm <- data.table::fread(file = counts_data,
                               header = TRUE,check.names = FALSE)

    ct <- dplyr::select(ct.nm, -1)
    ct <- 2^ct - 1
    ct <- ct |>
        dplyr::mutate(dplyr::across(dplyr::everything(), as.integer))
    rownames(ct) <- dplyr::pull(ct.nm[,1])

    rowData <- subset(rowData, gene_id %in% rownames(ct))

    # keep same order
    rowData <- rowData[match(rownames(ct),rowData$gene_id),]

    # identical(rownames(ct),rowData$gene_id)

    # ==========================================================================
    # load pheno data
    if(is.null(pheno_data)){
        stop("Please supply pheno data for samples!")
    }

    ph <- read.delim(file = pheno_data,header = TRUE,check.names = FALSE) |>
        dplyr::mutate(sample_name = sample)

    # filter
    ph <- subset(ph, sample %in% colnames(ct))
    rownames(ph) <- ph$sample

    ph <- ph[colnames(ct), ]

    # identical(colnames(ct), rownames(ph))

    nm <- sapply(strsplit(ph$sample,split = "-"),"[",4)
    nm <- substr(nm, start = 1, stop = 2)

    ph$id2 <- nm
    ph <- ph |>
        dplyr::mutate(group = ifelse(id2 >= 10,"Normal","Tumor"),
                      batch = "TCGA")

    # ==========================================================================
    # load survival data
    if(!is.null(survival_data)){
        sur <- read.delim(file = survival_data,header = TRUE,check.names = FALSE) |>
            dplyr::filter(sample %in% rownames(ph))

        mt <- ph |>
            dplyr::left_join(y = sur, by = "sample")

        rownames(mt) <- mt$sample

        mt <- mt[colnames(ct), ]

        # identical(colnames(ct), rownames(mt))
    }else{
        mt <- ph
    }

    # ==========================================================================
    # whether supply gtex data
    if(!is.null(gtex_counts_data)){
        # load gtex data
        gtex <- data.table::fread(gtex_counts_data,skip = 2,header = T)

        gtex.ct <- gtex[,c(-1,-2)]
        rownames(gtex.ct) <- gtex$Name

        # intersect gene name
        gov <- intersect(rownames(gtex.ct), rownames(ct))

        # filter genes
        ct <- data.frame(ct, check.names = FALSE,row.names = rownames(ct))
        ct <- ct[gov,]
        gtex.ct <- data.frame(gtex.ct, check.names = FALSE,row.names = rownames(gtex.ct))
        gtex.ct <- gtex.ct[gov,]

        # identical(rownames(ct), rownames(gtex.ct))
        ct2 <- cbind(ct,gtex.ct)

        rowData <- rowData[gov,]

        gtex.ph <- data.frame(sample = colnames(gtex.ct)) |>
            dplyr::mutate(group = "Normal",batch = "GTEX")
        rownames(gtex.ph) <- gtex.ph$sample

        mt2 <- dplyr::bind_rows(mt, gtex.ph)

        mt2$group2 <- paste(mt2$batch, mt2$group,sep = "-")

        # identical(colnames(ct2), rownames(mt2))
    }else{
        ct2 <- ct
        mt2 <- mt
    }

    # ==========================================================================
    # sce obj
    sce <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = methods::as(as.matrix(ct2), "dgCMatrix")),
        rowData = S4Vectors::DataFrame(rowData,
                                       row.names = rownames(rowData)),
        colData = S4Vectors::DataFrame(mt2,
                                       row.names = rownames(mt2))
    )

    # return
    return(.omicscope(sce,
                      gtfAnno = gtf,
                      gtfPath = gtf_anno))
}







#' Convert TCGAbiolinks Object to OmicScope Object
#'
#' @description
#' This function converts a TCGAbiolinks SummarizedExperiment object into an
#' OmicScope object format. It processes TCGA RNA-seq data and optionally
#' integrates GTEx normal tissue data for comparative analysis. The function
#' handles metadata formatting, counts matrix extraction, survival data
#' processing, and gene annotation.
#'
#' @param gtf_anno Character string. Path to a GTF annotation file. Required
#'   when \code{gtex_counts_data} is provided. The GTF file should contain
#'   gene_id, gene_name, and gene_biotype (or gene_type) columns. Default is NULL.
#' @param tcgabiolinks_obj A SummarizedExperiment object obtained from
#'   TCGAbiolinks (e.g., using \code{GDCprepare()}). Must contain:
#'   \itemize{
#'     \item colData: sample metadata including tissue_type, sample, barcode,
#'       vital_status, days_to_death, days_to_last_follow_up
#'     \item assay "unstranded": raw counts matrix
#'     \item assay "tpm_unstrand": normalized TPM values (optional)
#'     \item rowData: gene annotation with gene_id, gene_name, gene_type
#'   }
#' @param gtex_counts_data Character string or NULL. Path to a GTEx counts
#'   data file (tab-separated format with gene names in columns 1-2 and sample
#'   counts in subsequent columns). When provided, GTEx normal tissue data will
#'   be integrated with TCGA data. Default is NULL.
#'
#' @return An OmicScope object (SummarizedExperiment) containing:
#'   \item{assays}{
#'     \itemize{
#'       \item counts: Combined sparse matrix (dgCMatrix) of TCGA and GTEx
#'         counts (if GTEx provided)
#'       \item normed_counts: TPM normalized counts (only when GTEx not provided)
#'     }
#'   }
#'   \item{colData}{
#'     Sample metadata including:
#'     \itemize{
#'       \item group: Tissue type (Tumor/Normal)
#'       \item batch: Data source (TCGA/GTEX)
#'       \item group2: Combined batch-group identifier
#'       \item OS.time: Overall survival time in days
#'       \item OS: Overall survival status (1=dead, 0=alive)
#'       \item All original TCGA metadata columns
#'     }
#'   }
#'   \item{rowData}{Gene annotation with gene_id, gene_name, gene_biotype}
#'   \item{gtfAnno}{GenomicRanges object from GTF file (empty if GTEx not used)}
#'
#'
#'
#'
#' @examples
#' \dontrun{
#' # Example 1: Convert TCGA data only
#' library(TCGAbiolinks)
#'
#' # Download and prepare TCGA data
#' query <- GDCquery(
#'   project = "TCGA-LUAD",
#'   data.category = "Transcriptome Profiling",
#'   data.type = "Gene Expression Quantification",
#'   workflow.type = "STAR - Counts"
#' )
#' GDCdownload(query)
#' tcga_data <- GDCprepare(query)
#'
#' # Convert to OmicScope object
#' omics_obj <- TCGAbiolinksToObj(tcgabiolinks_obj = tcga_data)
#'
#' # Example 2: Integrate with GTEx data
#' omics_obj_with_gtex <- TCGAbiolinksToObj(
#'   gtf_anno = "gencode.v22.annotation.gtf",
#'   tcgabiolinks_obj = tcga_data,
#'   gtex_counts_data = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
#' )
#'
#' # Access different components
#' counts <- assay(omics_obj_with_gtex, "counts")
#' metadata <- colData(omics_obj_with_gtex)
#' gene_info <- rowData(omics_obj_with_gtex)
#'
#' # Check survival data
#' table(metadata$OS)
#' summary(metadata$OS.time)
#' }
#'
#'
#'
#' @importFrom SummarizedExperiment colData assay rowData SummarizedExperiment
#' @importFrom dplyr rename mutate select distinct bind_rows
#' @importFrom data.table fread
#' @importFrom rtracklayer import.gff
#' @importFrom GenomicRanges GRanges
#' @importFrom methods as
#' @importFrom S4Vectors DataFrame
#'
#' @export
TCGAbiolinksToObj <- function(gtf_anno = NULL,
                              tcgabiolinks_obj = NULL,
                              gtex_counts_data = NULL){
    # ============================================================================
    # format metadata
    mt <- SummarizedExperiment::colData(tcgabiolinks_obj) |>
        data.frame() |>
        dplyr::rename(group = tissue_type,
                      sample2 = sample,
                      sample = barcode) |>
        dplyr::mutate(batch = "TCGA")

    # get counts assay
    ct <- SummarizedExperiment::assay(tcgabiolinks_obj, "unstranded")

    # normalized data
    norm <- SummarizedExperiment::assay(tcgabiolinks_obj, "tpm_unstrand")

    # ===============
    if(is.null(gtex_counts_data)){
        # gene annotation
        rowData <- SummarizedExperiment::rowData(tcgabiolinks_obj) |>
            data.frame(check.names = FALSE) |>
            dplyr::select(gene_id,gene_name,gene_type) |>
            dplyr::rename(gene_biotype = gene_type)

        gtf <- GenomicRanges::GRanges()
    }else{
        if(!exists("rowData.rda")){
            gtf <- rtracklayer::import.gff(gtf_anno, format = "gtf")

            gtf.tmp <- gtf |> data.frame(check.names = FALSE)

            # feature
            if("gene_id" %in% colnames(gtf.tmp)){
                if(!("gene_name" %in% colnames(gtf.tmp))){
                    gtf.tmp$gene_name <- gtf.tmp$gene_id
                }

                if(!("gene_biotype" %in% colnames(gtf.tmp))){
                    if("gene_type" %in% colnames(gtf.tmp)){
                        gtf.tmp$gene_biotype <- gtf.tmp$gene_type
                    }else{
                        gtf.tmp$gene_biotype <- "gene"
                    }
                }

                rowData <- gtf.tmp[,c("gene_id","gene_name","gene_biotype")] |>
                    dplyr::distinct()
            }else{
                stop("gene_id column not in your gtf file!")
            }

            rownames(rowData) <- rowData$gene_id

            save(rowData, file = "rowData.rda")
        }else{
            load("rowData.rda")
        }
    }


    # add survival time
    mt <- mt |>
        dplyr::mutate(OS.time = ifelse(vital_status == "Dead",
                                       days_to_death,days_to_last_follow_up)) |>
        dplyr::mutate(OS = ifelse(vital_status == "Dead", 1, 0))

    # ============================================================================
    # load getx data
    # whether supply gtex data
    if(!is.null(gtex_counts_data)){
        # load gtex data
        gtex <- data.table::fread(gtex_counts_data,skip = 2,header = T)

        gtex.ct <- gtex[,c(-1,-2)]
        rownames(gtex.ct) <- gtex$Name

        # intersect gene name
        gov <- intersect(rownames(gtex.ct), rownames(ct))

        # filter genes
        ct <- data.frame(ct, check.names = FALSE,row.names = rownames(ct))
        ct <- ct[gov,]
        gtex.ct <- data.frame(gtex.ct, check.names = FALSE,row.names = rownames(gtex.ct))
        gtex.ct <- gtex.ct[gov,]

        # identical(rownames(ct), rownames(gtex.ct))
        ct2 <- cbind(ct,gtex.ct)

        rowData <- rowData[gov,]

        gtex.ph <- data.frame(sample = colnames(gtex.ct)) |>
            dplyr::mutate(group = "Normal",batch = "GTEX")
        rownames(gtex.ph) <- gtex.ph$sample

        mt2 <- dplyr::bind_rows(mt, gtex.ph)

        mt2$group2 <- paste(mt2$batch, mt2$group,sep = "-")

        # identical(colnames(ct2), rownames(mt2))

        asy <- list(counts = methods::as(as.matrix(ct2), "dgCMatrix"))
    }else{
        ct2 <- ct
        mt2 <- mt
        gtf_anno <- ""

        asy <- list(counts = methods::as(as.matrix(ct2), "dgCMatrix"),
                    normed_counts = methods::as(norm, "dgCMatrix"))
    }

    # identical(colnames(ct2), rownames(mt2))
    # ==========================================================================
    # sce obj
    sce <- SummarizedExperiment::SummarizedExperiment(
        assays = asy,
        rowData = S4Vectors::DataFrame(rowData,
                                       row.names = rownames(rowData)),
        colData = S4Vectors::DataFrame(mt2,
                                       row.names = rownames(mt2))
    )

    # return
    return(.omicscope(sce,
                      gtfAnno = gtf,
                      gtfPath = gtf_anno))
}
