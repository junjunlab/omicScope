
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
#' @param counts Character string specifying the path to a count matrix file,
#'   or NULL. The file should contain gene expression counts with genes as rows
#'   and samples as columns
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
        gtf <- rtracklayer::import.gff(gtfAnno, format = "gtf")

        # feature
        rowData <- data.frame(gtf)[,c("gene_id","gene_name","gene_biotype")] |>
            dplyr::distinct()

        rownames(rowData) <- rowData$gene_id
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
        ct <- read.table(file = counts,skip = "#",
                         header = TRUE,
                         row.names = 1,check.names = FALSE)

        # check featurecounts output?
        if(featureCountOutput == TRUE){
            ct <- ct[,c(8:ncol(ct))]
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
            if(all(c("sample","group") %in% colnames(metadata))){
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
                stop("sample and group must included in metadata!")
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
            if(all(c("sample","group") %in% colnames(metadata))){
                mt <- mt |>
                    dplyr::inner_join(y = metadata, by = "sample")

                rownames(mt) <- mt$sample

                # check sample name
                if(!("sample_name" %in% colnames(mt))){
                    mt$sample_name <- mt$sample
                }
            }else{
                stop("sample and group must included in metadata!")
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
