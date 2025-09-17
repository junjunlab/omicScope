
#' Count RNA-seq reads using featureCounts
#'
#' This function performs read counting for RNA-seq data using the Rsubread::featureCounts
#' function. It takes BAM files specified in the omicscope object and counts reads
#' overlapping with genomic features defined in a GTF annotation file.
#'
#' @param object An omicscope object containing sample information and GTF annotation path
#' @param rna_pairedEnd Logical value indicating whether the RNA-seq data is paired-end.
#'   Default is TRUE
#' @param nThreads Integer specifying the number of threads to use for counting.
#'   Default is 1
#' @param featureCounts_params List of additional parameters to pass to
#'   \code{\link[Rsubread]{featureCounts}}. These will override default parameters
#' @param ... Additional arguments (currently not used)
#'
#' @return An omicscope object with updated assays containing the count matrix and
#'   updated rowData containing only genes present in the count matrix
#'
#' @details
#' The function uses the following default parameters for featureCounts:
#' \itemize{
#'   \item \code{isGTFAnnotationFile = TRUE}: Use GTF format annotation
#'   \item \code{GTF.featureType = "exon"}: Count reads overlapping exons
#'   \item \code{GTF.attrType = "gene_id"}: Use gene_id as primary identifier
#'   \item \code{GTF.attrType.extra = c("gene_name","gene_biotype")}: Extract additional attributes
#' }
#'
#' The resulting count matrix is stored as a sparse matrix (dgCMatrix) to save memory.
#' Column names of the count matrix correspond to sample names from the colData.
#' The rowData is filtered to include only genes present in the final count matrix.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' os <- count_data(os)
#'
#' # For single-end data with multiple threads
#' os <- count_data(os, rna_pairedEnd = FALSE, nThreads = 4)
#'
#' # With custom featureCounts parameters
#' os <- count_data(os,
#'                  nThreads = 8,
#'                  featureCounts_params = list(
#'                    minMQS = 10,
#'                    strandSpecific = 2
#'                  ))
#' }
#'
#' @seealso
#' \code{\link[Rsubread]{featureCounts}} for detailed information about the underlying
#' counting function and its parameters
#'
#' @importFrom Rsubread featureCounts
#' @importFrom SummarizedExperiment colData rowData assays rowData<- assays<-
#' @importFrom methods as
#'
#' @export
setGeneric("count_data",function(object,...){
    standardGeneric("count_data")
})







#' @rdname count_data
#' @export
setMethod("count_data",
          signature(object = "omicscope"),
          function(object,
                   rna_pairedEnd = TRUE,
                   nThreads = 1,
                   featureCounts_params = list()){
              mt <- SummarizedExperiment::colData(object)

              # ==================================================================
              # get counts
              ct <- do.call(Rsubread::featureCounts,
                            modifyList(list(
                                files = mt$sample,
                                isGTFAnnotationFile = TRUE,
                                isPairedEnd = rna_pairedEnd,
                                GTF.featureType = "exon",
                                annot.ext = object@gtfPath,
                                GTF.attrType = "gene_id",
                                GTF.attrType.extra = c("gene_name","gene_biotype"),
                                nthreads = nThreads),
                                featureCounts_params))

              cts <- ct$counts
              colnames(cts) <- rownames(mt)
              cts <- methods::as(as.matrix(cts), "dgCMatrix")

              # reassign rowdata
              rowdata <- SummarizedExperiment::rowData(object)
              rowdata <- subset(rowdata, gene_id %in% rownames(cts))
              rowdata <- rowdata[rownames(cts),]

              # identical(rownames(rowdata),rownames(cts))

              SummarizedExperiment::rowData(object) <- rowdata

              # return
              SummarizedExperiment::assays(object) <- list(counts = cts)


              return(object)
          }
)
