
#' Correct Batch Effect for omicscope Objects
#'
#' @description
#' This method corrects batch effects in gene expression count data using the
#' ComBat-seq algorithm. It operates on omicscope objects containing count data
#' and applies batch effect correction based on batch and group information
#' stored in the object's column metadata.
#'
#' @param object An \code{omicscope} object containing count data in the "counts" assay
#' @param ComBat_seq_params A list of parameters to be passed to the \code{ComBat_seq_devel}
#'   function. Default parameters include counts data, batch information, and group
#'   information from the object's metadata. Additional parameters can be specified
#'   to customize the batch correction process.
#' @param ... Additional arguments (currently not used)
#'
#' @return An \code{omicscope} object with batch-corrected count data. The corrected
#'   counts replace the original counts in the "counts" assay slot.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks that the object contains a "counts" assay
#'   \item Extracts count data and column metadata
#'   \item Applies ComBat-seq batch correction using batch and group information
#'   \item Returns the object with corrected counts
#' }
#'
#' The batch and group information must be present in the object's \code{colData}
#' as columns named "batch" and "group" respectively.
#'
#' @section Required Data:
#' \itemize{
#'   \item The object must contain a "counts" assay
#'   \item Column metadata must include "batch" and "group" variables
#'   \item Count data should be raw count data suitable for ComBat-seq
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' corrected_obj <- correct_batch_effect(my_omicscope_object)
#'
#' # With custom ComBat-seq parameters
#' corrected_obj <- correct_batch_effect(
#'   my_omicscope_object,
#'   ComBat_seq_params = list(shrink = TRUE, shrink.disp = TRUE)
#' )
#' }
#'
#' @seealso
#' \code{\link[sva]{ComBat_seq}} for details on the underlying batch correction algorithm
#'
#'
#' @export
setGeneric("correct_batch_effect",function(object,...){
    standardGeneric("correct_batch_effect")
})






#' @rdname correct_batch_effect
#' @export
setMethod("correct_batch_effect",
          signature(object = "omicscope"),
          function(object,
                   ComBat_seq_params = list()) {
              # get counts
              ck <- "counts" %in% SummarizedExperiment::assayNames(object)

              if(!ck){
                  stop("Please supply counts data!")
              }

              asy <- SummarizedExperiment::assay(object,"counts") |>
                  as.matrix()

              # metadata
              coldata <- data.frame(SummarizedExperiment::colData(object),
                                    check.names = FALSE,
                                    stringsAsFactors = TRUE)

              # identical(colnames(asy),rownames(coldata))

              # batch effect removal
              # btrm <- ComBat_seq_devel(counts = asy,
              #                          batch = coldata$batch,
              #                          group = coldata$group)

              btrm <- do.call(ComBat_seq_devel,
                              modifyList(
                                  list(counts = asy,
                                       batch = coldata$batch,
                                       group = coldata$group),
                                  ComBat_seq_params))

              # return
              SummarizedExperiment::assays(object) <- S4Vectors::SimpleList(
                  counts = methods::as(as.matrix(btrm), "dgCMatrix"))

              return(object)

          })
