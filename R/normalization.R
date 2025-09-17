
#' Normalize expression data for omicscope objects
#'
#' This function performs normalization on expression count data using different
#' methods including TPM (Transcripts Per Million), RPKM (Reads Per Kilobase Million),
#' and log1p transformation.
#'
#' @param object An \code{omicscope} object containing expression count data
#' @param norm_type Character string specifying the normalization method.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"tpm"}: Transcripts Per Million normalization
#'     \item \code{"rpkm"}: Reads Per Kilobase Million normalization
#'     \item \code{"log1p"}: CPM followed by log1p transformation
#'   }
#' @param ... Additional arguments (currently not used)
#'
#'
#' @return An \code{omicscope} object with updated assays containing both
#'   original counts and normalized counts
#'
#' @details
#' The function performs different normalizations based on the selected method:
#' \itemize{
#'   \item \strong{TPM}: First calculates RPK (reads per kilobase), then scales
#'     by the sum of all RPK values per sample and multiplies by 10^6
#'   \item \strong{RPKM}: First calculates CPM (counts per million), then
#'     divides by gene length in kilobases
#'   \item \strong{log1p}: Calculates CPM values, then applies log1p transformation
#'     commonly used in single-cell analysis
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Normalize using TPM
#' normalized_obj <- normalize_data(omics_obj, norm_type = "tpm")
#'
#' # Normalize using log1p transformation
#' normalized_obj <- normalize_data(omics_obj, norm_type = "log1p")
#'
#' # Access normalized counts
#' norm_counts <- SummarizedExperiment::assay(normalized_obj, "normed_counts")
#' }
#'
#'
#' @export
setGeneric("normalize_data",function(object,...){
    standardGeneric("normalize_data")
})




#' @rdname normalize_data
#' @export
setMethod("normalize_data",
          signature(object = "omicscope"),
          function(object,
                   norm_type = c("tpm","rpkm","log1p")){
              norm_type <- match.arg(norm_type,choices = c("tpm","rpkm","log1p"))
              # ================================================================
              # get counts
              count <- SummarizedExperiment::assay(object)

              # get gene length
              gl <- reduced_gene_length(gtf = object@gtfAnno)
              rownames(gl) <- gl$gene_id

              # keep same order
              gl <- gl[match(rownames(count),gl$gene_id),]

              # identical(rownames(count),gl$gene_id)

              # ================================================================
              if(norm_type == "tpm"){
                  rpk <- as.matrix(count)/(gl$length/1000)
                  norm <- t(t(rpk)/colSums(rpk)*10^6)
              }else if(norm_type == "rpkm"){
                  norm <- t(t(as.matrix(count))/colSums(as.matrix(count))*10^6)
                  norm <- norm/(gl$length/1000)
              }else if(norm_type == "log1p"){
                  cpm <- t(t(as.matrix(count))/colSums(as.matrix(count))*10^4)
                  norm <- log1p(cpm + 1)
              }

              norm <- methods::as(norm, "dgCMatrix")

              # return
              SummarizedExperiment::assays(object) <- S4Vectors::SimpleList(
                  counts = count,
                  normed_counts = norm)

              return(object)
          }
)






#' Extract and format normalized expression data
#'
#' This function extracts normalized expression data from an \code{omicscope}
#' object and formats it into both wide and long formats for downstream analysis
#' and visualization.
#'
#' @param object An \code{omicscope} object containing normalized expression data.
#'   The object must have been processed with \code{\link{normalize_data}} first.
#' @param ... Additional arguments (currently not used)
#'
#'
#' @return An \code{omicscope} object with updated \code{normalizedData} slot
#'   containing:
#'   \itemize{
#'     \item \code{wider}: A data frame in wide format with genes as rows,
#'       samples as columns, plus gene annotation columns
#'     \item \code{longer}: A data frame in long format with columns for
#'       sample, value, and all gene annotation information
#'   }
#'
#'
#'
#' @examples
#' \dontrun{
#' # First normalize the data
#' omics_obj <- normalize_data(omics_obj, norm_type = "tpm")
#'
#' # Then extract formatted normalized data
#' omics_obj <- get_normalized_data(omics_obj)
#'
#' # Access wide format data
#' wide_data <- omics_obj@normalizedData$wider
#'
#' # Access long format data for plotting
#' long_data <- omics_obj@normalizedData$longer
#'
#' }
#'
#' @seealso \code{\link{normalize_data}}
#'
#' @importFrom tidyr pivot_longer
#'
#' @export
setGeneric("get_normalized_data",function(object,...){
    standardGeneric("get_normalized_data")
})




#' @rdname get_normalized_data
#' @export
setMethod("get_normalized_data",
          signature(object = "omicscope"),
          function(object){
              # ================================================================
              # get normalized counts
              ck <- "normed_counts" %in%
                  SummarizedExperiment::assayNames(object)

              if(!ck){
                  stop("Please run get_normalized_data function first!")
              }

              asy <- SummarizedExperiment::assay(object,"normed_counts")


              # get anno
              ga <- SummarizedExperiment::rowData(object) |>
                  data.frame(check.names = FALSE)


              asy.anno <- cbind(as.matrix(asy),ga)

              # tolong format
              asy.anno.lg <- asy.anno |>
                  tidyr::pivot_longer(cols = colnames(asy),
                                      names_to = "sample",
                                      values_to = "value")

              cold <- SummarizedExperiment::colData(object) |>
                  data.frame(check.names = FALSE)

              asy.anno.lg <- asy.anno.lg |>
                  dplyr::inner_join(y = cold, by = "sample")


              object@normalizedData <- list(wider = asy.anno,
                                            longer = asy.anno.lg)

              return(object)
          }
)


