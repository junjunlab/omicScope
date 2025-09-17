

#' Infer Pathway or Transcription Factor Activity
#'
#' This function infers pathway or transcription factor (TF) activity from gene expression data
#' using the decoupleR framework. It supports both count data and differential expression data
#' as input, and can analyze either pathway activities or transcription factor activities.
#'
#' @param object An \code{omicscope} object containing gene expression data and metadata.
#' @param input_type Character string specifying the type of input data. Options are:
#'   \itemize{
#'     \item \code{"counts"}: Use raw count data from the \code{counts} assay
#'     \item \code{"diff_data"}: Use differential expression statistics
#'   }
#' @param diff_data_obj A differential expression data object. Required when
#'   \code{input_type = "diff_data"}. Should be from \code{object@diffExpData}
#'   (e.g., \code{object@diffExpData$deseq2$treat_vs_control}).
#' @param organism Character string specifying the organism. Supported options:
#'   \code{"human"}, \code{"mouse"}, \code{"rat"}. Default is \code{"human"}.
#' @param statistics Character string specifying the statistical method to use.
#'   Default is \code{"mlm"} (multivariate linear model). Other options include
#'   \code{"ulm"} (univariate linear model), \code{"wsum"} (weighted sum), etc.
#' @param infer_type Character string specifying what to infer. Options are:
#'   \itemize{
#'     \item \code{"pathway"}: Infer pathway activities using PROGENy networks
#'     \item \code{"tf"}: Infer transcription factor activities using CollecTRI networks
#'   }
#' @param decouple_params A list of additional parameters to pass to
#'   \code{decoupleR::decouple()}. Default is an empty list.
#' @param use_local_netdata Logical indicating whether to use local network data
#'   bundled with the package (\code{TRUE}) or download from online databases
#'   (\code{FALSE}). Default is \code{FALSE}.
#' @param ... Additional arguments (currently not used)
#'
#'
#' @return An \code{omicscope} object with activity inference results stored in
#'   the \code{@activityData} slot as an \code{activitydata} object.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts gene annotations and validates input data
#'   \item Processes input data based on \code{input_type}:
#'     \itemize{
#'       \item For \code{"counts"}: Converts to CPM (counts per million)
#'       \item For \code{"diff_data"}: Extracts appropriate statistics based on the
#'             differential expression method (DESeq2: stat, edgeR: -log10(p)*log2FC,
#'             limma: t-statistic)
#'     }
#'   \item Loads appropriate network data:
#'     \itemize{
#'       \item For \code{infer_type = "pathway"}: Uses PROGENy pathway-gene networks
#'       \item For \code{infer_type = "tf"}: Uses CollecTRI transcription factor-gene networks
#'     }
#'   \item Runs activity inference using \code{decoupleR::decouple()}
#'   \item Returns the updated \code{omicscope} object with results
#' }
#'
#' @section Network Data Sources:
#' \itemize{
#'   \item \strong{PROGENy}: Pathway RespOnsive GENes for activity inference
#'   \item \strong{CollecTRI}: A comprehensive resource of transcriptional regulation
#' }
#'
#' @section Statistical Methods:
#' \itemize{
#'   \item \strong{mlm}: Multivariate Linear Model - fits all regulators simultaneously
#'   \item \strong{ulm}: Univariate Linear Model - fits each regulator independently
#'   \item \strong{wsum}: Weighted Sum - simple weighted aggregation of target genes
#' }
#'
#' @note
#' \itemize{
#'   \item When using \code{input_type = "counts"}, the function requires a \code{"counts"}
#'     assay in the \code{omicscope} object
#'   \item For \code{input_type = "diff_data"}, the function automatically selects
#'     appropriate statistics based on the differential expression method
#'   \item The \code{mlm} method may fail with highly correlated regulators. Consider
#'     using \code{statistics = "ulm"} or \code{statistics = "wsum"} as alternatives
#'   \item Gene names are made unique using \code{make.unique()} to avoid duplicates
#' }
#'
#' @examples
#' \dontrun{
#' # Infer pathway activities from count data
#' os_pathway <- infer_activity(os,
#'                             input_type = "counts",
#'                             infer_type = "pathway",
#'                             organism = "mouse",
#'                             use_local_netdata = TRUE)
#'
#' # Infer pathway activities from differential expression data
#' os_pathway <- infer_activity(os,
#'                             input_type = "diff_data",
#'                             diff_data_obj = os@diffExpData$deseq2$treat_vs_control,
#'                             infer_type = "pathway",
#'                             organism = "mouse")
#'
#' # Infer transcription factor activities using ULM method
#' os_tf <- infer_activity(os,
#'                        input_type = "counts",
#'                        infer_type = "tf",
#'                        organism = "mouse",
#'                        statistics = "ulm",
#'                        use_local_netdata = TRUE)
#'
#' # Use custom parameters for decoupleR
#' os_custom <- infer_activity(os,
#'                            input_type = "counts",
#'                            infer_type = "pathway",
#'                            organism = "human",
#'                            decouple_params = list(minsize = 10,
#'                                                  consensus_score = TRUE))
#' }
#'
#'
#'
#' @importFrom SummarizedExperiment rowData assayNames assay
#' @importFrom dplyr filter mutate select
#' @importFrom tibble column_to_rownames
#' @importFrom edgeR cpm
#' @importFrom Matrix rowSums
#' @importFrom stats na.omit
#' @importFrom methods new
#' @export
setGeneric("infer_activity",function(object,...){
    standardGeneric("infer_activity")
})






#' @rdname infer_activity
#' @export
setMethod("infer_activity",
          signature(object = "omicscope"),
          function(object,
                   input_type = c("counts","diff_data"),
                   diff_data_obj = NULL,
                   organism = c("human", "mouse","rat"),
                   statistics = "mlm",
                   infer_type = c("pathway","tf"),
                   decouple_params = list(),
                   use_local_netdata = FALSE){
              organism <- match.arg(organism,choices = c("human", "mouse","rat"))
              infer_type <- match.arg(infer_type,choices = c("pathway","tf"))
              input_type <- match.arg(input_type,choices = c("counts","diff_data"))
              # ================================================================
              # get anno
              ga <- SummarizedExperiment::rowData(object) |>
                  data.frame(check.names = FALSE) |>
                  na.omit()

              ga$gene_name <- make.unique(ga$gene_name)

              # check input data
              if(input_type == "counts"){
                  # get counts
                  ck <- "counts" %in% SummarizedExperiment::assayNames(object)

                  if(!ck){
                      stop("Please supply counts data!")
                  }

                  asy <- SummarizedExperiment::assay(object,"counts")

                  asy <- asy[Matrix::rowSums(asy) > 0,]
                  com <- intersect(rownames(asy),rownames(ga))
                  asy <- asy[com,]
                  ga <- ga[com,]

                  rownames(asy) <- ga$gene_name
                  colnames(asy) <- SummarizedExperiment::colData(object)$sample_name

                  input <- edgeR::cpm(y = asy)
                  diffd <- data.frame()
              }else if(input_type == "diff_data"){
                  # check method
                  input <- diff_data_obj@data |>
                      dplyr::filter(!is.na(gene_name)) |>
                      dplyr::mutate(gene_name = make.unique(gene_name))

                  diffd <- input

                  if(diff_data_obj@method == "deseq2"){
                      input <- input |>
                          dplyr::select(gene_name, stat)
                  }else if(diff_data_obj@method == "edger"){
                      input <- input |>
                          dplyr::mutate(stat = -log10(pvalue)*log2FoldChange) |>
                          dplyr::select(gene_name, stat)
                  }else if(diff_data_obj@method == "limma"){
                      input <- input |>
                          dplyr::select(gene_name, t)
                  }

                  input <- input |>
                      tibble::column_to_rownames(var = "gene_name") |>
                      as.matrix()
              }

              # ================================================================
              # PROGENy model

              OmnipathR::omnipath_set_cachedir(tempdir())
              OmnipathR::omnipath_set_console_loglevel('trace')

              # check infer_type
              if(infer_type == "pathway"){
                  if(use_local_netdata == TRUE){
                      if(organism == "human"){
                          net <- readRDS(file = system.file("extdata/human_pathway.rds", package = "omicScope"))
                      }else{
                          net <- readRDS(file = system.file("extdata/mouse_pathway.rds", package = "omicScope"))
                      }

                  }else{
                      net <- decoupleR::get_progeny(organism = organism,
                                                    top = 500)
                  }
              }else{
                  if(use_local_netdata == TRUE){
                      if(organism == "human"){
                          net <- readRDS(file = system.file("extdata/human_tf.rds", package = "omicScope"))
                      }else{
                          net <- readRDS(file = system.file("extdata/mouse_tf.rds", package = "omicScope"))
                      }

                  }else{
                      net <- decoupleR::get_collectri(organism = organism,
                                                      split_complexes = FALSE)
                  }
              }

              # =================================================================

              # Run mlm
              # Run mlm
              if(infer_type == "pathway"){
                  args = list(mlm = list(.mor = "weight"))
              }else{
                  args = list(mlm = list(.mor = "mor"))
              }

              dc.data <- do.call(
                  decoupleR::decouple,modifyList(
                      list(
                          mat = input,
                          network  = net,
                          .source = 'source',
                          .target = 'target',
                          statistics = statistics,
                          consensus_score = FALSE,
                          minsize = 5,
                          args = args
                      ),
                      decouple_params
                  )
              )

              # return
              res <- methods::new(Class = "activitydata",
                                  inferType = infer_type,
                                  netData = net,
                                  inputData = input,
                                  diffData = data.frame(diffd,check.names = FALSE),
                                  resData = dc.data)

              object@activityData <- res

              return(object)
          }
)
