
#' Perform Dimensionality Reduction on omicscope Object
#'
#' This function performs dimensionality reduction analysis (PCA, UMAP, or t-SNE)
#' on normalized expression data stored in an omicscope object. The function applies
#' z-score normalization before reduction and stores the results in the object's
#' reduction slot.
#'
#' @param object An \code{omicscope} S4 object containing normalized expression data.
#'   The object must have "normed_counts" in its assay slots (run \code{normalize_data} first).
#' @param reduction Character string specifying the dimensionality reduction method to use.
#'   Options are:
#'   \itemize{
#'     \item \code{"pca"} - Principal Component Analysis
#'     \item \code{"umap"} - Uniform Manifold Approximation and Projection
#'     \item \code{"tsne"} - t-Distributed Stochastic Neighbor Embedding
#'   }
#'   Default is \code{c("pca","umap","tsne")}, which will use "pca" as the first choice.
#' @param umapParams A list of parameters to pass to \code{\link[umap]{umap}} function.
#'   Common parameters include:
#'   \itemize{
#'     \item \code{n_neighbors} - Number of nearest neighbors (default: 15)
#'     \item \code{min_dist} - Minimum distance between points (default: 0.1)
#'     \item \code{n_epochs} - Number of training epochs (default: 200)
#'     \item \code{metric} - Distance metric (default: "euclidean")
#'   }
#' @param tsneParams A list of parameters to pass to \code{\link[Rtsne]{Rtsne}} function.
#'   Common parameters include:
#'   \itemize{
#'     \item \code{perplexity} - Perplexity parameter (default: 30)
#'     \item \code{max_iter} - Maximum number of iterations (default: 1000)
#'     \item \code{theta} - Speed/accuracy trade-off (default: 0.5)
#'     \item \code{dims} - Output dimensionality (default: 2)
#'   }
#' @param ... Additional arguments (currently not used).
#'
#' @return Returns the input \code{omicscope} object with dimensionality reduction
#'   results stored in the \code{reduction} slot. The reduction results can be
#'   accessed using \code{object@reduction[[reduction_method]]}.
#'
#'
#' For PCA, the function uses \code{\link[stats]{prcomp}}.
#' For UMAP, the function uses \code{\link[umap]{umap}}.
#' For t-SNE, the function uses \code{\link[Rtsne]{Rtsne}}.
#'
#'
#' @examples
#' \dontrun{
#' # Basic PCA
#' obj <- run_reduction(obj, reduction = "pca")
#'
#' # UMAP with custom parameters
#' obj <- run_reduction(obj,
#'                      reduction = "umap",
#'                      umapParams = list(n_neighbors = 15,
#'                                        min_dist = 0.5,
#'                                        n_epochs = 500))
#'
#' # t-SNE with custom parameters
#' obj <- run_reduction(obj,
#'                      reduction = "tsne",
#'                      tsneParams = list(perplexity = 50,
#'                                        max_iter = 1500,
#'                                        theta = 0.3))
#'
#' # Access results
#' pca_result <- obj@reduction[["pca"]]
#' umap_result <- obj@reduction[["umap"]]
#' tsne_result <- obj@reduction[["tsne"]]
#' }
#'
#'
#' @importFrom stats prcomp
#' @importFrom utils modifyList
#'
#' @export
setGeneric("run_reduction",function(object,...){
    standardGeneric("run_reduction")
})






#' @rdname run_reduction
#' @export
setMethod("run_reduction",
          signature(object = "omicscope"),
          function(object,
                   reduction = c("pca","umap","tsne"),
                   umapParams = list(),
                   tsneParams = list()){
              reduction <- match.arg(reduction, choices = reduction)
              # ==================================================================
              # get normalized counts
              ck <- "normed_counts" %in% SummarizedExperiment::assayNames(object)

              if(!ck){
                  stop("Please run normalize_data function first!")
              }

              asy <- data.frame(as.matrix(SummarizedExperiment::assay(object,"normed_counts")),
                                check.names = FALSE)


              asy <- asy[rowSums(asy) > 0, ]

              # zscore
              asy.zs <- t(scale(t(asy), center = TRUE,scale = TRUE))

              # do reduction
              if(reduction == "pca"){
                  reduc <- stats::prcomp(x = t(asy.zs))
              }else if(reduction == "umap"){
                  if (!requireNamespace("umap", quietly = TRUE)) {
                      stop("Package 'umap' is required for UMAP reduction.\n",
                           "Please install it with: install.packages('umap')")
                  }

                  # reduc <- umap::umap(d = t(asy.zs))
                  reduc <- do.call(umap::umap,
                                   modifyList(list(d = t(asy.zs)),
                                              umapParams))
              }else if(reduction == "tsne"){
                  if (!requireNamespace("Rtsne", quietly = TRUE)) {
                      stop("Package 'Rtsne' is required for t-SNE reduction.\n",
                           "Please install it with: install.packages('Rtsne')")
                  }

                  # reduc <- Rtsne::Rtsne(X = t(asy.zs))
                  reduc <- do.call(Rtsne::Rtsne,
                                   modifyList(list(X = t(asy.zs)),
                                              tsneParams))
              }


              # return
              object@reduction[[reduction]] <- reduc

              return(object)
          }
)
