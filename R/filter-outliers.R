
#' @title Identify and Remove Outlier Samples Using Multi-dimensional PCA
#' @description This function identifies outlier samples based on their positions
#' in a multi-dimensional PCA space. It uses the Mahalanobis distance to find samples
#' that are statistically distant from the center of the data distribution across
#' the specified principal components.
#'
#' @details
#' The function operates as follows:
#' \enumerate{
#'   \item It requires that PCA has already been run on the `omicscope` object via `run_reduction(..., reduction = "pca")`.
#'   \item It calculates the Mahalanobis distance for each sample in the space defined by the selected principal components. This distance metric accounts for the covariance between the PCs.
#'   \item It determines a distance threshold based on a p-value from the Chi-squared distribution. The degrees of freedom for the test are equal to the number of PCs used (`length(dims_to_use)`).
#'   \item Samples with a distance exceeding this threshold are considered outliers.
#'   \item A message is printed identifying the outliers, and a new `omicscope` object is returned with these samples removed.
#' }
#'
#' This function is designed to be used as a data cleaning step. A typical workflow is:
#' run PCA → use `filter_outliers` → run PCA again on the cleaned data → visualize the improved sample clustering.
#'
#' **Key advantages of multi-dimensional filtering:**
#' - Can detect outliers hidden in higher-dimensional space that appear normal in 2D
#' - More statistically rigorous than 2D-only approaches
#' - Automatically adjusts statistical test for the number of dimensions used
#'
#' @param object An `omicscope` object on which `run_reduction` with `reduction = "pca"` has been run.
#' @param dims_to_use A numeric vector specifying which principal components to
#' use for calculating the distance. Defaults to `1:3` (i.e., PC1, PC2, and PC3).
#' Using more dimensions can capture outliers that are hidden in lower-dimensional projections.
#' @param p_value_threshold A numeric value (between 0 and 1) for the significance
#' level. Samples with a Mahalanobis distance corresponding to a p-value less
#' than this threshold will be flagged as outliers. A lower value makes the
#' filtering more stringent. Defaults to `0.05`.
#' @param ... Additional arguments (currently not used).
#'
#'
#' @return A new, filtered `omicscope` object with outlier samples removed.
#' If no outliers are found, the original object is returned unmodified.
#' **Important:** The `@reduction` slot is cleared in the returned object, requiring
#' re-running of `run_reduction()` to generate valid dimensionality reduction results.
#'
#'
#'
#' @examples
#' \dontrun{
#' # === Complete real-world workflow example ===
#' # This example shows the typical analysis pipeline with visualization
#'
#' # Assume 'obj' is your preprocessed omicscope object with normalized counts
#'
#' # Step 1: Run initial PCA for outlier detection
#' obj <- run_reduction(object = obj)
#'
#' # Step 2: Create before-filtering visualization
#' p1 <- dim_plot(obj)
#'
#' # Step 3: Apply multi-dimensional outlier filtering
#' # Using 3D space (PC1, PC2, PC3) with 5% significance threshold
#' obj <- filter_outliers(obj,
#'                        dims_to_use = 1:3,
#'                        p_value_threshold = 0.05)
#'
#' # Step 4: Re-run PCA on the cleaned dataset
#' obj <- run_reduction(object = obj)
#'
#' # Step 5: Create after-filtering visualization
#' p2 <- dim_plot(obj)
#'
#' # Step 6: Compare before and after plots side-by-side
#' library(patchwork)
#' comparison_plot <- p1 + p2
#' print(comparison_plot)
#'
#' # The left plot shows the original data with potential outliers
#' # The right plot shows the cleaned data with improved sample clustering
#' }
#'
#'
#' @importFrom stats cov mahalanobis qchisq
#'
#' @export
setGeneric("filter_outliers",function(object,...){
    standardGeneric("filter_outliers")
})






#' @rdname filter_outliers
#' @export
setMethod("filter_outliers",
          signature(object = "omicscope"),
          function(object,
                   dims_to_use = 1:3,
                   p_value_threshold = 0.05) {

              # --- 1. Sanity Checks ---
              if (!"pca" %in% names(object@reduction)) {
                  stop("PCA has not been run. Please run run_reduction(object, reduction = 'pca') first.")
              }

              # --- 2. Extract PCA Scores ---
              pca_data <- object@reduction[["pca"]]
              pc_scores <- pca_data$x

              # Validate PC indices
              max_available_dims <- ncol(pc_scores)
              if (max(dims_to_use) > max_available_dims) {
                  warning(sprintf("Requested PCs exceed available range, adjusting to first %d PCs", max_available_dims))
                  dims_to_use <- 1:max_available_dims
              }

              # Ensure at least 2 PCs are used
              if (length(dims_to_use) < 2) {
                  dims_to_use <- 1:min(2, max_available_dims)
              }

              # Get scores for the selected dimensions
              scores_subset <- pc_scores[, dims_to_use, drop = FALSE]

              # --- 3. Calculate Mahalanobis Distance ---
              # The mahalanobis function calculates the squared Mahalanobis distance.
              # It requires the data, the center, and the covariance matrix.
              data_center <- colMeans(scores_subset)
              data_cov <- stats::cov(scores_subset)

              distances <- stats::mahalanobis(x = scores_subset,
                                              center = data_center,
                                              cov = data_cov)

              # --- 4. Determine the Cutoff Threshold ---
              # The squared Mahalanobis distance for a p-dimensional normal data follows
              # a Chi-squared distribution with p degrees of freedom (here p=2).
              # We find the value of the Chi-squared distribution for a given p-value.
              # `df = 2` because we are using two dimensions (PC1 and PC2).
              cutoff <- stats::qchisq(p = 1 - p_value_threshold, df = length(dims_to_use))

              # --- 5. Identify and Report Outliers ---
              is_outlier <- distances > cutoff
              outlier_samples <- rownames(scores_subset)[is_outlier]
              samples_to_keep <- rownames(scores_subset)[!is_outlier]

              if (length(outlier_samples) > 0) {
                  message(paste("Identified", length(outlier_samples), "outlier(s) based on PCA (p <", p_value_threshold, "):"))
                  message(paste(outlier_samples, collapse = ", "))

                  message(sprintf("\nRemoved %d outlier samples, retained %d normal samples",
                                  length(outlier_samples), length(samples_to_keep)))
              } else {
                  message("No outliers identified with the current threshold.")
                  return(object) # Return the original object if no outliers
              }

              # --- 6. Filter the original object ---

              # Subset the SummarizedExperiment object
              object_filtered <- object[, samples_to_keep]

              # Clear invalidated reduction results
              object_filtered@reduction <- list()

              message("Reduction results cleared, please re-run run_reduction()")
              message("\nReturning a new object with outliers removed.")

              return(object_filtered)
          })
