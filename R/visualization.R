
#' Create Dimensionality Reduction Scatter Plot
#'
#' This function creates a scatter plot visualization of dimensionality reduction
#' results stored in an omicscope object. Currently only supports PCA plotting.
#'
#' @param object An \code{omicscope} S4 object containing dimensionality reduction results.
#'   The object must have reduction data in the \code{reduction} slot
#'   (run \code{run_reduction} first).
#' @param reduction Character string specifying which dimensionality reduction to plot.
#'   Options are:
#'   \itemize{
#'     \item \code{"pca"} - Principal Component Analysis scatter plot (PC1 vs PC2)
#'     \item \code{"umap"} - UMAP scatter plot (currently not implemented)
#'     \item \code{"tsne"} - t-SNE scatter plot (currently not implemented)
#'   }
#'   Default is \code{c("pca","umap","tsne")}, which will use "pca" as the first choice.
#' @param color_by A character string specifying the column name in `colData(object)`
#' to be used for coloring the points on the plot. Defaults to `"group"`.
#' @param ... Additional arguments (currently not used).
#'
#' @return A \code{ggplot2} object showing the PCA scatter plot with samples colored
#'   by sample names.
#'
#'
#'
#' @examples
#' \dontrun{
#' # Run dimensionality reduction first
#' obj <- normalize_data(obj)
#' obj <- run_reduction(obj, reduction = "pca")
#'
#' # Create PCA plot
#' p <- dim_plot(obj, reduction = "pca")
#' print(p)
#'
#' # The plot can be further customized with ggplot2
#' library(ggplot2)
#' p + labs(title = "My PCA Analysis") +
#'   theme(legend.position = "bottom")
#' }
#'
#' @seealso
#' \code{\link{run_reduction}}, \code{\link{normalize_data}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline theme_bw theme
#' @importFrom ggplot2 element_blank element_text scale_color_brewer coord_equal xlab ylab
#' @export
setGeneric("dim_plot",function(object,...){
    standardGeneric("dim_plot")
})






#' @rdname dim_plot
#' @export
setMethod("dim_plot",
          signature(object = "omicscope"),
          function(object,
                   reduction = c("pca","umap","tsne"),
                   color_by = "group"){
              reduction <- match.arg(reduction, choices = reduction)
              # ==================================================================
              # check data
              ck <- length(object@reduction) == 0

              if(ck){
                  stop("Please run run_reduction function first!")
              }

              # Check if specific reduction exists
              if (!reduction %in% names(object@reduction)) {
                  stop("Reduction '", reduction, "' not found. Available reductions: ",
                       paste(names(object@reduction), collapse = ", "))
              }

              # do reduction
              if(reduction == "pca"){
                  data <- object@reduction[["pca"]]

                  # The PC scores are stored in the "x" value of the prcomp object
                  pc_scores <- data.frame(data$x, check.names = FALSE)
                  pc_scores$sample <- rownames(pc_scores)

                  # metadata
                  coldata <- data.frame(SummarizedExperiment::colData(object),
                                        check.names = FALSE,
                                        stringsAsFactors = T)

                  pc_scores <- coldata |>
                      dplyr::inner_join(y = pc_scores, by = "sample")

              }else if(reduction == "umap"){

              }else if(reduction == "tsne"){

              }

              # plot
              ggplot(pc_scores,aes(x = PC1, y = PC2,color = .data[[color_by]])) +
                  geom_point(size = 1) +
                  geom_hline(yintercept = 0,lty = "dashed") +
                  geom_vline(xintercept = 0,lty = "dashed") +
                  theme_bw() +
                  theme(panel.grid = element_blank(),
                        aspect.ratio = 1,
                        axis.text = element_text(colour = "black")) +
                  scale_color_brewer(palette = "Set2") +
                  xlab("Dim 1") + ylab("Dim 2") +
                  guides(color = guide_legend(override.aes = list(size = 4)))
          }
)





#' Create Volcano Plot for Differential Expression Results
#'
#' This function creates volcano plots to visualize differential expression results
#' stored in an omicscope object. It displays log2 fold change vs -log10(p-value)
#' with color-coded significance levels and optional gene labeling.
#'
#' @param object An \code{omicscope} S4 object containing differential expression results.
#'   The object must have differential expression data in the \code{diffExpData} slot
#'   (run \code{differential_expression} first).
#' @param method Character string specifying which differential expression method
#'   results to plot. Options are:
#'   \itemize{
#'     \item \code{"deseq2"} - Plot DESeq2 results
#'     \item \code{"edger"} - Plot edgeR results
#'     \item \code{"limma"} - Plot limma results
#'   }
#'   Default is \code{c("deseq2","edger","limma")}, which will use "deseq2" as first choice.
#' @param marker_gene A character vector of gene names to highlight on the plot.
#'   If \code{NULL} (default), the top significant genes will be labeled instead.
#' @param top_gene An integer specifying how many top significant genes to label
#'   when \code{marker_gene} is NULL. Default is 10. Set to 0 to disable gene labeling.
#'   Top genes are selected based on lowest p-values within each significance category.
#' @param gene_label_size Numeric value specifying the size of gene name labels.
#'   Default is 3.
#' @param gene_number_size Numeric value specifying the size of the gene count
#'   labels (e.g., "sigUp: 150"). Default is 4.
#' @param gene_number_label_pos A numeric vector of length 2 specifying the position
#'   of gene count labels as normalized coordinates (0-1). Default is \code{c(0.95, 0.95)}
#'   for top-right corner positioning.
#' @param point_size Numeric value specifying the size of points in the scatter plot.
#'   Default is 1.
#' @param color A character vector of length 3 specifying colors for significantly
#'   upregulated, non-significant, and significantly downregulated genes respectively.
#'   Default is \code{c("#AF1740", "grey", "#074799")} (red, grey, blue).
#' @param nrow An integer specifying the number of rows when multiple comparisons
#'   are plotted. If \code{NULL} (default), \code{cowplot} will determine the layout automatically.
#' @param ... Additional arguments (currently not used).
#'
#' @return A \code{ggplot2} object (for single comparison) or a \code{cowplot} grid
#'   object (for multiple comparisons) showing volcano plots with the following features:
#'   \itemize{
#'     \item Points colored by significance status
#'     \item Horizontal dashed line at p-value threshold
#'     \item Vertical dashed lines at log2FC thresholds
#'     \item Optional gene labels for top significant or specified marker genes
#'     \item Gene count annotations in corners
#'   }
#'
#'
#' @examples
#' \dontrun{
#' # Basic volcano plot with default settings
#' p1 <- volcano_plot(obj, method = "deseq2")
#' print(p1)
#'
#' # Plot without gene labels
#' p2 <- volcano_plot(obj, method = "edger", top_gene = 0)
#' print(p2)
#'
#' # Highlight specific marker genes
#' markers <- c("GAPDH", "ACTB", "TP53", "MYC")
#' p3 <- volcano_plot(obj,
#'                    method = "limma",
#'                    marker_gene = markers)
#' print(p3)
#'
#' # Customize appearance
#' p4 <- volcano_plot(obj,
#'                    method = "deseq2",
#'                    top_gene = 5,
#'                    point_size = 1.5,
#'                    gene_label_size = 4,
#'                    color = c("red", "darkgrey", "blue"))
#' print(p4)
#'
#' # Multiple comparisons in custom layout
#' p5 <- volcano_plot(obj, nrow = 2)
#' print(p5)
#' }
#'
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline theme_bw theme
#' @importFrom ggplot2 element_blank element_text guides guide_legend scale_color_manual labs
#' @importFrom dplyr slice_min
#' @importFrom cowplot plot_grid
#' @export
setGeneric("volcano_plot",function(object,...){
    standardGeneric("volcano_plot")
})






#' @rdname volcano_plot
#' @export
setMethod("volcano_plot",
          signature(object = "omicscope"),
          function(object,
                   method = c("deseq2","edger","limma"),
                   marker_gene = NULL,
                   top_gene = 10,
                   gene_label_size = 3,
                   gene_number_size = 4,
                   gene_number_label_pos = c(0.95,0.95),
                   point_size = 1,
                   color = c("#AF1740", "grey", "#074799"),
                   nrow = NULL){
              method <- match.arg(method,choices = c("deseq2","edger","limma"))
              # ==================================================================
              # check data
              ck <- length(object@diffExpData) == 0

              if(ck){
                  stop("Please run differential_expression function first!")
              }

              # plot

              # get data
              diff.data <- object@diffExpData[[method]]

              exp.names <- names(diff.data)

              # loop plot
              # x = 1
              lapply(seq_along(diff.data),function(x){
                  tmp.data <- diff.data[[x]]

                  # ============================================================
                  # mark gene
                  if(is.null(marker_gene)){
                      # filter top genes
                      up <- subset(tmp.data@data,type == "sigUp") |>
                          dplyr::slice_min(order_by = pvalue,n = top_gene)

                      dn <- subset(tmp.data@data,type == "sigDown") |>
                          dplyr::slice_min(order_by = pvalue,n = top_gene)
                  }else{
                      mk <- subset(tmp.data@data, gene_name %in% marker_gene)

                      up <- subset(mk,type == "sigUp")
                      dn <- subset(mk,type == "sigDown")
                  }

                  # label pos
                  posx <- round(range(tmp.data@data[,"log2FoldChange"]),digits = 0)

                  # count sig genes
                  upn <- paste("sigUp:",tmp.data@sigUp)
                  dnn <- paste("sigDown:",tmp.data@sigDown)


                  # ============================================================
                  p <-
                      ggplot(tmp.data@data) +
                      geom_point(aes(x = log2FoldChange,y = -log10(pvalue), color = type), size = point_size) +
                      geom_hline(yintercept = -log10(tmp.data@pvalueThreshold),lty = "dashed") +
                      geom_vline(xintercept = c(-tmp.data@log2FCthreshold,tmp.data@log2FCthreshold),lty = "dashed") +
                      theme_bw() +
                      theme(panel.grid = element_blank(),
                            plot.title = element_text(face = "bold"),
                            axis.text = element_text(colour = "black")) +
                      guides(color = guide_legend(override.aes = list(size = 4))) +
                      scale_color_manual(values = c(sigUp = color[1],nonSig = color[2],sigDown = color[3])) +
                      labs(title = exp.names[x])


                  if(requireNamespace(c("ggrepel","ggpp"), quietly = TRUE)){
                      p2 <- p +
                          # ================ sig gene numbers
                          ggpp::geom_text_npc(data = data.frame(),
                                              aes(npcx = gene_number_label_pos[1],npcy = gene_number_label_pos[2],label = upn),
                                              color = color[1],size = gene_number_size) +
                          ggpp::geom_text_npc(data = data.frame(),
                                              aes(npcx = 1-gene_number_label_pos[1],npcy = gene_number_label_pos[2],label = dnn),
                                              color = color[3],size = gene_number_size) +
                          # ================ sig gene labels
                          ggrepel::geom_text_repel(data = up,
                                                   aes(x = log2FoldChange,y = -log10(pvalue),label = gene_name),
                                                   direction = "y",
                                                   nudge_x = posx[2] - up[,"log2FoldChange"],
                                                   hjust = 1,
                                                   size = gene_label_size,
                                                   fontface = "italic",
                                                   segment.linetype = "dashed",
                                                   segment.colour = "grey",
                                                   max.overlaps = Inf) +
                          ggrepel::geom_text_repel(data = dn,
                                                   aes(x = log2FoldChange,y = -log10(pvalue),label = gene_name),
                                                   direction = "y",
                                                   nudge_x = -posx[2] - dn[,"log2FoldChange"],
                                                   hjust = 0,
                                                   size = gene_label_size,
                                                   fontface = "italic",
                                                   segment.linetype = "dashed",
                                                   segment.colour = "grey",
                                                   max.overlaps = Inf)
                  }else{
                      warning("Package 'ggrepel', 'ggpp' is needed for this function to work.")
                  }

                  return(p2)
              }) -> plist

              cowplot::plot_grid(plotlist = plist,align = "hv",nrow = nrow)
          }
)






#' Visualize Biological Activity Results
#'
#' @description
#' Create visualizations for biological activity inference results, including pathway
#' activities and transcription factor activities. The function automatically adapts
#' the plot type based on the number of samples and inference type.
#'
#' @param object An \code{omicscope} object containing activity data in the
#'   \code{@activityData} slot. Activity data must be generated first using
#'   \code{infer_activity()}.
#' @param top_tf Integer. Number of top variable transcription factors to display
#'   when \code{infer_type = "tf"}. For multi-sample data, selects TFs with highest
#'   standard deviation across samples. For single-sample data, selects top and
#'   bottom scoring TFs. Default is 20.
#' @param target_pathway Character vector. Specific pathway names to visualize.
#'   When provided for single-sample data, creates scatter plots showing pathway
#'   weights vs. gene statistics. Only used when \code{infer_type = "pathway"}.
#' @param target_tf Character vector. Specific transcription factor names to visualize.
#'   When provided for single-sample data, creates scatter plots showing TF regulation
#'   modes vs. differential expression. Only used when \code{infer_type = "tf"}.
#' @param complexHeatmap_params List. Additional parameters to pass to
#'   \code{ComplexHeatmap::Heatmap()} function when creating heatmaps for
#'   multi-sample data. Default is an empty list.
#' @param ... Additional arguments (currently not used)
#'
#' @return
#' Returns different plot objects depending on the data type and parameters:
#' \itemize{
#'   \item \strong{Multi-sample data}: A \code{ComplexHeatmap} object showing activity
#'     scores across samples (z-score normalized)
#'   \item \strong{Single-sample data (overview)}: A \code{ggplot} object showing
#'     activity scores as horizontal bar chart
#'   \item \strong{Single-sample data (specific targets)}: A \code{cowplot} grid
#'     combining multiple scatter plots for specified pathways or TFs
#' }
#'
#'
#'
#' @examples
#' \dontrun{
#' # First perform activity inference
#' os <- infer_activity(os, input_type = "counts", infer_type = "pathway")
#'
#' # Basic pathway activity visualization
#' activity_plot(os)
#'
#' # TF activity with more top factors
#' os_tf <- infer_activity(os, infer_type = "tf", statistics = "ulm")
#' activity_plot(os_tf, top_tf = 40)
#'
#' # Specific pathway detailed view (single-sample data)
#' activity_plot(os, target_pathway = c("p53", "MAPK"))
#'
#' # Specific TF detailed view (single-sample data)
#' activity_plot(os_tf, target_tf = c("Pou5f1", "Mef2a"))
#'
#' # Custom heatmap parameters for multi-sample data
#' activity_plot(os, complexHeatmap_params = list(
#'   cluster_rows = FALSE,
#'   show_row_names = TRUE,
#'   row_names_gp = gpar(fontsize = 10)
#' ))
#' }
#'
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_hline geom_vline theme_bw theme element_blank element_text scale_fill_viridis_b scale_color_manual xlab ylab ggtitle guides guide_legend
#' @importFrom dplyr group_by summarise arrange pull slice filter mutate case_when
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom stats sd
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom utils head
#' @importFrom cowplot plot_grid
#'
#'
#'
#' @export
setGeneric("activity_plot",function(object,...){
    standardGeneric("activity_plot")
})






#' @rdname activity_plot
#' @export
setMethod("activity_plot",
          signature(object = "omicscope"),
          function(object,
                   top_tf = 20,
                   target_pathway = NULL,
                   target_tf = NULL,
                   complexHeatmap_params = list()){

              # ==================================================================
              # get plot data
              if(length(object@activityData) == 0){
                  stop("Please supply activity data by running infer_activity function first!")
              }

              pdt <- object@activityData

              # check input data
              if(ncol(pdt@inputData) > 1){
                  pdf <- pdt@resData

                  # Transform to wide matrix
                  pdf_mat <- pdf |>
                      tidyr::pivot_wider(id_cols = 'condition',
                                         names_from = 'source',
                                         values_from = 'score') |>
                      tibble::column_to_rownames('condition')

                  # check inferType
                  if(pdt@inferType == "tf"){
                      # Get top tfs with more variable means across clusters

                      tfs <- pdf |>
                          dplyr::group_by(source) |>
                          dplyr::summarise(std = stats::sd(score)) |>
                          dplyr::arrange(-abs(std)) |>
                          head(top_tf) |>
                          dplyr::pull(source)

                      pdf_mat <- pdf_mat[,tfs]
                  }

                  # Scale per feature
                  pdf_mat <- as.matrix(scale(pdf_mat,center = TRUE,
                                             scale = TRUE))


                  pres <- do.call(ComplexHeatmap::Heatmap,
                                  modifyList(list(matrix = pdf_mat,
                                                  name = "zscore",
                                                  rect_gp = grid::gpar(col = "white"),
                                                  border = TRUE),
                                             complexHeatmap_params))
              }else if(ncol(pdt@inputData) == 1){
                  if(pdt@inferType == "tf"){
                      if(is.null(target_tf)){
                          pdf <- pdt@resData |>
                              dplyr::arrange(score)

                          tfs <- pdf |>
                              dplyr::slice(c(1:top_tf, (nrow(pdf)-top_tf+1):nrow(pdf))) |>
                              dplyr::pull(source)

                          pdf <- subset(pdf, source %in% tfs)

                          # order
                          pdf$source <- factor(pdf$source, levels = pdf$source)

                          # plot
                          pres <-
                              ggplot(pdf) +
                              geom_col(aes(x = score, y = source, fill = -log10(p_value)),
                                       width = 0.75,color = "black") +
                              geom_vline(xintercept = 0,lty = "dashed", color = "black") +
                              theme_bw() +
                              theme(panel.grid = element_blank(),
                                    axis.text = element_text(colour = "black")) +
                              scale_fill_viridis_b(option = "mako",direction = 1) +
                              ylab("Pathway") + xlab("Activity")
                      }else{
                          # loop for target_tf
                          # x = 1
                          lapply(seq_along(target_tf),function(x){
                              net <-  subset(pdt@netData, source %in% target_tf[x]) |>
                                  tibble::column_to_rownames('target')

                              dif <- pdt@diffData |>
                                  dplyr::filter(gene_name %in% rownames(net)) |>
                                  tibble::column_to_rownames(var = "gene_name")

                              ov <- intersect(rownames(net),rownames(dif))
                              dif <- dif[ov,]
                              net <- net[ov,]

                              plf <- cbind(net,dif) |>
                                  dplyr::mutate(reg = ifelse(mor == 1,
                                                             "Activation",
                                                             "Inhibition"))

                              # plot
                              p <-
                                  ggplot(plf) +
                                  geom_point(aes(x = log2FoldChange, y = -log10(pvalue), color = reg)) +
                                  geom_hline(yintercept = 0,lty = "dashed", color = "black") +
                                  geom_vline(xintercept = 0,lty = "dashed", color = "black") +
                                  theme_bw() +
                                  theme(panel.grid = element_blank(),
                                        aspect.ratio = 1,
                                        axis.text = element_text(colour = "black")) +
                                  scale_color_manual(values = c("Activation" = "#F97A00",
                                                                "Inhibition" = "#239BA7"),
                                                     name = "Regulation") +
                                  ggtitle(label = target_tf[x]) +
                                  guides(color = guide_legend(override.aes = list(size = 4)))

                              return(p)
                          }) -> plist

                          # combine
                          pres <- cowplot::plot_grid(plotlist = plist)
                      }
                  }else{
                      if(is.null(target_pathway)){
                          pdf <- pdt@resData |>
                              dplyr::arrange(score)

                          # order
                          pdf$source <- factor(pdf$source, levels = pdf$source)

                          # plot
                          pres <-
                              ggplot(pdf) +
                              geom_col(aes(x = score, y = source, fill = -log10(p_value)),
                                       width = 0.75,color = "black") +
                              geom_vline(xintercept = 0,lty = "dashed", color = "black") +
                              theme_bw() +
                              theme(panel.grid = element_blank(),
                                    axis.text = element_text(colour = "black")) +
                              scale_fill_viridis_b(option = "mako",direction = 1) +
                              ylab("Pathway") + xlab("Activity")

                      }else{
                          # loop for target_pathway
                          lapply(seq_along(target_pathway),function(x){
                              net <-  subset(pdt@netData, source %in% target_pathway[x]) |>
                                  tibble::column_to_rownames('target')

                              dif <- pdt@inputData

                              ov <- intersect(rownames(dif),rownames(net))

                              net <- net[ov,]
                              dif <- data.frame(dif[ov,], check.names = FALSE)

                              net$y <- dif[,c(1)]
                              net <- net |>
                                  dplyr::mutate(type = dplyr::case_when(weight > 0 & y > 0 ~ "Up",
                                                                        weight < 0 & y < 0 ~ "Down",
                                                                        .default = "non"))

                              # plot
                              p <-
                                  ggplot(net) +
                                  geom_point(aes(x = weight, y = y, color = type)) +
                                  geom_hline(yintercept = 0,lty = "dashed", color = "black") +
                                  geom_vline(xintercept = 0,lty = "dashed", color = "black") +
                                  theme_bw() +
                                  theme(panel.grid = element_blank(),
                                        aspect.ratio = 1,
                                        axis.text = element_text(colour = "black")) +
                                  scale_color_manual(values = c(Up = "#DC143C", Down = "#19183B", non = "grey")) +
                                  ylab("Statistic value") + xlab("Weight") +
                                  ggtitle(label = target_pathway[x]) +
                                  guides(color = guide_legend(override.aes = list(size = 4)))

                              return(p)
                          }) -> plist

                          # combine
                          pres <- cowplot::plot_grid(plotlist = plist)

                      }
                  }

              }


              return(pres)
          }
)





#' Generate Expression Heatmap from Normalized Data
#'
#' @description
#' Create a heatmap visualization of gene expression data using normalized counts.
#' The function displays z-score normalized expression values across samples with
#' optional gene and sample filtering capabilities.
#'
#' @param object An \code{omicscope} object containing normalized expression data.
#'   The normalized data must be available in \code{object@normalizedData$longer}.
#'   Run \code{get_normalized_data()} first if this data is not available.
#' @param selected_gene Character vector. Specific gene names to include in the heatmap.
#'   If \code{NULL} (default), randomly selects 25 genes from the dataset for visualization.
#' @param selected_sample Character vector of sample names to include in the heatmap.
#'   If NULL (default), all samples will be included.
#' @param color_by Character string specifying the column name in colData to use
#'   for sample annotation coloring. Default is "group".
#' @param complexHeatmap_params List. Additional parameters to pass to
#'   \code{ComplexHeatmap::Heatmap()} function for customizing the heatmap appearance.
#'   Default is an empty list.
#' @param ... Additional arguments (currently not used)
#'
#'
#' @return
#' A \code{ComplexHeatmap} object displaying the expression heatmap.
#'
#'
#'
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom grid gpar
#'
#' @examples
#' \dontrun{
#' # Basic heatmap with random 25 genes
#' heatmap1 <- exp_heatmap_plot(omics_obj)
#'
#' # Heatmap with specific genes
#' selected_genes <- c("GENE1", "GENE2", "GENE3")
#' heatmap2 <- exp_heatmap_plot(omics_obj, selected_gene = selected_genes)
#'
#' # Heatmap with specific samples and custom grouping
#' selected_samples <- c("Sample1", "Sample2", "Sample3")
#' heatmap3 <- exp_heatmap_plot(omics_obj,
#'                              selected_gene = selected_genes,
#'                              selected_sample = selected_samples,
#'                              color_by = "treatment")
#'
#' # Custom ComplexHeatmap parameters
#' custom_params <- list(
#'   show_column_dend = FALSE,
#'   column_title = "Gene Expression",
#'   heatmap_legend_param = list(title = "Expression Z-score")
#' )
#' heatmap4 <- exp_heatmap_plot(omics_obj,
#'                              complexHeatmap_params = custom_params)
#'
#' # Display the heatmap
#' ComplexHeatmap::draw(heatmap1)
#' }
#'
#'
#'
#' @export
setGeneric("exp_heatmap_plot",function(object,...){
    standardGeneric("exp_heatmap_plot")
})







#' @rdname exp_heatmap_plot
#' @export
setMethod("exp_heatmap_plot",
          signature(object = "omicscope"),
          function(object,
                   selected_gene = NULL,
                   selected_sample = NULL,
                   color_by = "group",
                   complexHeatmap_params = list()){
              # ==================================================================
              # get normalized counts
              ck <- is.null(object@normalizedData)

              if(ck){
                  stop("Please run get_normalized_data function first!")
              }

              pdf <- object@normalizedData

              # check gene
              if(is.null(selected_gene)){
                  pdf <- subset(pdf, gene_name %in% sample(unique(pdf$gene_name),
                                                           size = 25,
                                                           replace = FALSE))
              }else{
                  pdf <- subset(pdf, gene_name %in% selected_gene)
              }

              # tolong format
              pdf.lg <- pdf |>
                  tidyr::pivot_longer(cols = colnames(pdf)[1:(ncol(pdf) - 3)],
                                      names_to = "sample",
                                      values_to = "value")

              cold <- SummarizedExperiment::colData(object) |>
                  data.frame(check.names = FALSE)

              pdf.anno.lg <- pdf.lg |>
                  dplyr::inner_join(y = cold, by = "sample")

              # check sample
              if(!is.null(selected_sample)){
                  pdf.anno.lg <- subset(pdf.anno.lg, sample %in% selected_sample)
              }

              pdf_mat <- pdf.anno.lg |>
                  dplyr::select(gene_name, sample, value)

              # Transform to wide matrix
              pdf_mat <- pdf_mat |>
                  tidyr::pivot_wider(id_cols = 'gene_name',
                                     names_from = 'sample',
                                     values_from = 'value') |>
                  tibble::column_to_rownames('gene_name') |>
                  t()


              # Scale per feature
              pdf_mat <- as.matrix(scale(pdf_mat,center = TRUE,scale = TRUE))

              gp <- pdf.anno.lg[,c("sample",color_by)] |> unique() |> data.frame(check.names = FALSE)
              rownames(gp) <- gp$sample
              gp <- gp[rownames(pdf_mat),]

              raw.anno <- ComplexHeatmap::rowAnnotation(Sample = gp$group,
                                                        show_annotation_name = FALSE)

              # plot
              pres <- do.call(ComplexHeatmap::Heatmap,
                              modifyList(list(matrix = pdf_mat,
                                              name = "Z-score",
                                              cluster_columns = TRUE,
                                              # rect_gp = grid::gpar(col = "white"),
                                              na_col = "black",
                                              left_annotation = raw.anno,
                                              show_row_names = F,
                                              column_names_gp = grid::gpar(fontface = "italic"),
                                              border = TRUE),
                                         complexHeatmap_params))

              return(pres)
          }
)







#' Generate Correlation Plots for Expression Data
#'
#' @description
#' Create correlation visualizations for normalized expression data. The function can
#' generate either a correlation matrix heatmap for all samples or a scatter plot
#' comparing two specific samples.
#'
#' @param object An \code{omicscope} object containing normalized expression data.
#'   The normalized data must be available in \code{object@normalizedData}.
#'   Run \code{get_normalized_data()} first if this data is not available.
#' @param x Character string. Name of the sample to use for x-axis in scatter plot.
#'   If \code{NULL} (default), generates a correlation matrix instead of scatter plot.
#' @param y Character string. Name of the sample to use for y-axis in scatter plot.
#'   If \code{NULL} (default), generates a correlation matrix instead of scatter plot.
#' @param color Character string. Color for points in scatter plot. Default is "grey20".
#'   Only used when both \code{x} and \code{y} are specified.
#' @param point_size Numeric. Size of points in scatter plot. Default is 0.5.
#'   Only used when both \code{x} and \code{y} are specified.
#' @param corrplot_params List. Additional parameters to pass to
#'   \code{corrplot::corrplot()} function when generating correlation matrix.
#'   Default is an empty list.
#' @param ... Additional arguments (currently not used)
#'
#'
#' @return
#' Returns different plot objects depending on the input parameters:
#' \itemize{
#'   \item \strong{Correlation Matrix} (when \code{x} or \code{y} is \code{NULL}):
#'     A correlation matrix heatmap showing pairwise correlations between all samples
#'   \item \strong{Scatter Plot} (when both \code{x} and \code{y} are specified):
#'     A \code{ggplot} object showing sample-to-sample expression correlation with
#'     optional correlation statistics
#' }
#'
#'
#' @importFrom stats cor
#' @importFrom corrplot corrplot COL1
#' @importFrom ggplot2 ggplot aes geom_point geom_abline theme_bw theme element_blank element_text coord_equal scale_x_continuous scale_y_continuous xlab ylab
#' @importFrom scales label_log
#' @importFrom SummarizedExperiment colData
#'
#' @examples
#' \dontrun{
#' # First ensure normalized data is available
#' os <- get_normalized_data(os)
#'
#' # Generate correlation matrix for all samples
#' correlation_plot(os)
#'
#' # Compare two specific samples
#' correlation_plot(os,
#'                  x = "day0-rep1",
#'                  y = "day0-rep2")
#'
#' # Scatter plot with custom appearance
#' correlation_plot(os,
#'                  x = "day0-rep1",
#'                  y = "day10-rep1",
#'                  color = "red",
#'                  point_size = 1.0)
#'
#' # Custom correlation matrix parameters
#' correlation_plot(os,
#'                  corrplot_params = list(
#'                    method = "circle",
#'                    type = "upper",
#'                    order = "hclust",
#'                    addCoef.col = "white"
#'                  ))
#' }
#'
#'
#'
#' @export
setGeneric("correlation_plot",function(object,...){
    standardGeneric("correlation_plot")
})







#' @rdname correlation_plot
#' @export
setMethod("correlation_plot",
          signature(object = "omicscope"),
          function(object,
                   x = NULL, y = NULL,
                   color = "grey20", point_size = 0.5,
                   corrplot_params = list()){
              # ==================================================================
              # get normalized counts
              ck <- length(object@normalizedData$longer) == 0

              if(ck){
                  stop("Please run get_normalized_data function first!")
              }

              cd <- SummarizedExperiment::colData(object)

              pdf.w <- object@normalizedData$wider[,1:nrow(cd)]

              colnames(pdf.w) <- cd$sample_name

              # check type
              if(is.null(x) | is.null(y)){
                  # calculate cor
                  m <- stats::cor(pdf.w)

                  # plot
                  do.call(corrplot::corrplot,
                          modifyList(
                              list(
                                  corr = m,
                                  method = "square",
                                  type = "lower",
                                  is.corr = FALSE,
                                  col = corrplot::COL1(sequential = "Blues"),
                                  cl.ratio = 0.25,
                                  addCoef.col = "black"
                              ),corrplot_params
                          ))
              }else{
                  p <-
                      ggplot(pdf.w,
                             aes(x = .data[[x]] + 1,y = .data[[y]] + 1)) +
                      geom_point(color = color,size = point_size) +
                      geom_abline(intercept = 0,slope = 1,lty = "dashed") +
                      theme_bw(base_size = 12) +
                      theme(panel.grid = element_blank(),
                            axis.text = element_text(colour = "black")) +
                      coord_equal() +
                      scale_x_continuous(transform = "log2", guide = "axis_logticks",
                                         labels = scales::label_log(base = 2,digits = 1)) +
                      scale_y_continuous(transform = "log2", guide = "axis_logticks",
                                         labels = scales::label_log(base = 2,digits = 1)) +
                      xlab(x) + ylab(y)

                  if(requireNamespace("ggpubr", quietly = TRUE)){
                      p + ggpubr::stat_cor()
                  }else{
                      warning("Package 'ggpubr' is needed for this function to work.")
                  }
              }

          }
)


