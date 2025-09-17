

#' Perform differential expression analysis
#'
#' This function performs differential expression analysis using DESeq2, edgeR, or limma-voom
#' methods. It processes RNA-seq count data and identifies significantly differentially
#' expressed genes based on specified thresholds.
#'
#' @param object An omicscope object containing count data and sample information
#' @param selectedSample Character vector of sample names to include in the analysis.
#'   If NULL (default), all samples are used
#' @param contrastName Character string specifying the name for this contrast.
#'   If NULL, defaults to "treat_vs_control"
#' @param deseq2Design Formula specifying the design matrix for DESeq2 analysis.
#'   Default is ~group
#' @param deseq2Contrast Character vector of length 3 specifying the contrast for DESeq2:
#'   c(factor_name, numerator_level, denominator_level).
#'   Default is c('group', 'treat', 'control')
#' @param edgerDesign Design matrix for edgeR analysis. If NULL (default),
#'   uses ~group formula based on colData. Should have rownames matching sample names
#' @param limmaDesign Design matrix for limma analysis. If NULL (default),
#'   uses ~group formula based on colData. Should have rownames matching sample names
#' @param limmaApproach Character specifying the limma approach to use.
#'   One of "voomLmFit" (default), "trend", or "voom"
#' @param log2FCthreshold Numeric value for log2 fold change threshold.
#'   Default is 1 (corresponding to 2-fold change)
#' @param pvalueThreshold Numeric value for p-value significance threshold.
#'   Default is 0.05
#' @param method Character specifying the differential expression method.
#'   One of "deseq2" (default), "edger", or "limma"
#' @param ... Additional arguments (currently not used)
#'
#' @return An omicscope object with differential expression results stored in
#'   the diffExpData slot under the specified method. Results include:
#'   \itemize{
#'     \item Statistical test results (log2FoldChange, pvalue, padj)
#'     \item Gene annotations merged with results
#'     \item Significance classification (sigUp, sigDown, nonSig)
#'     \item Analysis parameters and summary statistics
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts count data and sample metadata from the omicscope object
#'   \item Filters samples if selectedSample is specified
#'   \item Performs differential expression analysis using the selected method:
#'     \itemize{
#'       \item \strong{DESeq2}: DESeqDataSetFromMatrix -> DESeq -> results pipeline
#'       \item \strong{edgeR}: DGEList -> filtering -> normalization -> glmQLFit -> glmQLFTest
#'       \item \strong{limma}: DGEList -> normalization -> voom/trend -> lmFit -> eBayes
#'     }
#'   \item Merges results with gene annotations from rowData
#'   \item Classifies genes based on significance thresholds
#'   \item Stores results in a diffdata object within the omicscope object
#' }
#'
#' For \strong{limma analysis}, three approaches are available:
#' \itemize{
#'   \item \code{voomLmFit}: Uses edgeR::voomLmFit with sample weights (recommended)
#'   \item \code{voom}: Standard limma::voom transformation
#'   \item \code{trend}: Uses log-CPM with eBayes trend=TRUE
#' }
#'
#' For \strong{edgeR and limma}, the function performs:
#' \itemize{
#'   \item Gene filtering using filterByExpr
#'   \item TMM normalization via calcNormFactors
#'   \item Robust dispersion estimation (edgeR only)
#'   \item Quasi-likelihood F-tests (edgeR) or empirical Bayes (limma)
#' }
#'
#' The significance classification is based on:
#' \itemize{
#'   \item \code{sigUp}: log2FoldChange >= log2FCthreshold AND pvalue < pvalueThreshold
#'   \item \code{sigDown}: log2FoldChange <= -log2FCthreshold AND pvalue < pvalueThreshold
#'   \item \code{nonSig}: All other genes
#' }
#'
#' @note
#' \itemize{
#'   \item The function requires 'counts' assay to be present in the omicscope object
#'   \item Sample filtering is based on the 'sample_name' column in colData
#'   \item All character columns in colData are automatically converted to factors
#'   \item For custom design matrices, ensure rownames match sample names after filtering
#' }
#'
#' @examples
#' \dontrun{
#' # Basic DESeq2 analysis with default parameters
#' os <- differential_expression(os)
#'
#' # edgeR analysis with custom contrast name and thresholds
#' os <- differential_expression(os,
#'                              method = "edger",
#'                              contrastName = "treatment_vs_control",
#'                              log2FCthreshold = 1.5,
#'                              pvalueThreshold = 0.01)
#'
#' # DESeq2 analysis with specific samples and contrast
#' os <- differential_expression(os,
#'                              method = "deseq2",
#'                              selectedSample = c("day0-rep1", "day0-rep2",
#'                                               "day10-rep1", "day10-rep2"),
#'                              deseq2Contrast = c('group', 'day10', 'day0'))
#'
#' # edgeR with custom design matrix
#' group <- factor(c("control", "control", "treated", "treated"))
#' custom_design <- model.matrix(~group)
#' os <- differential_expression(os,
#'                              method = "edger",
#'                              edgerDesign = custom_design)
#'
#' # limma-voom analysis
#' os <- differential_expression(os,
#'                              method = "limma",
#'                              limmaApproach = "voom",
#'                              contrastName = "day10_vs_day0")
#'
#' # Access results
#' diff_results <- os@diffExpData$deseq2$treat_vs_control@data
#' }
#'
#'
#'
#' @importFrom SummarizedExperiment assayNames assay rowData colData
#' @importFrom dplyr filter rename mutate mutate_if case_when group_by summarise n
#' @importFrom methods new
#' @importFrom stats model.matrix
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom limma voom lmFit eBayes topTable
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest voomLmFit cpm
#'
#' @export
setGeneric("differential_expression",function(object,...){
    standardGeneric("differential_expression")
})






#' @rdname differential_expression
#' @export
setMethod("differential_expression",
          signature(object = "omicscope"),
          function(object,
                   selectedSample = NULL,
                   contrastName = NULL,
                   deseq2Design = ~group,
                   deseq2Contrast = c('group', 'treat', 'control'),
                   edgerDesign = NULL,
                   limmaDesign = NULL,
                   limmaApproach = c("voomLmFit","trend","voom"),
                   log2FCthreshold = 1,
                   pvalueThreshold = 0.05,
                   method = c("deseq2","edger","limma")){
              limmaApproach <- match.arg(
                  limmaApproach,
                  choices = c("voomLmFit","trend","voom"))
              method <- match.arg(method,choices = c("deseq2","edger","limma"))
              # ================================================================
              # get counts
              ck <- "counts" %in% SummarizedExperiment::assayNames(object)

              if(!ck){
                  stop("Please supply counts data!")
              }

              asy <- SummarizedExperiment::assay(object,"counts")


              # get anno
              ga <- SummarizedExperiment::rowData(object) |>
                  data.frame(check.names = FALSE)

              if(is.null(contrastName)){
                  contrastName <- "treat_vs_control"
              }

              # metadata
              coldata <- data.frame(SummarizedExperiment::colData(object),
                                    check.names = FALSE,
                                    stringsAsFactors = T)

              # filter sample
              if(!is.null(selectedSample)){
                  coldata <- subset(coldata, sample_name %in% selectedSample)
                  asy <- asy[,rownames(coldata)]
              }


              coldata <- coldata |>
                  dplyr::mutate_if(is.character, as.factor)
              # ================================================================
              # diff anslysis
              if(method == "deseq2"){
                  cat("Performing deseq2 differential analysis now...")

                  design <- deseq2Design

                  # metadata
                  col.data <- coldata

                  # DESeqDataSet
                  dds <- DESeq2::DESeqDataSetFromMatrix(
                      countData = asy,
                      colData = col.data,
                      design = deseq2Design)

                  # diff
                  dds <- DESeq2::DESeq(dds)

                  # get results
                  res <- DESeq2::results(dds,
                                         contrast = deseq2Contrast)

                  # filter NA
                  res_all <- data.frame(res, check.names = FALSE) |>
                      dplyr::filter(log2FoldChange != 'NA')

              }else if(method == "edger"){
                  # ============================================================
                  cat("Performing edger differential analysis now...")

                  # edger
                  dge <- edgeR::DGEList(counts = asy, group = coldata$group)

                  gde.fl <- dge[edgeR::filterByExpr(dge), ,
                                keep.lib.sizes = FALSE]

                  gde.fl <- edgeR::calcNormFactors(gde.fl)

                  # design matrix
                  if(is.null(edgerDesign)){
                      group <- coldata$group

                      design <- stats::model.matrix(~group)
                  }else{
                      design <- edgerDesign
                  }

                  rownames(design) <- colnames(gde.fl)

                  gde.fl <- edgeR::estimateDisp(gde.fl, design, robust=TRUE)

                  fit <- edgeR::glmQLFit(y = gde.fl,
                                         design = design, robust = TRUE)

                  qlf <- edgeR::glmQLFTest(glmfit = fit)

                  res_all <- data.frame(edgeR::topTags(object = qlf,n = Inf),
                                        check.names = FALSE) |>
                      dplyr::rename(log2FoldChange = logFC,
                                    pvalue = PValue,
                                    padj = FDR)

              }else if(method == "limma"){
                  # ============================================================
                  cat("Performing limma differential analysis now...")

                  # limma
                  dge <- edgeR::DGEList(counts = asy, group = coldata$group)

                  gde.fl <- dge[edgeR::filterByExpr(dge), ,
                                keep.lib.sizes = FALSE]

                  gde.fl <- edgeR::calcNormFactors(gde.fl)

                  # design matrix
                  if(is.null(limmaDesign)){
                      group <- coldata$group

                      design <- stats::model.matrix(~group)
                  }else{
                      design <- limmaDesign
                  }

                  rownames(design) <- colnames(gde.fl)

                  # check approach
                  if(limmaApproach == "trend"){
                      input <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
                      fit <- limma::lmFit(input, design)
                      fit <- limma::eBayes(fit, trend = TRUE)
                  }else if(limmaApproach == "voom"){
                      input <- limma::voom(dge, design)
                      fit <- limma::lmFit(input, design)
                      fit <- limma::eBayes(fit)
                  }else if(limmaApproach == "voomLmFit"){
                      fit <- edgeR::voomLmFit(dge, design = design,
                                              sample.weights = TRUE)
                      fit <- limma::eBayes(fit)
                  }


                  res_all <- limma::topTable(fit, coef = ncol(design),
                                             number = Inf) |>
                      dplyr::rename(log2FoldChange = logFC,
                                    pvalue = P.Value,
                                    padj = adj.P.Val)

              }

              # add annotation
              res_all$gene_id <- rownames(res_all)
              res_all_anno <- merge(res_all,ga, by = 'gene_id')

              # add sig type
              res_all_anno <- res_all_anno |>
                  dplyr::mutate(
                      type = dplyr::case_when(
                          log2FoldChange >= log2FCthreshold &
                              pvalue < pvalueThreshold ~ "sigUp",
                          log2FoldChange <= -log2FCthreshold &
                              pvalue < pvalueThreshold ~ "sigDown",
                          .default = "nonSig"))

              # ================================================================
              # return
              # summary
              res.sm <- res_all_anno |>
                  dplyr::group_by(type) |>
                  dplyr::summarise(nm = dplyr::n())

              diff.res <- methods::new("diffdata",
                                       contrastName = contrastName,
                                       method = method,
                                       design = list(design),
                                       log2FCthreshold = log2FCthreshold,
                                       pvalueThreshold = pvalueThreshold,
                                       sigUp = subset(res.sm, type == "sigUp")$nm,
                                       sigDown = subset(res.sm, type == "sigDown")$nm,
                                       nonSig = subset(res.sm, type == "nonSig")$nm,
                                       data = res_all_anno)

              sdata <- list(diff.res)
              names(sdata) <- contrastName

              object@diffExpData[[method]] <- sdata

              return(object)
          }
)

