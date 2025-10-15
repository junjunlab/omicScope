globalVariables(c("FDR", "P.Value", "PC1", "PC2", "PValue", "adj.P.Val", "gene_name",
                  "log2FoldChange", "logFC", "gene_id","mor", "p_value", "reg", "score",
                  "std", "weight", "y","group", "value",".data", "id2","OS", "OS.time",
                  "modifyList", "na.omit", "pvalue", "sample_name", "stat", "type","alpha",
                  "barcode", "days_to_death", "days_to_last_follow_up", "gene_type",
                  "tissue_type","vital_status", "signature_tme","cell_type", "prop",
                  ".", "end", "exon_length", "label", "pos", "rpm", "strand", "target_genes",
                  "transcript_id", "transcript_length", "transcript_rank", "which_label",
                  "x", "xend", "y_max","gene_span", "gene_start", "track", "trans_end",
                  "trans_span", "trans_start"))



#' S4 Class for Storing Activity Inference Results
#'
#' @description
#' A container class that stores results from biological activity inference
#' (pathway or transcription factor activities) using decoupleR methods.
#'
#' @slot inferType Character. Type of inference performed ("pathway" or "tf").
#' @slot netData Data.frame. The regulatory network used for inference, containing
#'   source, target, and regulatory information.
#' @slot inputData Matrix. Processed input data used for activity calculation
#'   (genes as rows, samples as columns).
#' @slot diffData Data.frame. Original differential expression results when using
#'   diff_data input type. Empty data.frame when using counts input.
#' @slot resData Data.frame. Activity inference results from decoupleR containing
#'   source, condition, score, and statistical method information.
#'
#' @details
#' This class is automatically created by the \code{infer_activity} method and
#' stored in the \code{@activityData} slot of an \code{omicscope} object.
#' The \code{resData} slot contains the main results for downstream analysis.
#'
#'
#'
#' @seealso
#' \code{\link{infer_activity}}, \code{\link[methods]{setClass}},
#' \code{\link[decoupleR]{run_mlm}}
#'
#' @name activitydata-class
#' @rdname activitydata-class
#' @exportClass activitydata
activitydata <- setClass("activitydata",
                         slots = list(inferType = "character",
                                      netData = "data.frame",
                                      inputData = "matrix",
                                      diffData = "data.frame",
                                      resData = "data.frame"))






#' The omicscope class for comprehensive RNA-seq analysis
#'
#' The omicscope class extends SummarizedExperiment to provide a unified framework
#' for RNA-seq data analysis workflows, including read counting, normalization,
#' differential expression, dimensionality reduction, and functional enrichment analysis.
#'
#' @slot gtfAnno A \code{\link[GenomicRanges]{GRanges}} object containing gene
#'   annotations parsed from GTF file, including gene coordinates, biotypes, and metadata
#' @slot gtfPath Character vector specifying the file path to the GTF annotation file
#'   used for gene annotation and read counting
#' @slot normalizedData List containing normalized expression matrices with different
#'   methods. Names typically include "TPM", "FPKM", "CPM", "vst", "rlog", etc.
#' @slot reduction List storing dimensionality reduction results from methods like
#'   PCA, t-SNE, UMAP. Each element contains coordinates and method parameters
#' @slot diffExpData List organized by analysis method (deseq2, edger, limma),
#'   containing differential expression results for different contrasts
#' @slot enrichmentData List storing functional enrichment analysis results from
#'   various databases (GO, KEGG, Reactome, etc.)
#' @slot activityData An \code{activitydata} object containing pathway activity
#'   analysis results and related metadata
#' @slot tmeData List storing tumor microenvironment analysis results from IOBR
#'   for various methods and contrasts.
#'
#'
#' @author Jun Zhang
#' @name omicscope
#' @rdname omicscope
#' @aliases omicscope
#'
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom methods setClass new
#' @importFrom S4Vectors SimpleList
#'
#'
#' @export omicscope
.omicscope <- setClass("omicscope",
                       slots = list(gtfAnno = "GRanges",
                                    gtfPath = "character",
                                    normalizedData = "data.frame",
                                    reduction = "SimpleList",
                                    diffExpData = "SimpleList",
                                    enrichmentData = "SimpleList",
                                    activityData = "activitydata",
                                    tmeData = "SimpleList"),
                       contains = "SummarizedExperiment")








#' An S4 Class to Store Differential Expression Analysis Results
#'
#' The \code{diffdata} class is designed to store and organize results from
#' differential expression analysis, including statistical parameters,
#' significance thresholds, and summary statistics.
#'
#' @slot contrastName A character string specifying the name of the contrast
#'   or comparison being analyzed (e.g., "Treatment_vs_Control", "GroupA_vs_GroupB").
#' @slot method A character string indicating the statistical method used for
#'   differential expression analysis (e.g., "DESeq2", "edgeR", "limma", "t-test").
#' @slot design A list containing the experimental design information, including
#'   factors, covariates, and design matrix specifications used in the analysis.
#' @slot log2FCthreshold A numeric value specifying the log2 fold change threshold
#'   used to determine statistical significance (typically 0, 0.5, or 1).
#' @slot pvalueThreshold A numeric value specifying the p-value or adjusted p-value
#'   threshold used to determine statistical significance (typically 0.01, 0.05).
#' @slot sigUp An integer indicating the number of significantly upregulated
#'   genes/features (positive log2FC above threshold and p-value below threshold).
#' @slot sigDown An integer indicating the number of significantly downregulated
#'   genes/features (negative log2FC below threshold and p-value below threshold).
#' @slot nonSig An integer indicating the number of non-significant genes/features
#'   that do not meet the specified thresholds.
#' @slot data A data.frame containing the complete differential expression results,
#'   typically including columns such as log2FoldChange, pvalue, padj, baseMean, etc.
#'
#' @details
#' The \code{diffdata} class provides a structured way to store differential
#' expression analysis results along with the parameters and summary statistics.
#' This allows for easy tracking of analysis parameters and quick access to
#' summary information without recalculating from the raw results.
#'
#' The \code{data} slot should contain a data.frame with gene/feature identifiers
#' as row names and statistical results as columns. Common columns include:
#' \itemize{
#'   \item \code{log2FoldChange} - Log2 fold change values
#'   \item \code{pvalue} - Raw p-values
#'   \item \code{padj} - Adjusted p-values (e.g., FDR, Bonferroni)
#'   \item \code{baseMean} - Base mean expression levels
#'   \item \code{stat} - Test statistics
#' }
#'
#'
#'
#' @name diffdata-class
#' @rdname diffdata-class
#' @exportClass diffdata
diffdata <- setClass("diffdata",
                     slots = list(contrastName = "character",
                                  method = "character",
                                  design = "list",
                                  log2FCthreshold = "numeric",
                                  pvalueThreshold = "numeric",
                                  sigUp = "integer",
                                  sigDown = "integer",
                                  nonSig = "integer",
                                  data = "data.frame"))




