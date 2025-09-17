
#' Perform Enrichment Analysis on Differential Expression Results
#'
#' This function performs functional enrichment analysis (Gene Ontology or KEGG pathway)
#' on differential expression data stored in an omicscope object. It supports both
#' over-representation analysis (ORA) and Gene Set Enrichment Analysis (GSEA).
#'
#' @param object An \code{omicscope} S4 object containing differential expression results.
#'   The object must have differential expression data in the \code{diffExpData} slot
#'   (run \code{differential_expression} first).
#' @param diff_data_selected A character vector specifying which differential expression
#'   comparisons to analyze. If \code{NULL} (default), all available comparisons will
#'   be analyzed. Should match names in \code{object@diffExpData}.
#' @param enrich_type Character string specifying the type of enrichment analysis. Options are:
#'   \itemize{
#'     \item \code{"go"} - Gene Ontology over-representation analysis
#'     \item \code{"gsea_go"} - Gene Ontology Gene Set Enrichment Analysis
#'     \item \code{"kegg"} - KEGG pathway over-representation analysis
#'     \item \code{"gsea_kegg"} - KEGG pathway Gene Set Enrichment Analysis
#'   }
#' @param OrgDb An \code{OrgDb} annotation package object for gene ID conversion and
#'   annotation. Required for all analysis types. Examples: \code{org.Hs.eg.db} for human,
#'   \code{org.Mm.eg.db} for mouse. If \code{NULL}, the function will attempt to detect
#'   the appropriate package but may fail.
#' @param organism A character string specifying the KEGG organism code for KEGG pathway
#'   analysis. Default is \code{"hsa"} for human. Common codes include:
#'   \itemize{
#'     \item \code{"hsa"} - Homo sapiens (human)
#'     \item \code{"mmu"} - Mus musculus (mouse)
#'     \item \code{"rno"} - Rattus norvegicus (rat)
#'     \item \code{"dme"} - Drosophila melanogaster (fruit fly)
#'   }
#'   Only used for KEGG analysis types.
#' @param pvalueCutoff A numeric value specifying the p-value cutoff for enrichment
#'   significance. Default is 1 (no filtering). Common values are 0.05 or 0.01.
#' @param ... Additional arguments (currently not used).
#'
#' @return Returns the input \code{omicscope} object with enrichment results stored
#'   in the \code{enrichmentData} slot. Results can be accessed using
#'   \code{object@enrichmentData[[comparison_name]]}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that differential expression data exists in the object
#'   \item For each selected comparison, extracts gene data and removes duplicates by selecting
#'         the gene with the lowest p-value
#'   \item Converts gene symbols to ENTREZ IDs using the provided OrgDb
#'   \item Performs the specified enrichment analysis:
#'   \itemize{
#'     \item For ORA methods ("go", "kegg"): Analyzes significantly up- and down-regulated genes separately
#'     \item For GSEA methods ("gsea_go", "gsea_kegg"): Uses ranked gene lists based on log2FoldChange
#'   }
#'   \item Stores results in the object's enrichmentData slot with descriptive names
#' }
#'
#' For GO analysis, all three ontologies (BP, CC, MF) are analyzed together (ont = "ALL").
#' For KEGG analysis, results are made readable by converting ENTREZ IDs back to gene symbols.
#'
#' The differential expression data must contain columns: \code{gene_name}, \code{pvalue},
#' \code{log2FoldChange}, and \code{type} (indicating "sigUp" or "sigDown").
#'
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(org.Hs.eg.db)  # for human
#'
#' # Basic GO over-representation analysis
#' obj <- run_enrichment(obj,
#'                       enrich_type = "go",
#'                       OrgDb = org.Hs.eg.db,
#'                       pvalueCutoff = 0.05)
#'
#' # KEGG pathway analysis for specific comparisons
#' obj <- run_enrichment(obj,
#'                       diff_data_selected = c("Treatment_vs_Control"),
#'                       enrich_type = "kegg",
#'                       OrgDb = org.Hs.eg.db,
#'                       organism = "hsa",
#'                       pvalueCutoff = 0.05)
#'
#' # Gene Set Enrichment Analysis for GO
#' obj <- run_enrichment(obj,
#'                       enrich_type = "gsea_go",
#'                       OrgDb = org.Hs.eg.db,
#'                       pvalueCutoff = 0.25)
#'
#' # Mouse analysis
#' library(org.Mm.eg.db)
#' obj <- run_enrichment(obj,
#'                       enrich_type = "gsea_kegg",
#'                       OrgDb = org.Mm.eg.db,
#'                       organism = "mmu")
#'
#' # Access results
#' enrich_results <- obj@enrichmentData
#' go_results <- obj@enrichmentData[["Treatment_vs_Control"]]
#' }
#'
#' @importFrom clusterProfiler bitr gseGO gseKEGG enrichGO enrichKEGG setReadable
#'
#' @export
setGeneric("run_enrichment",function(object,...){
    standardGeneric("run_enrichment")
})








#' @rdname run_enrichment
#' @export
setMethod("run_enrichment",
          signature(object = "omicscope"),
          function(object,
                   diff_data_selected = NULL,
                   enrich_type = c("go","gsea_go","kegg","gsea_kegg"),
                   OrgDb = NULL,
                   organism = "hsa",
                   pvalueCutoff = 1){
              enrich_type <- match.arg(
                  enrich_type,
                  choices = c("go","gsea_go","kegg","gsea_kegg"))
              # ================================================================
              # check diff data
              if(length(object@diffExpData) == 0){
                  stop("Please run differential_expression first!")
              }

              # get diff data
              diff.list <- object@diffExpData

              if(!is.null(diff_data_selected)){
                  diff.list.ft <- diff.list[diff_data_selected]
              }else{
                  diff.list.ft <- diff.list
              }

              group.name <- names(diff.list.ft)
              # ================================================================
              # do enrichment analysis

              # x = 1
              lapply(seq_along(diff.list.ft),function(x){
                  cat(paste0("Start ",group.name[x]," enrichment analysis..."))

                  tmp <- diff.list.ft[[x]]@data |>
                      dplyr::group_by(gene_name) |>
                      dplyr::arrange(pvalue) |>
                      dplyr::slice_head(n = 1) |>
                      dplyr::ungroup()


                  # id transform
                  id <- clusterProfiler::bitr(geneID = tmp$gene_name,
                                              fromType = "SYMBOL",
                                              toType = "ENTREZID",
                                              OrgDb = OrgDb)

                  id.tmp <- id |>
                      dplyr::inner_join(y = tmp,
                                        by = c("SYMBOL" = "gene_name")) |>
                      dplyr::arrange(dplyr::desc(log2FoldChange))

                  # check
                  if(enrich_type %in% c("gsea_go","gsea_kegg")){
                      glist <- id.tmp$log2FoldChange
                      names(glist) <- id.tmp$ENTREZID

                      if(enrich_type == "gsea_go"){
                          ego <- clusterProfiler::gseGO(
                              geneList = glist,
                              ont = "ALL",
                              OrgDb = OrgDb,
                              keyType = "ENTREZID",
                              pvalueCutoff = pvalueCutoff)
                      }else if(enrich_type == "gsea_kegg"){
                          ego <- clusterProfiler::gseKEGG(
                              geneList = glist,
                              organism = organism,
                              pvalueCutoff = pvalueCutoff)

                          # to readble
                          ego <- clusterProfiler::setReadable(
                              x = ego,
                              OrgDb = OrgDb,keyType = "ENTREZID")
                      }

                      golist <- list(ego)

                      # assign names
                      names(golist) <- paste0(group.name[x],"|",enrich_type)

                  }else{
                      # ========================================================
                      # loop for each type
                      tp <- c("sigUp", "sigDown")

                      # j = 1
                      lapply(seq_along(tp),function(j){
                          tmp.data <- subset(id.tmp, type == tp[j])

                          # check enrich type
                          if(enrich_type == "go"){
                              ego <- clusterProfiler::enrichGO(
                                  gene = tmp.data$ENTREZID,
                                  OrgDb = OrgDb,
                                  keyType = "ENTREZID",
                                  ont = "ALL",
                                  pvalueCutoff = pvalueCutoff,
                                  readable = TRUE
                              )

                          }else if(enrich_type == "kegg"){
                              ego <- clusterProfiler::enrichKEGG(
                                  gene = tmp.data$ENTREZID,
                                  organism = organism,
                                  pvalueCutoff = pvalueCutoff)


                              # to readble
                              ego <- clusterProfiler::setReadable(
                                  x = ego,
                                  OrgDb = OrgDb,keyType = "ENTREZID")
                          }
                      }) -> golist

                      # assign names
                      names(golist) <- paste0(group.name[x],"|",
                                              enrich_type,"|",tp)
                  }

                  return(golist)
              }) -> enrich.list

              # assign names
              names(enrich.list) <- group.name

              # ================================================================
              # return
              object@enrichmentData <- enrich.list

              return(object)
          }
)
