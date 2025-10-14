
#' Run Functional Enrichment Analysis
#'
#' @description
#' Performs functional enrichment analysis on differential expression results
#' using Gene Ontology (GO) or KEGG pathway databases. Supports both
#' over-representation analysis (ORA) and gene set enrichment analysis (GSEA).
#'
#' @param object An \code{omicscope} object containing differential expression
#'   results from \code{differential_expression()} method.
#' @param enrich_type Character string specifying the enrichment analysis type.
#'   Options are:
#'   \itemize{
#'     \item \code{"go"} - Gene Ontology over-representation analysis
#'     \item \code{"gsea_go"} - Gene Ontology gene set enrichment analysis
#'     \item \code{"kegg"} - KEGG pathway over-representation analysis
#'     \item \code{"gsea_kegg"} - KEGG pathway gene set enrichment analysis
#'   }
#'   Default is \code{"go"}.
#' @param OrgDb An organism database object (e.g., \code{org.Hs.eg.db} for human,
#'   \code{org.Mm.eg.db} for mouse) used for gene ID conversion and annotation.
#'   Required parameter.
#' @param organism Character string specifying the KEGG organism code
#'   (e.g., "hsa" for human, "mmu" for mouse). Only used for KEGG analysis.
#'   Default is \code{"hsa"}.
#' @param pvalueCutoff Numeric value specifying the p-value cutoff for
#'   significance testing. Default is \code{1} (no filtering).
#' @param enrich_fun_keytype Character string. The `keyType` argument passed to
#'   `clusterProfiler` functions, specifying the type of gene IDs used for enrichment.
#'   Defaults to \code{"ENTREZID"}. Should not typically be changed unless you are using
#'   non-standard gene identifiers.
#' @param id_trans Logical. If \code{TRUE} (default), the function attempts to
#'   translate gene SYMBOLs to ENTREZ IDs. If \code{FALSE}, no translation is performed,
#'   and the `gene_name` column is assumed to already contain the correct IDs.
#' @param kegg_setReadable Logical, applicable only to KEGG analyses (`"kegg"` or `"gsea_kegg"`).
#'   If \code{TRUE} (default), ENTREZ IDs in the results are converted back to readable
#'   gene SYMBOLs.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' This method performs functional enrichment analysis on previously computed
#' differential expression results. The analysis workflow includes:
#'
#' \enumerate{
#'   \item Gene ID conversion from gene symbols to Entrez IDs using \code{clusterProfiler::bitr()}
#'   \item For ORA methods (\code{"go"}, \code{"kegg"}): separate analysis for
#'         up-regulated and down-regulated genes
#'   \item For GSEA methods (\code{"gsea_go"}, \code{"gsea_kegg"}): ranked gene
#'         list analysis using log2 fold change values
#'   \item Results are stored in the \code{enrichmentData} slot of the object
#' }
#'
#' The method requires that differential expression analysis has been performed
#' first using the \code{differential_expression()} method.
#'
#' For GO analysis, all three ontologies (Biological Process, Cellular Component,
#' Molecular Function) are analyzed simultaneously using \code{ont = "ALL"}.
#'
#' @return An \code{omicscope} object with enrichment results stored in the
#'   \code{enrichmentData} slot. The results include:
#'   \itemize{
#'     \item For ORA: separate results for up-regulated and down-regulated genes
#'     \item For GSEA: ranked gene set enrichment results
#'     \item Gene symbols are made readable in the final results
#'   }
#'
#'
#'
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(org.Mm.eg.db)  # for mouse
#' library(org.Hs.eg.db)  # for human
#'
#' # Assuming 'os' is an omicscope object with differential expression results
#'
#' # GO over-representation analysis (mouse)
#' os_go <- run_enrichment(os,
#'                         enrich_type = "go",
#'                         OrgDb = org.Mm.eg.db)
#'
#' # KEGG gene set enrichment analysis (human)
#' os_gsea <- run_enrichment(os,
#'                           enrich_type = "gsea_kegg",
#'                           OrgDb = org.Hs.eg.db,
#'                           organism = "hsa",
#'                           pvalueCutoff = 0.05)
#'
#' # Access enrichment results
#' enrichment_results <- os_go@enrichmentData
#' }
#'
#'
#' @importFrom clusterProfiler bitr enrichGO enrichKEGG gseGO gseKEGG setReadable
#' @importFrom dplyr group_by arrange slice_head ungroup inner_join desc
#'
#'
#' @rdname run_enrichment
#' @export
setGeneric("run_enrichment",function(object,...){
    standardGeneric("run_enrichment")
})








#' @rdname run_enrichment
#' @export
setMethod("run_enrichment",
          signature(object = "omicscope"),
          function(object,
                   enrich_type = c("go","gsea_go","kegg","gsea_kegg"),
                   OrgDb = NULL,
                   organism = "hsa",
                   pvalueCutoff = 1,
                   enrich_fun_keytype = "ENTREZID",
                   id_trans = TRUE,
                   kegg_setReadable = TRUE){
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

              method.name <- names(diff.list)

              # ================================================================
              # do enrichment analysis

              # j = 1
              lapply(seq_along(diff.list),function(j){

                  tmp <- diff.list[[j]]

                  group.name <- names(tmp)

                  cat(paste0("Start for method of ",method.name[j], " and contrast for ",group.name," enrichment analysis..."))

                  # loop for each contrast
                  lapply(seq_along(tmp),function(x){
                      tmp <- tmp[[x]]@data |>
                          dplyr::group_by(gene_name) |>
                          dplyr::arrange(pvalue) |>
                          dplyr::slice_head(n = 1) |>
                          dplyr::ungroup()


                      # id transform
                      if(id_trans == TRUE){
                          id <- clusterProfiler::bitr(geneID = tmp$gene_name,
                                                      fromType = "SYMBOL",
                                                      toType = "ENTREZID",
                                                      OrgDb = OrgDb)

                          id.tmp <- id |>
                              dplyr::inner_join(y = tmp,
                                                by = c("SYMBOL" = "gene_name")) |>
                              dplyr::arrange(dplyr::desc(log2FoldChange))
                      }else{
                          id.tmp <-  tmp |>
                              dplyr::mutate(ENTREZID = gene_name)
                      }


                      # check
                      if(enrich_type %in% c("gsea_go","gsea_kegg")){
                          glist <- id.tmp$log2FoldChange
                          names(glist) <- id.tmp$ENTREZID

                          if(enrich_type == "gsea_go"){
                              ego <- clusterProfiler::gseGO(
                                  geneList = glist,
                                  ont = "ALL",
                                  OrgDb = OrgDb,
                                  keyType = enrich_fun_keytype,
                                  pvalueCutoff = pvalueCutoff)
                          }else if(enrich_type == "gsea_kegg"){
                              ego <- clusterProfiler::gseKEGG(
                                  geneList = glist,
                                  organism = organism,
                                  pvalueCutoff = pvalueCutoff)

                              # to readble
                              if(kegg_setReadable == TRUE){
                                  ego <- clusterProfiler::setReadable(
                                      x = ego,
                                      OrgDb = OrgDb,keyType = enrich_fun_keytype)
                              }

                          }

                          golist <- list(ego)

                          # assign names
                          names(golist) <- paste0(method.name[j],"|",group.name[x],"|",enrich_type)

                      }else{
                          # ========================================================
                          # loop for each type
                          tp <- c("sigUp", "sigDown")

                          # z = 1
                          lapply(seq_along(tp),function(z){
                              tmp.data <- subset(id.tmp, type == tp[z])

                              # check enrich type
                              if(enrich_type == "go"){
                                  ego <- clusterProfiler::enrichGO(
                                      gene = tmp.data$ENTREZID,
                                      OrgDb = OrgDb,
                                      keyType = enrich_fun_keytype,
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
                                  if(kegg_setReadable == TRUE){
                                      ego <- clusterProfiler::setReadable(
                                          x = ego,
                                          OrgDb = OrgDb,keyType = enrich_fun_keytype)
                                  }

                              }
                          }) -> golist

                          # assign names
                          names(golist) <- paste0(method.name[j],"|",
                                                  group.name[x],"|",
                                                  enrich_type,"|",tp)
                      }

                      return(golist)
                  }) -> reslist

                  return(reslist)
              }) -> enrich.list


              # ================================================================
              # return
              if(enrich_type == "go"){
                  object@enrichmentData[["GO"]] <- S4Vectors::SimpleList(enrich.list[[1]][[1]])
              }else if(enrich_type == "gsea_go"){
                  object@enrichmentData[["GSEA_GO"]] <- S4Vectors::SimpleList(enrich.list[[1]][[1]])
              }else if(enrich_type == "kegg"){
                  object@enrichmentData[["KEGG"]] <- S4Vectors::SimpleList(enrich.list[[1]][[1]])
              }else if(enrich_type == "gsea_kegg"){
                  object@enrichmentData[["GSEA_KEGG"]] <- S4Vectors::SimpleList(enrich.list[[1]][[1]])
              }

              return(object)
          }
)
