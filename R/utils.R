#' Calculate the Reduced (Non-redundant) Length of Genes
#'
#' @description
#' This function calculates the total non-redundant length of all exons
#' associated with each gene. It first groups exons by gene, then merges any
#' overlapping or adjacent exons into a single feature, and finally sums the
#' widths of these merged features to get a single length value per gene.
#'
#'
#' @param gtf A \code{GRanges} object, typically imported from a GTF or GFF file.
#'   It must contain metadata columns named 'type' (to filter for "exon") and
#'   'gene_id' (for grouping).
#'
#' @return A \code{data.frame} with two columns:
#'   \item{gene_id}{Character, the unique gene identifier.}
#'   \item{length}{Integer, the calculated non-redundant exonic length for that gene.}
#'
#' @importFrom GenomicRanges reduce width
#' @importFrom S4Vectors split
#' @export
reduced_gene_length <- function(gtf = NULL){
    # Filter for exons
    exons <- subset(gtf, type == "exon")

    # Group exons by gene ID
    exons_by_gene <- split(exons, ~gene_id)

    # Merge overlapping exons for each gene
    reduced_exons_list <- GenomicRanges::reduce(exons_by_gene, ignore.strand = T)

    # Sum the lengths of the merged exons for each gene
    nonredundant_lengths <- sum(GenomicRanges::width(reduced_exons_list))

    # Create a final data frame with gene lengths
    gene_lengths <- data.frame(gene_id = names(nonredundant_lengths),
                               length = unlist(nonredundant_lengths),
                               row.names = NULL)

    return(gene_lengths)
}
