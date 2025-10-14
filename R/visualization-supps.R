
#' @title Generate Coverage Profile Plots with Gene Models
#' @description
#' This function calculates read coverage (in Reads Per Million, RPM) for specified
#' gene regions from BAM files and visualizes it alongside the corresponding gene/transcript
#' structures derived from a GTF file. It uses `ggh4x` for flexible faceting to
#' compare profiles across different samples/groups and genes.
#'
#' @param bam_file A character vector of paths to one or more BAM files.
#' @param sample_name A character vector providing sample names for each BAM file.
#'   If NULL, the BAM filenames are used.
#' @param group_name A character vector assigning a group to each sample.
#'   If NULL, the BAM filenames are used.
#' @param gtf_file Path to a GTF file or a `GRanges` object containing gene annotations.
#' @param target_gene A character vector of one or more target gene names (e.g., from the `gene_name`
#' attribute in the GTF) to visualize.
#' @param sample_order A character vector to specify the display order of samples in the facets.
#' @param group_order A character vector to specify the display order of groups in the facets.
#' @param gene_order A character vector to specify the display order of genes in the facets.
#' @param merge_group Logical. If `TRUE`, RPM values for samples within the same group are
#'   averaged and displayed as a single group track. Defaults to `FALSE`.
#' @param add_range_label Logical. If `TRUE`, adds a text label indicating the Y-axis
#'   range to each coverage panel. Defaults to `FALSE`.
#' @param range_digit Numeric. The number of decimal places to use for the range label.
#'   Defaults to 1.
#' @param remove_labelY Logical. If `TRUE`, removes Y-axis tick labels and marks.
#'   Defaults to `FALSE`.
#' @param collapse_exon Logical. If `TRUE`, all transcripts of a gene are collapsed
#'   into a single, flattened gene model for display. Defaults to `FALSE`.
#' @param arrow_length Numeric. Controls the size of the arrowhead for the transcript
#'   direction arrows. Defaults to 1.2.
#' @param arrow_linewidth Numeric. Controls the line width of the transcript direction
#'   arrows. Defaults to 0.3.
#' @param arrow_col Character string. The color of the transcript direction arrows.
#'   Defaults to "black".
#' @param range_pos A numeric vector of length 2 `c(x, y)` specifying the relative position
#'   (from 0 to 1) of the range label within the panel. Defaults to `c(0.9, 0.9)`.
#' @param exon_col Character string. The color of the exon segments in the gene model.
#'   Defaults to "black".
#' @param exon_linewidth Numeric. The line width of the exon segments. Defaults to 3.
#'
#' @return A `ggplot` object, which can be printed to display the plot or modified
#'   further with additional ggplot2 layers.
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges GRangesList reduce seqnames start
#' @importFrom Rsamtools indexBam idxstatsBam pileup BamFile PileupParam ScanBamParam
#' @importFrom dplyr filter select mutate group_by summarise arrange row_number ungroup left_join bind_rows distinct n reframe first
#' @importFrom ggplot2 ggplot aes geom_col geom_segment geom_blank geom_text scale_y_continuous scale_x_continuous theme_bw theme element_blank element_text xlab ylab expansion
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This example requires actual BAM and GTF files to run.
#' # The following shows the basic function call structure.
#'
#' # Assume you have the following files:
#' my_bams <- c("path/to/sample1.bam", "path/to/sample2.bam")
#' my_gtf <- "path/to/genes.gtf"
#' target_genes <- c("Actb", "Nanog")
#'
#' # Basic usage
#' p <- coverage_plot(
#'   bam_file = my_bams,
#'   sample_name = c("Sample1", "Sample2"),
#'   group_name = c("Control", "Treatment"),
#'   gtf_file = my_gtf,
#'   target_gene = target_genes
#' )
#' print(p)
#'
#' # Customized plot with merged groups
#' p_merged <- coverage_plot(
#'   bam_file = c("path/to/ctrl_rep1.bam", "path/to/ctrl_rep2.bam"),
#'   sample_name = c("Control_rep1", "Control_rep2"),
#'   group_name = c("Control", "Control"),
#'   gtf_file = my_gtf,
#'   target_gene = target_genes,
#'   merge_group = TRUE,
#'   add_range_label = TRUE,
#'   collapse_exon = TRUE,
#'   arrow_col = "blue",
#'   exon_col = "darkblue"
#' )
#' print(p_merged)
#' }
coverage_plot <- function(bam_file = NULL,
                          sample_name = NULL,
                          group_name = NULL,
                          gtf_file = NULL,
                          target_gene = NULL,
                          sample_order = NULL, group_order = NULL, gene_order = NULL,
                          merge_group = FALSE,
                          add_range_label = FALSE, range_digit = 1,
                          remove_labelY = FALSE,
                          collapse_exon = FALSE,
                          arrow_length = 1.2, arrow_linewidth = 0.3, arrow_col = "black",
                          range_pos = c(0.9, 0.9),
                          exon_col = "black", exon_linewidth = 3){
    # ==========================================================================
    # check package needed
    if (!requireNamespace(c("ggh4x","ggarrow"), quietly = TRUE)) {
        stop("Package 'ggh4x' and 'ggarrow' are required for coverage track plot.\n",
             "Please install them with: install.packages('ggh4x') \n
             install.packages('ggarrow')")
    }
    # ==========================================================================
    # get gene coordinate
    if(inherits(gtf_file,"GRanges")){
        gtf <- data.frame(gtf_file, check.names = FALSE)
    }else if(is.character(gtf_file)){
        gtf <- rtracklayer::import(gtf_file) |>
            data.frame(check.names = FALSE)
    }


    # filter exons
    exon <- subset(gtf, type == "exon")

    # filter genes
    # x = 1
    lapply(seq_along(target_gene),function(x){
        query_region <- exon |>
            dplyr::filter(gene_name == target_gene[x]) |>
            GenomicRanges::GRanges() |>
            GenomicRanges::reduce()

        query_region$gene_name <- target_gene[x]

        return(query_region)
    }) |> GenomicRanges::GRangesList() |>
        unlist(,use.names = FALSE) -> query_region

    # ==========================================================================
    # get coverage
    genes <- unique(query_region$gene_name)

    # check sample_name
    if(is.null(sample_name)){
        sample_name <- bam_file
    }

    if(is.null(group_name)){
        group_name <- bam_file
    }

    # x = 1
    # g = 1
    lapply(seq_along(bam_file),function(x){
        lapply(seq_along(genes),function(g){
            fts <- subset(query_region, gene_name == genes[g])

            # check index file for bam
            if(!file.exists(paste(bam_file[x],"bai",sep = "."))){
                Rsamtools::indexBam(bam_file[x], nThreads = parallel::detectCores())
            }

            total_reads <- sum(Rsamtools::idxstatsBam(bam_file[x])$mapped)

            pileup_result <- Rsamtools::pileup(
                file = Rsamtools::BamFile(bam_file[x]),
                pileupParam = Rsamtools::PileupParam(distinguish_nucleotides = FALSE,
                                                     distinguish_strands = FALSE,
                                                     max_depth = 10^8),
                scanBamParam = Rsamtools::ScanBamParam(which = fts)) |>
                dplyr::select(-which_label)

            # check data
            if(nrow(pileup_result) == 0){
                pileup_result <- data.frame(seqnames = GenomicRanges::seqnames(fts) |> unique(),
                                            pos = GenomicRanges::start(fts)[1],
                                            count = 0)
            }

            pileup_result$sample <- sample_name[x]
            pileup_result$group <- group_name[x]
            pileup_result$gene_name <- genes[g]
            pileup_result$rpm <- (pileup_result$count/total_reads)*10^6

            return(pileup_result)
        }) %>% do.call("rbind",.) |> data.frame(check.names = FALSE) -> res

        return(res)
    }) %>% do.call("rbind",.) |> data.frame(check.names = FALSE) -> cov.res

    # average replicates
    if(merge_group == TRUE){
        cov.res <- cov.res |>
            dplyr::group_by(group,gene_name,seqnames,pos) %>%
            dplyr::summarise(rpm = mean(rpm)) |>
            dplyr::mutate(sample = group)

        facet_layer <- ggh4x::facet_grid2(group ~ gene_name, scales = "free", independent = "y", switch = "y")
    }else{
        facet_layer <- ggh4x::facet_grid2(sample ~ gene_name, scales = "free", independent = "y", switch = "y")
    }

    # add max value
    cov.res <- cov.res |>
        dplyr::group_by(gene_name) |>
        dplyr::mutate(y_max = max(rpm, na.rm = TRUE))

    # ==========================================================================
    # prepare plot data
    strc <- exon |>
        dplyr::filter(gene_name %in% target_gene) |>
        dplyr::select(seqnames,start,end,gene_name,transcript_id) |>
        dplyr::mutate(sample = "Gene structure",group = "Gene structure")

    tid <- gtf |>
        dplyr::filter(gene_name %in% target_genes) |>
        dplyr::filter(type == "transcript")

    strc.info <- strc |>
        dplyr::mutate(exon_length = end - start) |>
        dplyr::group_by(gene_name, transcript_id) |>
        dplyr::summarise(transcript_length = sum(exon_length), .groups = 'drop') |>
        dplyr::arrange(gene_name, transcript_length) |>
        dplyr::group_by(gene_name) |>
        dplyr::mutate(transcript_rank = dplyr::row_number()) |>
        dplyr::ungroup() |>
        dplyr::group_by(gene_name) |>
        dplyr::mutate(y_max = dplyr::n() + 1)

    if(collapse_exon == TRUE){
        strc.info$y_max <- 2
    }

    strc <- strc |>
        dplyr::left_join(y = strc.info,by = c("gene_name","transcript_id"))

    # add arrow direction
    tid <- tid |>
        left_join(y = strc.info,by = c("gene_name","transcript_id")) |>
        dplyr::mutate(sample = "Gene structure",group = "Gene structure") |>
        dplyr::mutate(x = ifelse(strand == "+", start, end),
                      xend = ifelse(strand == "+", end, start))

    cov.strc <- dplyr::bind_rows(cov.res,strc)

    # rpm range text label
    rg <- cov.res |>
        dplyr::select(sample, group, gene_name, y_max) |>
        dplyr::distinct() |>
        dplyr::mutate(label = paste0("[0-",round(y_max,digits = range_digit),"]"))

    # ==========================================================================
    # orders
    if(!is.null(sample_order)){
        cov.strc$sample <- factor(cov.strc$sample, levels = c(sample_order,"Gene structure"))
        rg$sample <- factor(rg$sample, levels = sample_order)
        tid$sample <- factor(tid$sample, levels = c(sample_order,"Gene structure"))
        strc$sample <- factor(strc$sample, levels = c(sample_order,"Gene structure"))
    }

    if(!is.null(group_order)){
        cov.strc$group <- factor(cov.strc$group, levels = c(group_order,"Gene structure"))
        rg$group <- factor(rg$group, levels = group_order)
        tid$group <- factor(tid$group, levels = c(sample_order,"Gene structure"))
        strc$group <- factor(strc$group, levels = c(sample_order,"Gene structure"))
    }

    if(!is.null(gene_order)){
        cov.strc$gene_name <- factor(cov.strc$gene_name, levels = gene_order)
        rg$gene_name <- factor(rg$gene_name, levels = gene_order)
        tid$gene_name <- factor(tid$gene_name, levels = gene_order)
        strc$gene_name <- factor(strc$gene_name, levels = gene_order)
    }
    # ==========================================================================
    # plot
    p <-
        ggplot(cov.strc) +
        geom_col(aes(x = pos, y = rpm,fill = gene_name,color = gene_name),
                 show.legend = FALSE)

    if(collapse_exon == TRUE){
        tid2 <- tid |>
            dplyr::group_by(gene_name,sample,group) |>
            dplyr::reframe(x = ifelse(strand == "+",min(x),max(x)),
                           xend = ifelse(strand == "+",max(xend),min(xend)))

        p <- p +
            # arrow layer
            ggarrow::geom_arrow_segment(data = tid2,
                                        aes(x = x, xend = xend,y = 1, yend = 1),
                                        color = arrow_col,
                                        arrow_mid = ggarrow::arrow_head_line(),
                                        mid_place = seq(0,1,0.1),
                                        length = arrow_length,
                                        linewidth = arrow_linewidth) +
            # exon layer
            geom_segment(data = strc,
                         aes(x = start, xend = end,y = 1, yend = 1),
                         color = exon_col,
                         linewidth = exon_linewidth)
    }else{
        p <- p +
            # arrow layer
            ggarrow::geom_arrow_segment(data = tid,
                                        aes(x = x, xend = xend,y = transcript_rank, yend = transcript_rank),
                                        color = arrow_col,
                                        arrow_mid = ggarrow::arrow_head_line(),
                                        mid_place = seq(0,1,0.1),
                                        length = arrow_length,
                                        linewidth = arrow_linewidth) +
            # exon layer
            geom_segment(data = strc,
                         aes(x = start, xend = end,y = transcript_rank, yend = transcript_rank),
                         color = exon_col,
                         linewidth = exon_linewidth)
    }


    p <- p +
        facet_layer +
        geom_blank(aes(y = y_max)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                           position = "right") +
        scale_x_continuous(labels = function(x) x/1000) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              strip.text.x = element_text(face = "bold.italic"),
              strip.background.y = element_blank(),
              strip.text.y = element_text(face = "bold",hjust = 1),
              strip.text.y.left = element_text(angle = 0),
              axis.text = element_text(colour = "black")) +
        xlab("Genomic coordinate (kb)") + ylab("Reads coverage (rpm)")

    # check range label
    if(add_range_label == TRUE){
        p <- p +
            # range layer
            geom_text(data = rg,
                      aes(x = I(range_pos[1]),y = I(range_pos[2]),label = label))
    }

    # remove y label
    if(remove_labelY == TRUE){
        p <- p +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())
    }

    return(p)
}
