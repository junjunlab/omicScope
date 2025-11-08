


#' Plot Single-Cell Coverage Tracks with Gene Annotations
#'
#' @description
#' This function generates a composite plot displaying single-cell coverage tracks
#' (e.g., from scATAC-seq) for specified genomic regions or genes.
#' It integrates gene/transcript structures from a GTF file and can optionally
#' overlay other genomic features like peaks and interaction links. The plot is
#' highly customizable, using `ggplot2` and `ggh4x` for faceting.
#'
#' @details
#' The function fetches coverage data using an internal helper function
#' (conceptually similar to `Signac::CoveragePlot`) and processes a GTF file to
#' draw gene models. It aligns these tracks, allowing for comparison across
#' different samples, groups, or genomic loci.
#'
#' Either `target_gene` or `target_region` must be specified, but not both.
#'
#' @param object A single-cell object from which to pull coverage data (e.g., a Seurat object
#'   with an Assay containing fragment file information).
#' @param gtf_file A `GRanges` object or a character string specifying the path to a
#'   GTF/GFF file for gene annotations.
#' @param links_data An optional `data.frame` for plotting genomic interaction links
#'   (e.g., from Hi-C or Cicero). It must contain at least `seqnames`, `start`, `end`,
#'   and `sample` columns. **Crucially, if plotting is faceted by gene or region
#'   (i.e., using the `target_gene` or `target_region` arguments), this data frame
#'   must also include a corresponding column (`gene_name` or `target_region`)
#'   with values matching the targets to ensure links are mapped to the correct plot panel.**
#' @param peaks_data An optional `data.frame` for plotting genomic peaks (e.g., from
#'   MACS2). If provided, it must contain at least the following columns: `seqnames`,
#'   `start`, `end`, `sample_name`, and `sample`. The `sample` column is used to
#'   group all peaks into a single dedicated track and its value must be set to `"Peaks"`.
#'   **Crucially, if plotting is faceted by gene or region (using `target_gene` or
#'   `target_region`), this data frame must also include a corresponding column
#'   (`gene_name` or `target_region`) with values matching the targets to ensure
#'   peaks are mapped to the correct plot panel.** If `NULL` (the default), peaks are
#'   fetched automatically from the provided single-cell `object`.
#' @param peaks_width Numeric. The linewidth aesthetics for drawing peaks. Default is 3.
#' @param target_gene A character vector of gene names to plot. The function will
#'   determine the genomic coordinates automatically from the `gtf_file`.
#' @param target_region A character vector of genomic regions to plot, formatted as
#'   "chr:start-end" (e.g., "chr1:1000-2000").
#' @param SingleCoveragePlot2_params A list of additional parameters passed to the
#'   internal coverage fetching function (e.g., `Signac::SingleCoveragePlot`).
#' @param sample_col A `ggplot2` scale object (e.g., `scale_fill_manual`) for coloring
#'   the coverage tracks. Defaults to `scale_fill_discrete()`.
#' @param peaks_sample_col A `ggplot2` scale object (e.g., `scale_color_manual`) for
#'   coloring the peak tracks by sample. Defaults to `scale_color_discrete()`.
#' @param links_col A `ggplot2` scale object (e.g., `scale_color_gradient`) for
#'   coloring the interaction links. Defaults to a gradient.
#' @param extend_up Numeric. Number of base pairs to extend the plotting region
#'   upstream of the gene/region bounds. Default is 2000.
#' @param extend_down Numeric. Number of base pairs to extend the plotting region
#'   downstream. Default is 2000.
#' @param sample_order,group_order,gene_order,region_order,links_sample_order Character
#'   vectors specifying the desired order of items in the plot facets. If `NULL`,
#'   default ordering is used.
#' @param links_col_var Character. The name of the column in `links_data` to use for
#'   coloring the links (e.g., "score", "pvalue"). Default is "score".
#' @param merge_group Logical. If `TRUE`, coverage signals from samples within the same
#'   group are averaged, and one track per group is shown. Default is `FALSE`.
#' @param add_range_label Logical. If `TRUE`, adds a text label indicating the
#'   y-axis range (e.g., `[0-10]`) for each coverage track. Default is `FALSE`.
#' @param range_digit Integer. The number of decimal places to round the range label to.
#'   Default is 0.
#' @param remove_labelY Logical. If `TRUE`, remove y-axis text and ticks for a cleaner plot.
#'   Default is `FALSE`.
#' @param collapse_exon Logical. If `TRUE`, all transcripts for a given gene are
#'   collapsed into a single, representative gene model. Default is `FALSE`.
#' @param show_utr Logical. If `TRUE`, UTRs and CDS regions are drawn with different
#'   thicknesses. If `FALSE`, all exonic parts are drawn with the same thickness.
#'   Default is `TRUE`.
#' @param arrow_length,arrow_linewidth,arrow_col Aesthetics for the arrows on the
#'   gene models indicating transcription direction. Passed to `ggarrow::geom_arrow_segment`.
#' @param range_pos A numeric vector of length 2, `c(x, y)`, specifying the position of
#'   the range label in normalized plot coordinates (from 0 to 1). Default is `c(0.9, 0.9)`.
#' @param exon_col,exon_linewidth Aesthetics for the exons in the gene models.
#' @param add_backsqure Logical. If `TRUE`, use `geom_label` for gene names (with a
#'   background box); if `FALSE`, use `geom_text`. Default is `TRUE`.
#' @param add_gene_label Logical. If `TRUE`, adds gene name labels to the gene structure
#'   track when using `target_region`. Defaults to `TRUE`.
#' @param gene_label_size Numeric. The size of the gene name labels. Defaults to 1.
#' @param gene_label_aes Character. Column name to use from GTF for gene labeling (like "gene_name", "gene_id").
#' Defaults to "gene_name".
#' @param highlight_region An optional `data.frame` with `seqnames`, `start`, `end`
#'   columns to add a shaded background highlighting specific genomic regions.
#' @param highlight_col_aes Character. The column name in `highlight_region` to map to the
#'   fill color. Default is "seqnames".
#' @param highlight_col A named vector or `ggplot2` scale to control the fill colors
#'   of the highlighted regions.
#' @param highlight_alpha Numeric. The alpha transparency for the highlighted regions.
#'   Default is 0.2.
#'
#' @return A `ggplot` object.
#'
#' @importFrom dplyr filter select mutate group_by summarise left_join arrange row_number n ungroup distinct reframe cur_group_id bind_rows
#' @importFrom GenomicRanges GRanges GRangesList reduce seqnames
#' @importFrom IRanges start end
#' @importFrom ggplot2 ggplot aes geom_area geom_segment geom_blank scale_y_continuous expansion scale_x_continuous theme_bw theme element_blank element_text xlab ylab geom_text geom_label scale_fill_discrete scale_color_discrete scale_color_gradient scale_fill_manual
#' @importFrom ggbio geom_arch
#' @importFrom scales rescale
#' @importFrom utils modifyList
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'pbmc' is a Seurat object with a scATAC-seq assay, and
#' # 'hm' is a GRanges object loaded from a GTF file.
#' hm <- rtracklayer::import.gff("gencode.v49.annotation.gtf.gz")
#'
#' # Example 1: Plot coverage for a single gene
#' scCoverage_plot(object = pbmc,
#'                 gtf_file = hm,
#'                 target_gene = "CD8A")
#'
#' # Example 2: Plot coverage for a specific genomic region
#' scCoverage_plot(object = pbmc,
#'                 gtf_file = hm,
#'                 target_region = "chr2:86784605-86808396")
#'
#' # Example 3: Plot multiple genes, collapse exons, and customize labels
#' scCoverage_plot(object = pbmc,
#'                 gtf_file = hm,
#'                 target_gene = c('MS4A1', 'CD3D', 'LEF1'),
#'                 add_range_label = TRUE,
#'                 remove_labelY = TRUE,
#'                 range_pos = c(0.9, 0.8),
#'                 collapse_exon = TRUE) +
#'   theme(panel.spacing.y = unit(0, "mm"))
#'
#' # Example 4: Add custom peaks and interaction links data
#' # First, create some dummy data for demonstration
#' gfl <- dplyr::filter(as.data.frame(hm), gene_name == "CD8A" & type == "gene")
#' links <- data.frame(
#'   seqnames = gfl$seqnames[1],
#'   start = sample(gfl$start:gfl$end, 20, replace = TRUE),
#'   end = sample(gfl$start:gfl$end, 20, replace = TRUE),
#'   score = runif(20),
#'   sample = rep(c("links-1", "links-2"), 10)
#' )
#' links <- dplyr::filter(links, end > start)
#'
#' peaks <- data.frame(
#'   seqnames = gfl$seqnames[1],
#'   start = sample(gfl$start:gfl$end, 20, replace = TRUE),
#'   end = sample(gfl$start:gfl$end, 20, replace = TRUE),
#'   sample_name = rep(c("peaks-1", "peaks-2"), 10),
#'   sample = "Peaks" # A general category for the track
#' )
#' peaks <- dplyr::filter(peaks, end > start)
#'
#' scCoverage_plot(object = pbmc,
#'                 gtf_file = hm,
#'                 target_gene = "CD8A",
#'                 links_data = links,
#'                 peaks_data = peaks)
#' }
#'
#'
scCoverage_plot <- function(object = NULL,
                            gtf_file = NULL,
                            links_data = NULL,
                            peaks_data = NULL,
                            peaks_width = 3,
                            target_gene = NULL,
                            target_region = NULL,
                            SingleCoveragePlot2_params = list(),
                            sample_col = scale_fill_discrete(),
                            peaks_sample_col = scale_color_discrete(),
                            links_col = scale_color_gradient(),
                            extend_up = 2000,extend_down = 2000,
                            sample_order = NULL, group_order = NULL,
                            gene_order = NULL, region_order = NULL,
                            links_sample_order = NULL,
                            links_col_var = "score",
                            merge_group = FALSE,
                            add_range_label = FALSE, range_digit = 0,
                            remove_labelY = FALSE,
                            collapse_exon = FALSE, show_utr = TRUE,
                            arrow_length = 1.2, arrow_linewidth = 0.3, arrow_col = "black",
                            range_pos = c(0.9, 0.9),
                            exon_col = "black", exon_linewidth = 3,
                            add_backsqure = TRUE,
                            add_gene_label = TRUE, gene_label_size = 1, gene_label_aes = "gene_name",
                            highlight_region = NULL,
                            highlight_col_aes = "seqnames", highlight_col = NULL, highlight_alpha = 0.2){
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

    # check gene names
    if(!is.null(target_gene)){
        gin <- target_gene %in% unique(gtf$gene_name)

        if("FALSE" %in% gin){
            message(paste0("Gene name: ",target_gene[!gin]," can't be found in gtf file, please check its symbol. This gene will not be used for plot."))
        }
    }

    # filter exons
    exon <- subset(gtf, type == "exon")

    # filter genes
    # x = 1
    if(!is.null(target_gene)){
        var <- "gene_name"

        lapply(seq_along(target_gene),function(x){
            query_region <- exon |>
                dplyr::filter(gene_name == target_gene[x]) |>
                GenomicRanges::GRanges() |>
                GenomicRanges::reduce()

            query_region$gene_name <- target_gene[x]

            return(query_region)
        }) |> GenomicRanges::GRangesList() |>
            unlist(,use.names = FALSE) -> query_region
    }else if(!is.null(target_region)){
        var <- "target_region"

        # x = 1
        lapply(seq_along(target_region),function(x){
            chr <- sapply(strsplit(target_region[x],split = ":"), "[",1)
            tmp <- sapply(strsplit(target_region[x],split = ":"), "[",2)
            st <- sapply(strsplit(tmp,split = "-"), "[",1) |> as.numeric()
            ed <- sapply(strsplit(tmp,split = "-"), "[",2) |> as.numeric()


            query_region <- exon |>
                dplyr::filter(seqnames == chr & start >= st & end <= ed) |>
                GenomicRanges::GRanges() |>
                GenomicRanges::reduce()

            query_region$target_region <- target_region[x]

            return(query_region)
        }) |> GenomicRanges::GRangesList() |>
            unlist(,use.names = FALSE) -> query_region
    }


    # ==========================================================================
    # get coverage
    if(!is.null(target_region)){
        genes <- unique(query_region$target_region)
    }else{
        genes <- unique(query_region$gene_name)
    }

    if(!is.null(target_gene)){
        region <- target_gene
    }else{
        region <- target_region
    }

    # loop to get coverage for clusters
    # x = 1
    lapply(seq_along(region),function(x){
        if(!is.null(target_region)){
            fts <- subset(query_region, target_region == genes[x])
        }else{
            fts <- subset(query_region, gene_name == genes[x])
        }

        cov <- do.call(SingleCoveragePlot2,
                       modifyList(list(
                           object = object,
                           region = region[x],
                           extend.upstream = extend_up,
                           extend.downstream = extend_down),
                           SingleCoveragePlot2_params))[[1]]

        if(!is.null(target_region)){
            cov$target_region <- genes[x]
        }else{
            cov$gene_name <- genes[x]
        }

        cov$seqnames <- as.character(GenomicRanges::seqnames(fts)[1])

        return(cov)
    }) %>% do.call("rbind",.) %>%
        data.frame(check.names = FALSE) -> cov.res

    # loop to get peaks
    # x = 1
    lapply(seq_along(region),function(x){
        if(!is.null(target_region)){
            fts <- subset(query_region, target_region == genes[x])
        }else{
            fts <- subset(query_region, gene_name == genes[x])
        }

        pk <- do.call(SingleCoveragePlot2,
                      modifyList(list(
                          object = object,
                          region = region[x],
                          extend.upstream = extend_up,
                          extend.downstream = extend_down),
                          SingleCoveragePlot2_params))[[2]]

        pk <- pk[,1:3]

        if(nrow(pk) == 0){
            pk <- data.frame(seqnames = "",
                             start = IRanges::start(fts)[1],
                             end = IRanges::start(fts)[1])
        }

        pk$sample <- "Peaks"
        pk$sample_name <- "Peaks"
        pk$group <- "Peaks"
        pk$y <- 1
        pk$y_max <- 2

        if(!is.null(target_region)){
            pk$target_region <- genes[x]
        }else{
            pk$gene_name <- genes[x]
        }

        pk$seqnames <- as.character(GenomicRanges::seqnames(fts)[1])

        return(pk)
    }) %>% do.call("rbind",.) %>%
        data.frame(check.names = FALSE) -> peak.res

    ylb <- "Normalized signal"


    # check prefix for seqname
    bm <- startsWith(as.character(cov.res$seqnames[1]),"chr")
    anno <- startsWith(as.character(gtf$seqnames[1]),"chr")

    if(!bm == anno){
        stop("The prefix of chromosome for gtf and bam file are not same! Please check!")
    }

    # average replicates
    if(merge_group == TRUE){
        cov.res <- cov.res |>
            dplyr::group_by(group,.data[[var]],seqnames,pos) %>%
            dplyr::summarise(rpm = mean(rpm)) |>
            dplyr::mutate(sample = group)

        facet_layer <- ggh4x::facet_grid2(group ~ .data[[var]], scales = "free", independent = "y", switch = "y")
    }else{
        facet_layer <- ggh4x::facet_grid2(sample ~ .data[[var]], scales = "free", independent = "y", switch = "y")
    }

    # add max value
    cov.res <- cov.res |>
        dplyr::group_by(.data[[var]]) |>
        dplyr::mutate(y_max = max(rpm, na.rm = TRUE))

    # ==========================================================================
    # prepare plot data
    if(is.null(target_region)){

        if(show_utr == TRUE){
            strc <- subset(gtf, type %in% c("CDS","five_prime_utr", "three_prime_utr","5UTR","3UTR"))
        }else{
            strc <- exon
        }


        strc <- strc |>
            dplyr::filter(gene_name %in% target_gene) |>
            dplyr::select(seqnames,start,end,gene_name,transcript_id,type) |>
            dplyr::mutate(sample = "Gene structure",group = "Gene structure")

        tid <- gtf |>
            dplyr::filter(gene_name %in% target_gene) |>
            dplyr::filter(type == "transcript")

        # extend region
        gid <- gtf |>
            dplyr::filter(gene_name %in% target_gene) |>
            dplyr::filter(type == "gene") |>
            dplyr::group_by(gene_name) |>
            dplyr::reframe(x_pos = c(min(start) - extend_up, max(end) + extend_down))


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

    }else{
        # loop
        lapply(seq_along(target_region),function(x){
            chr <- sapply(strsplit(target_region[x],split = ":"), "[",1)
            tmp <- sapply(strsplit(target_region[x],split = ":"), "[",2)
            st <- sapply(strsplit(tmp,split = "-"), "[",1) |> as.numeric()
            ed <- sapply(strsplit(tmp,split = "-"), "[",2) |> as.numeric()

            if(show_utr == TRUE){
                strc <- subset(gtf, type %in% c("CDS","five_prime_utr", "three_prime_utr","5UTR","3UTR"))
            }else{
                strc <- exon
            }

            strc <- strc |>
                dplyr::filter(seqnames == chr & start >= st & end <= ed) |>
                dplyr::select(seqnames,start,end,gene_name,transcript_id,type) |>
                dplyr::mutate(sample = "Gene structure",group = "Gene structure")
            strc$target_region <- target_region[x]

            tid <- gtf |>
                dplyr::filter(transcript_id %in% unique(strc$transcript_id)) |>
                dplyr::filter(type == "transcript")
            tid$target_region <- target_region[x]

            strc.info <- get_trans_pos(strc)
            strc.info$target_region <- target_region[x]
            strc.info$y_max <- max(strc.info$transcript_rank + 1)

            return(list(strc, tid, strc.info))
        }) -> tmp.lst


        strc <- do.call(rbind, lapply(tmp.lst, function(x) x[[1]]))
        tid <- do.call(rbind, lapply(tmp.lst, function(x) x[[2]]))
        strc.info <- do.call(rbind, lapply(tmp.lst, function(x) x[[3]]))

        if(collapse_exon == TRUE){
            strc.info$y_max <- 2
        }

        strc <- strc |>
            dplyr::left_join(y = strc.info,by = c(var,"gene_name","transcript_id"))

        # add arrow direction
        tid <- tid |>
            left_join(y = strc.info,by = c(var,"gene_name","transcript_id")) |>
            dplyr::mutate(sample = "Gene structure",group = "Gene structure") |>
            dplyr::mutate(x = ifelse(strand == "+", start, end),
                          xend = ifelse(strand == "+", end, start))
    }


    cov.strc <- dplyr::bind_rows(cov.res,strc)


    # rpm range text label
    rg <- cov.res |>
        dplyr::select(sample, group, .data[[var]], y_max) |>
        dplyr::distinct() |>
        dplyr::mutate(label = paste0("[0-",round(y_max,digits = range_digit),"]"))

    # ==========================================================================
    # orders
    if(!is.null(sample_order)){
        lvs <- sample_order
    }else{
        lvs <- levels(cov.res$sample)
    }

    if(is.null(peaks_data)){
        lvs <-  c(lvs,unique(peak.res$sample))
    }else{
        peak.res <- peaks_data
        peak.res <- peak.res |>
            dplyr::group_by(sample_name) |>
            dplyr::mutate(y = dplyr::cur_group_id())
        peak.res$y_max <- max(peak.res$y) + 1

        lvs <-  c(lvs,unique(peak.res$sample))
    }

    if(!is.null(links_data)){
        if(is.null(links_sample_order)){
            lvs <- c(lvs,unique(links_data$sample))
        }else{
            lvs <- c(lvs,links_sample_order)
        }
    }

    peak.res$sample <- factor(peak.res$sample,
                              levels = c(lvs,"Gene structure"))

    if(!is.null(links_data)){
        links_data$sample <- factor(links_data$sample,
                                    levels = c(lvs,"Gene structure"))
    }


    cov.strc$sample <- factor(cov.strc$sample, levels = c(lvs,"Gene structure"))
    rg$sample <- factor(rg$sample, levels = c(lvs,"Gene structure"))
    tid$sample <- factor(tid$sample, levels = c(lvs,"Gene structure"))
    strc$sample <- factor(strc$sample, levels = c(lvs,"Gene structure"))



    if(!is.null(group_order)){
        cov.strc$group <- factor(cov.strc$group, levels = c(group_order,"Gene structure"))
        rg$group <- factor(rg$group, levels = group_order)
        tid$group <- factor(tid$group, levels = c(group_order,"Gene structure"))
        strc$group <- factor(strc$group, levels = c(group_order,"Gene structure"))
    }

    if(!is.null(gene_order)){
        cov.strc$gene_name <- factor(cov.strc$gene_name, levels = gene_order)
        rg$gene_name <- factor(rg$gene_name, levels = gene_order)
        tid$gene_name <- factor(tid$gene_name, levels = gene_order)
        strc$gene_name <- factor(strc$gene_name, levels = gene_order)
    }

    if(!is.null(region_order)){
        cov.strc$target_region <- factor(cov.strc$target_region, levels = region_order)
        rg$target_region <- factor(rg$target_region, levels = region_order)
        tid$target_region <- factor(tid$target_region, levels = region_order)
        strc$target_region <- factor(strc$target_region, levels = region_order)
    }

    # ==========================================================================
    # plot
    if(length(unique(cov.res$group)) > 1){
        p <-
            ggplot(cov.strc) +
            geom_area(aes(x = pos, y = rpm, fill = group),
                      stat = "identity",alpha = 0.5,show.legend = FALSE)
    }else{
        p <-
            ggplot(cov.strc) +
            geom_area(aes(x = pos, y = rpm, fill = sample),
                      stat = "identity",alpha = 1,show.legend = FALSE)
    }

    p <- p + sample_col

    # check peaks data
    p <- p +
        ggnewscale::new_scale_color() +
        geom_segment(data = peak.res,
                     aes(x = start,xend = end,y = y,yend = y,color = sample_name),
                     linewidth = peaks_width) +
        geom_blank(data = peak.res,aes(y = y_max))


    p <- p + peaks_sample_col


    # check links data
    if(!is.null(links_data)){
        p <- p +
            ggnewscale::new_scale_color() +
            ggbio::geom_arch(data = links_data,
                             aes(x = start,xend = end,
                                 height = -scales::rescale(end - start,to = c(0,1)),
                                 color = .data[[links_col_var]]),linewidth = 0.5)

        p <- p + links_col
    }

    # check highlight region
    if(!is.null(highlight_region)){
        p <- p +
            ggnewscale::new_scale_fill() +
            geom_rect(data = highlight_region,
                      aes(xmin = start,xmax = end,ymin = -Inf,ymax = Inf,
                          fill = .data[[highlight_col_aes]]),alpha = highlight_alpha)

        if(!is.null(highlight_col)){
            p <- p +
                scale_fill_manual(values = highlight_col)
        }
    }



    if(collapse_exon == TRUE){
        tid2 <- tid |>
            dplyr::group_by(.data[[var]],sample,group) |>
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
                                        linewidth = arrow_linewidth)

        if(show_utr == TRUE){
            p <- p +
                # CDS layer
                geom_segment(data = subset(strc, type == "CDS"),
                             aes(x = start, xend = end,y = 1, yend = 1),
                             color = exon_col,
                             linewidth = exon_linewidth) +
                # UTR layer
                geom_segment(data = subset(strc, type != "CDS"),
                             aes(x = start, xend = end,y = 1, yend = 1),
                             color = exon_col,
                             linewidth = exon_linewidth/2)
        }else{
            p <- p +
                # exon layer
                geom_segment(data = strc,
                             aes(x = start, xend = end,y = 1, yend = 1),
                             color = exon_col,
                             linewidth = exon_linewidth)
        }

    }else{

        p <- p +
            # arrow layer
            ggarrow::geom_arrow_segment(data = tid,
                                        aes(x = x, xend = xend,y = transcript_rank, yend = transcript_rank),
                                        color = arrow_col,
                                        arrow_mid = ggarrow::arrow_head_line(),
                                        mid_place = seq(0,1,0.1),
                                        length = arrow_length,
                                        linewidth = arrow_linewidth)

        if(show_utr == TRUE){
            p <- p +
                # CDS layer
                geom_segment(data = subset(strc, type == "CDS"),
                             aes(x = start, xend = end,y = transcript_rank, yend = transcript_rank),
                             color = exon_col,
                             linewidth = exon_linewidth) +
                # UTR layer
                geom_segment(data = subset(strc, type != "CDS"),
                             aes(x = start, xend = end,y = transcript_rank, yend = transcript_rank),
                             color = exon_col,
                             linewidth = exon_linewidth/2)
        }else{
            p <- p +
                # exon layer
                geom_segment(data = strc,
                             aes(x = start, xend = end,y = transcript_rank, yend = transcript_rank),
                             color = exon_col,
                             linewidth = exon_linewidth)
        }

    }


    # add gene label
    if(add_gene_label == TRUE){
        if(collapse_exon == TRUE){
            tid$transcript_rank <- 1

            tid <- tid |>
                dplyr::group_by(gene_name) |>
                dplyr::slice_max(order_by = width)
        }

        if(add_backsqure == TRUE){
            p <- p +
                geom_label(data = tid,
                           aes(x = (start + end)/2,y = transcript_rank,
                               label = .data[[gene_label_aes]]),size = gene_label_size)
        }else{
            p <- p +
                geom_text(data = tid,
                          aes(x = (start + end)/2,y = transcript_rank,
                              label = .data[[gene_label_aes]]),size = gene_label_size)
        }

    }


    # extend region
    if(!is.null(target_gene)){
        p <- p +
            geom_blank(data = gid,aes(x = x_pos, y = 0))
    }else{
        p <- p +
            geom_blank(data = data.frame(target_region = target_region),
                       aes(y = 0))
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
        xlab("Genomic coordinate (kb)") + ylab(ylb)


    # check range label
    if(add_range_label == TRUE){
        p <- p +
            # range layer
            geom_text(data = rg,
                      aes(x = I(range_pos[1]),y = I(range_pos[2]),label = label),
                      show.legend = FALSE)
    }

    # remove y label
    if(remove_labelY == TRUE){
        p <- p +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())
    }

    return(p)
}
