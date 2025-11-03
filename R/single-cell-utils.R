globalVariables(c("Assay", "AverageCounts", "CellsPerGroup", "DefaultAssay", "Idents", "Idents<-",
                  "StringToGRanges", "UpdateChromatinObject", "WhichCells", "coverage", "granges",
                  "hue_pal", "mcols", "mcols<-", "norm.value", "position", "roll_sum", "slice_sample",
                  "subsetByOverlaps", "theme_browser"))

SetIfNull <- getFromNamespace("SetIfNull","Signac")
FindRegion <- getFromNamespace("FindRegion","Signac")
GetGroups <- getFromNamespace("GetGroups","Signac")
CutMatrix <- getFromNamespace("CutMatrix","Signac")
ApplyMatrixByGroup <- getFromNamespace("ApplyMatrixByGroup","Signac")
isRemote <- getFromNamespace("isRemote","Signac")
GetIndexFile <- getFromNamespace("GetIndexFile","Signac")
TabixOutputToDataFrame <- getFromNamespace("TabixOutputToDataFrame","Signac")
GetReadsInRegion <- getFromNamespace("GetReadsInRegion","Signac")
SingleFileCutMatrix <- getFromNamespace("SingleFileCutMatrix","Signac")
UpdateChromatinObject <- getFromNamespace("UpdateChromatinObject","Signac")



# slightly modified from Signac source code
PeakPlot2 <- function(
        object,
        region,
        assay = NULL,
        peaks = NULL,
        group.by = NULL,
        color = "dimgrey",
        sep = c("-", "-"),
        extend.upstream = 0,
        extend.downstream = 0) {
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
        stop("The requested assay is not a ChromatinAssay.")
    }

    if (!inherits(x = region, what = "GRanges")) {
        region <- Signac::StringToGRanges(regions = region)
    }
    if (is.null(x = peaks)) {
        peaks <- GenomicRanges::granges(x = object[[assay]])
        md <- object[[assay]][[]]
        S4Vectors::mcols(x = peaks) <- md
    }
    region <- FindRegion(
        object = object,
        region = region,
        sep = sep,
        assay = assay[[1]],
        extend.upstream = extend.upstream,
        extend.downstream = extend.downstream
    )
    # subset to covered range
    peak.intersect <- IRanges::subsetByOverlaps(x = peaks, ranges = region)
    peak.df <- as.data.frame(x = peak.intersect)

    start.pos <- GenomicRanges::start(x = region)
    end.pos <- GenomicRanges::end(x = region)
    chromosome <- GenomicRanges::seqnames(x = region)

    if (nrow(x = peak.df) > 0) {
        if (!is.null(x = group.by)) {
            if (!(group.by %in% colnames(x = peak.df))) {
                warning("Requested grouping variable not found")
                group.by <- NULL
            }
        }
        peak.df$start[peak.df$start < start.pos] <- start.pos
        peak.df$end[peak.df$end > end.pos] <- end.pos

        return(peak.df)

        # peak.plot <- ggplot(data = peak.df,
        #   aes_string(color = SetIfNull(x = group.by, y = "color"))) +
        #   geom_segment(aes(x = start, y = 0, xend = end, yend = 0),
        #                size = 2,
        #                data = peak.df)
    } else {
        # no peaks present in region, make empty panel
        # peak.plot <- ggplot(data = peak.df)
        return(data.frame())
    }

    # peak.plot <- peak.plot + theme_classic() +
    #   ylab(label = "Peaks") +
    #   theme(axis.ticks.y = element_blank(),
    #         axis.text.y = element_blank()) +
    #   xlab(label = paste0(chromosome, " position (bp)")) +
    #   xlim(c(start.pos, end.pos))
    # if (is.null(x = group.by)) {
    #   # remove legend, change color
    #   peak.plot <- peak.plot +
    #     scale_color_manual(values = color) +
    #     theme(legend.position = "none")
    # }
    # return(peak.plot)
}




# slightly modified from Signac source code
CoverageTrack2 <- function(
        cutmat,
        region,
        group.scale.factors,
        scale.factor,
        assay.scale,
        obj.groups,
        ymax,
        downsample.rate,
        split.assays = FALSE,
        region.highlight = NULL,
        window = 100,
        max.downsample = 3000,
        return.data = TRUE) {
    window.size <- IRanges::width(x = region)
    levels.use <- levels(x = obj.groups)
    chromosome <- as.character(x = GenomeInfoDb::seqnames(x = region))
    start.pos <- GenomicRanges::start(x = region)
    end.pos <- GenomicRanges::end(x = region)
    multicov <- length(x = cutmat) > 1

    cov.df <- data.frame()
    for (i in seq_along(along.with = cutmat)) {
        coverages <- ApplyMatrixByGroup(
            mat = cutmat[[i]],
            fun = colSums,
            groups = obj.groups,
            group.scale.factors = group.scale.factors[[i]],
            scale.factor = scale.factor[[i]],
            normalize = TRUE
        )
        if (!is.na(x = window)) {
            coverages <- dplyr::group_by(.data = coverages, group)
            coverages <- dplyr::mutate(.data = coverages, coverage = RcppRoll::roll_sum(
                x = norm.value, n = window, fill = NA, align = "center"
            ))
            coverages <- dplyr::ungroup(x = coverages)
        } else {
            coverages$coverage <- coverages$norm.value
        }

        coverages <- coverages[!is.na(x = coverages$coverage), ]
        coverages <- dplyr::group_by(.data = coverages, group)
        sampling <- min(max.downsample, window.size * downsample.rate)
        set.seed(seed = 1234)
        coverages <- dplyr::slice_sample(.data = coverages, n = as.integer(x = sampling))
        coverages$Assay <- names(x = cutmat)[[i]]
        if (multicov) {
            if (assay.scale == "separate") {
                # scale to fraction of max for each separately
                assay.max <- max(coverages$coverage, na.rm = TRUE)
                coverages$coverage <- coverages$coverage / assay.max
            }
        }
        cov.df <- rbind(cov.df, coverages)
    }
    coverages <- cov.df
    coverages$Assay <- factor(x = coverages$Assay, levels = names(x = cutmat))
    coverages$assay_group <- paste(coverages$group, coverages$Assay, sep = "_")

    # restore factor levels
    if (!is.null(x = levels.use)) {
        colors_all <- scales::hue_pal()(length(x = levels.use))
        names(x = colors_all) <- levels.use
        coverages$group <- factor(x = coverages$group, levels = levels.use)
    }
    covmax <- signif(x = max(coverages$coverage, na.rm = TRUE), digits = 2)
    if (is.null(x = ymax)) {
        ymax <- covmax
    } else if (is.character(x = ymax)) {
        if (!startsWith(x = ymax, prefix = "q")) {
            stop("Unknown ymax requested. Must be NULL, a numeric value, or
           a quantile denoted by 'qXX' with XX the desired quantile value,
           e.g. q95 for 95th percentile")
        }
        percentile.use <- as.numeric(
            x = sub(pattern = "q", replacement = "", x = as.character(x = ymax))
        ) / 100
        ymax <- covmax * percentile.use
    }
    ymin <- 0

    # perform clipping
    coverages$coverage[coverages$coverage > ymax] <- ymax

    # returns
    if(return.data == TRUE){
        coverages <- coverages[,c("group","position","coverage", "Assay")]
        colnames(coverages) <- c("sample","pos","rpm","group")

        return(coverages)
    }else{
        # gr <- GenomicRanges::GRanges(
        #     seqnames = chromosome,
        #     IRanges::IRanges(start = start.pos, end = end.pos)
        # )
        # if (multicov) {
        #     p <- ggplot(
        #         data = coverages,
        #         mapping = aes(x = position, y = coverage, fill = Assay)
        #     )
        # } else {
        #     p <- ggplot(
        #         data = coverages,
        #         mapping = aes(x = position, y = coverage, fill = group)
        #     )
        # }
        # p <- p +
        #     geom_area(
        #         stat = "identity",
        #         alpha = ifelse(test = !split.assays & multicov, yes = 0.5, no = 1)) +
        #     geom_hline(yintercept = 0, size = 0.1)
        # if (split.assays) {
        #     p <- p +
        #         facet_wrap(facets = ~assay_group, strip.position = "left", ncol = 1)
        # } else {
        #     p <- p + facet_wrap(facets = ~group, strip.position = "left", ncol = 1)
        # }
        # p <- p +
        #     xlab(label = paste0(chromosome, " position (bp)")) +
        #     ylab(label = paste0("Normalized signal \n(range ",
        #                         as.character(x = ymin), " - ",
        #                         as.character(x = ymax), ")")) +
        #     ylim(c(ymin, ymax)) +
        #     theme_browser(legend = multicov) +
        #     theme(panel.spacing.y = unit(x = 0, units = "line"))
        # if (!is.null(x = levels.use) & !multicov) {
        #     p <- p + scale_fill_manual(values = colors_all)
        # }
        # if (!is.null(x = region.highlight)) {
        #     if (!inherits(x = region.highlight, what = "GRanges")) {
        #         warning("region.highlight must be a GRanges object")
        #     } else {
        #         md <- mcols(x = region.highlight)
        #         if ("color" %in% colnames(x = md)) {
        #             color.use <- md$color
        #         } else {
        #             color.use <- rep(x = "grey", length(x = region.highlight))
        #         }
        #         df <- data.frame(
        #             "start" = start(x = region.highlight),
        #             "end" = end(x = region.highlight),
        #             "color" = color.use
        #         )
        #         df$start <- ifelse(
        #             test = df$start < start.pos,
        #             yes = start.pos,
        #             no = df$start
        #         )
        #         df$end <- ifelse(
        #             test = df$end > end.pos,
        #             yes = end.pos,
        #             no = df$end
        #         )
        #         p <- p +
        #             geom_rect(
        #                 data = df,
        #                 inherit.aes = FALSE,
        #                 aes_string(
        #                     xmin = "start",
        #                     xmax = "end",
        #                     ymin = 0,
        #                     ymax = ymax),
        #                 fill = rep(x = df$color, length(x = unique(x = coverages$group))),
        #                 color = "transparent",
        #                 alpha = 0.2
        #             )
        #     }
        # }
        # return(p)
    }
}




# slightly modified from Signac source code
SingleCoveragePlot2 <- function(
        object,
        region,
        features = NULL,
        assay = NULL,
        split.assays = FALSE,
        assay.scale = "common",
        show.bulk = FALSE,
        expression.assay = NULL,
        expression.slot = "data",
        annotation = TRUE,
        peaks = TRUE,
        peaks.group.by = NULL,
        ranges = NULL,
        ranges.group.by = NULL,
        ranges.title = "Ranges",
        region.highlight = NULL,
        links = TRUE,
        tile = FALSE,
        tile.size = 100,
        tile.cells = 100,
        bigwig = NULL,
        bigwig.type = "coverage",
        bigwig.scale = "common",
        group.by = NULL,
        split.by = NULL,
        window = 100,
        extend.upstream = 0,
        extend.downstream = 0,
        ymax = NULL,
        scale.factor = NULL,
        cells = NULL,
        idents = NULL,
        sep = c(":", "-"),
        heights = NULL,
        max.downsample = 3000,
        downsample.rate = 0.1,
        return.data = TRUE) {
    valid.assay.scale <- c("common", "separate")
    if (!(assay.scale %in% valid.assay.scale)) {
        stop(
            "Unknown assay.scale requested. Please choose from: ",
            paste(valid.assay.scale, collapse = ", ")
        )
    }
    cells <- SetIfNull(x = cells, y = colnames(x = object))
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    if (!inherits(x = assay, what = "list")) {
        assay <- list(assay)
    }
    lapply(X = assay, FUN = function(x) {
        if (!inherits(x = object[[x]], what = "ChromatinAssay")) {
            stop("Requested assay is not a ChromatinAssay.")
        }
    })
    if(length(colnames(object)) > length(colnames(object[[assay[[1]]]]))) {
        object <- UpdateChromatinObject(object = object,
                                        chromatin.assay = assay,
                                        expression.assay = expression.assay,
                                        features = features)
    }
    if (!is.null(x = group.by)) {
        Seurat::Idents(object = object) <- group.by
    }
    if (!is.null(x = idents)) {
        ident.cells <- SeuratObject::WhichCells(object = object, idents = idents)
        cells <- intersect(x = cells, y = ident.cells)
    }
    region <- FindRegion(
        object = object,
        region = region,
        sep = sep,
        assay = assay[[1]],
        extend.upstream = extend.upstream,
        extend.downstream = extend.downstream
    )
    if (!is.null(x = split.by)) {
        # combine split.by and group.by information
        grouping.var <- Seurat::Idents(object = object)
        combined.var <- paste0(object[[split.by]][, 1], "_", grouping.var)
        object$grouping_tmp <- combined.var
        Seurat::Idents(object = object) <- "grouping_tmp"
        group.by <- "grouping_tmp"
        if (!is.null(x = idents)) {
            # adjust idents parameter with new split.by information
            idents.keep <- combined.var[grouping.var %in% idents]
            idents <- unique(x = idents.keep)
        }
    }
    cells.per.group <- Signac::CellsPerGroup(
        object = object,
        group.by = group.by
    )

    obj.groups <- GetGroups(
        object = object,
        group.by = group.by,
        idents = idents
    )

    # subset to used cells
    obj.groups <- obj.groups[cells]

    cm.list <- list()
    sf.list <- list()
    gsf.list <- list()
    for (i in seq_along(along.with = assay)) {
        reads.per.group <- Signac::AverageCounts(
            object = object,
            assay = assay[[i]],
            group.by = group.by,
            verbose = FALSE
        )
        cutmat <- CutMatrix(
            object = object,
            region = region,
            assay = assay[[i]],
            cells = cells,
            verbose = FALSE
        )
        colnames(cutmat) <- IRanges::start(x = region):IRanges::end(x = region)
        group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
        scale.factor <- SetIfNull(
            x = scale.factor, y = median(x = group.scale.factors)
        )
        cm.list[[i]] <- cutmat
        sf.list[[i]] <- scale.factor
        gsf.list[[i]] <- group.scale.factors
    }
    names(x = cm.list) <- unlist(x = assay)

    # check return data
    if(return.data == TRUE){
        cov <- CoverageTrack2(
            cutmat = cm.list,
            region = region,
            group.scale.factors = gsf.list,
            scale.factor = sf.list,
            window = window,
            ymax = ymax,
            split.assays = split.assays,
            assay.scale = assay.scale,
            obj.groups = obj.groups,
            region.highlight = region.highlight,
            downsample.rate = downsample.rate,
            max.downsample = max.downsample,
            return.data = TRUE
        )

        peak.df <- PeakPlot2(object = object,
                             assay = assay[[1]],
                             region = region,
                             group.by = peaks.group.by)

        return(list(cov,peak.df))
    }else{
        # p <-
        #     CoverageTrack2(
        #         cutmat = cm.list,
        #         region = region,
        #         group.scale.factors = gsf.list,
        #         scale.factor = sf.list,
        #         window = window,
        #         ymax = ymax,
        #         split.assays = split.assays,
        #         assay.scale = assay.scale,
        #         obj.groups = obj.groups,
        #         region.highlight = region.highlight,
        #         downsample.rate = downsample.rate,
        #         max.downsample = max.downsample,
        #         return.data = FALSE
        #     )
        #
        # # create bigwig tracks
        # if (!is.null(x = bigwig)) {
        #     if (!inherits(x = bigwig, what = "list")) {
        #         warning("BigWig should be a list of file paths")
        #         bigwig <- list("bigWig" = bigwig)
        #     }
        #     if (length(x = bigwig.type) == 1) {
        #         bigwig.type <- rep(x = bigwig.type, length(x = bigwig))
        #     } else if (length(x = bigwig.type) != length(x = bigwig)) {
        #         stop("Must supply a bigWig track type for each bigWig file")
        #     }
        #     unique.types <- unique(x = bigwig.type)
        #     bw.all <- list()
        #     for (i in seq_along(unique.types)) {
        #         bw.use <- which(x = bigwig.type == unique.types[[i]])
        #         bw.all[[i]] <- BigwigTrack(
        #             region = region,
        #             bigwig = bigwig[bw.use],
        #             type = unique.types[[i]],
        #             bigwig.scale = bigwig.scale,
        #             ymax = ymax
        #         )
        #     }
        #     bigwig.tracks <- CombineTracks(
        #         plotlist = bw.all,
        #         heights = table(unlist(x = bigwig.type))
        #     )
        # } else {
        #     bigwig.tracks <- NULL
        # }
        # if (!is.null(x = features)) {
        #     ex.plot <- ExpressionPlot(
        #         object = object,
        #         features = features,
        #         assay = expression.assay,
        #         idents = idents,
        #         group.by = group.by,
        #         slot = expression.slot
        #     )
        #     widths <- c(10, length(x = features))
        # } else {
        #     ex.plot <- NULL
        #     widths <- NULL
        # }
        # if (is.logical(x = annotation)) {
        #     if (annotation) {
        #         gene.plot <- Signac::AnnotationPlot(
        #             object = object[[assay[[1]]]],
        #             region = region,
        #             mode = "gene"
        #         )
        #     } else {
        #         gene.plot <- NULL
        #     }
        # } else {
        #     gene.plot <- Signac::AnnotationPlot(
        #         object = object[[assay[[1]]]],
        #         region = region,
        #         mode = annotation
        #     )
        # }
        # if (is.character(x = links)) {
        #     # subset to genes in the desired list
        #     links.use <- Links(object = object)
        #     links.use <- links.use[links.use$gene %in% links]
        #     Links(object = object) <- links.use
        #     links <- TRUE
        # }
        # if (links) {
        #     link.plot <- LinkPlot(object = object[[assay[[1]]]], region = region)
        # } else {
        #     link.plot <- NULL
        # }
        # if (peaks) {
        #     peak.plot <- PeakPlot(
        #         object = object,
        #         assay = assay[[1]],
        #         region = region,
        #         group.by = peaks.group.by
        #     )
        # } else {
        #     peak.plot <- NULL
        # }
        # if (!is.null(x = ranges)) {
        #     range.plot <- PeakPlot(
        #         object = object,
        #         assay = assay[[1]],
        #         region = region,
        #         peaks = ranges,
        #         group.by = ranges.group.by,
        #         color = "brown3") +
        #         ylab(ranges.title)
        # } else {
        #     range.plot <- NULL
        # }
        # if (tile) {
        #     # reuse cut matrix
        #     # TODO implement for multi assay
        #     tile.df <- ComputeTile(
        #         cutmatrix = cm.list[[1]],
        #         groups = obj.groups,
        #         window = tile.size,
        #         n = tile.cells,
        #         order = "total"
        #     )
        #     tile.plot <- CreateTilePlot(
        #         df = tile.df,
        #         n = tile.cells
        #     )
        # } else {
        #     tile.plot <- NULL
        # }
        # if (show.bulk) {
        #     object$bulk <- "All cells"
        #     reads.per.group <- AverageCounts(
        #         object = object,
        #         assay = assay[[1]],
        #         group.by = "bulk",
        #         verbose = FALSE
        #     )
        #     cells.per.group <- CellsPerGroup(
        #         object = object,
        #         group.by = "bulk"
        #     )
        #     bulk.scale.factor <- suppressWarnings(reads.per.group * cells.per.group)
        #     bulk.groups <- rep(x = "All cells", length(x = obj.groups))
        #     names(x = bulk.groups) <- names(x = obj.groups)
        #     bulk.plot <- CoverageTrack(
        #         cutmat = cm.list,
        #         region = region,
        #         group.scale.factors = list(bulk.scale.factor),
        #         scale.factor = scale.factor,
        #         window = window,
        #         ymax = ymax,
        #         obj.groups = bulk.groups,
        #         downsample.rate = downsample.rate,
        #         max.downsample = max.downsample
        #     ) +
        #         scale_fill_manual(values = "grey") +
        #         ylab("")
        # } else {
        #     bulk.plot <- NULL
        # }
        # nident <- length(x = unique(x = obj.groups))
        # if (split.assays) {
        #     nident <- nident * length(x = assay)
        # }
        # bulk.height <- (1 / nident) * 10
        # bw.height <- 10
        # heights <- SetIfNull(
        #     x = heights, y = c(10, bulk.height, bw.height, 10, 3, 1, 1, 3)
        # )
        # p <- CombineTracks(
        #     plotlist = list(p, bulk.plot, bigwig.tracks, tile.plot, gene.plot,
        #                     peak.plot, range.plot, link.plot),
        #     expression.plot = ex.plot,
        #     heights = heights,
        #     widths = widths
        # ) & theme(
        #     legend.key.size = unit(x = 1/2, units = "lines"),
        #     legend.text = element_text(size = 7),
        #     legend.title = element_text(size = 8)
        # )
        # return(p)
    }
}
