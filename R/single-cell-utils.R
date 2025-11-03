globalVariables(c("Assay", "AverageCounts", "CellsPerGroup", "DefaultAssay", "Idents", "Idents<-",
                  "StringToGRanges", "UpdateChromatinObject", "WhichCells", "coverage", "granges",
                  "hue_pal", "mcols", "mcols<-", "norm.value", "position", "roll_sum", "slice_sample",
                  "subsetByOverlaps", "theme_browser"))

# SetIfNull <- getFromNamespace("SetIfNull","Signac")
# FindRegion <- getFromNamespace("FindRegion","Signac")
# GetGroups <- getFromNamespace("GetGroups","Signac")
# CutMatrix <- getFromNamespace("CutMatrix","Signac")
# ApplyMatrixByGroup <- getFromNamespace("ApplyMatrixByGroup","Signac")

SetIfNull <- function(x, y) {
    if (is.null(x = x)) {
        return(y)
    } else {
        return(x)
    }
}

FindRegion <- function(
        object,
        region,
        sep = c("-", "-"),
        assay = NULL,
        extend.upstream = 0,
        extend.downstream = 0) {
    if (!methods::is(object = region, class2 = "GRanges")) {
        # first try to convert to coordinates, if not lookup gene
        region <- tryCatch(
            expr = suppressWarnings(
                expr = StringToGRanges(regions = region, sep = sep)
            ),
            error = function(x) {
                region <- Signac::LookupGeneCoords(
                    object = object,
                    assay = assay,
                    gene = region
                )
                return(region)
            }
        )
        if (is.null(x = region)) {
            stop("Gene not found")
        }
    }
    region <- suppressWarnings(expr = Signac::Extend(
        x = region,
        upstream = extend.upstream,
        downstream = extend.downstream
    )
    )
    return(region)
}

GetGroups <- function(
        object,
        group.by,
        idents) {
    if (is.null(x = group.by)) {
        obj.groups <- Idents(object = object)
    } else {
        obj.md <- object[[group.by]]
        obj.groups <- obj.md[, 1]
        names(obj.groups) <- rownames(x = obj.md)
    }
    if (!is.null(idents)) {
        obj.groups <- obj.groups[obj.groups %in% idents]
    }
    return(obj.groups)
}

isRemote <- function(x) {
    return(grepl(pattern = "^http|^ftp", x = x))
}

GetIndexFile <- function(fragment, verbose = TRUE) {
    is.remote <- isRemote(x = fragment)
    index.filepaths <- c(paste0(fragment, ".tbi"),
                         paste0(fragment, ".csi"))
    index.file <- index.filepaths[file.exists(index.filepaths)]
    if (length(x = index.file) == 0 & !is.remote) {
        stop("Fragment file is not indexed.")
    } else if(length(x = index.file) == 0) {
        if (verbose) {
            message("Fragment file is on a remote server")
        }
        index.file = paste0(fragment, ".tbi")
    } else if (length(x = index.file) == 2) {
        if (verbose) {
            message("TBI and CSI index both present, using TBI index")
        }
        index.file <- index.file[1]
    }
    return(index.file)
}

TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
    if (record.ident) {
        nrep <- S4Vectors::elementNROWS(x = reads)
    }
    original_names = names(reads)
    reads <- unlist(x = reads, use.names = FALSE)
    if (length(x = reads) == 0 | is.null(x = original_names)) {
        df <- data.frame(
            "chr" = "",
            "start" = "",
            "end" = "",
            "cell" = "",
            "count" = ""
        )
        df <- df[-1, ]
        return(df)
    }
    reads <- stringr::stri_split_fixed(str = reads, pattern = "\t")
    n <- length(x = reads[[1]])
    unlisted <- unlist(x = reads)
    e1 <- unlisted[n * (seq_along(along.with = reads)) - (n - 1)]
    e2 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 2)])
    e3 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 3)])
    e4 <- unlisted[n * (seq_along(along.with = reads)) - (n - 4)]
    e5 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 5)])
    df <- data.frame(
        "chr" = e1,
        "start" = e2,
        "end" = e3,
        "cell" = e4,
        "count" = e5,
        stringsAsFactors = FALSE,
        check.rows = FALSE,
        check.names = FALSE
    )
    if (record.ident) {
        df$ident <- rep(x = seq_along(along.with = nrep), nrep)
    }
    return(df)
}

#' @importFrom Rsamtools scanTabix
GetReadsInRegion <- function(
        cellmap,
        region,
        tabix.file,
        cells = NULL,
        verbose = TRUE,
        ...) {
    file.to.object <- names(x = cellmap)
    names(x = file.to.object) <- cellmap

    if (verbose) {
        message("Extracting reads in requested region")
    }
    if (!methods::is(object = region, class2 = "GRanges")) {
        region <- StringToGRanges(regions = region, ...)
    }
    # remove regions that aren't in the fragment file
    common.seqlevels <- intersect(
        x = GenomeInfoDb::seqlevels(x = region),
        y = seqnamesTabix(file = tabix.file)
    )
    if (length(common.seqlevels) != 0) {
        region <- keepSeqlevels(
            x = region,
            value = common.seqlevels,
            pruning.mode = "coarse"
        )
        reads <- scanTabix(file = tabix.file, param = region)
        reads <- TabixOutputToDataFrame(reads = reads)
        reads <- reads[
            fastmatch::fmatch(x = reads$cell, table = cellmap, nomatch = 0L) > 0,
        ]
        # convert cell names to match names in object
        reads$cell <- file.to.object[reads$cell]
        if (!is.null(x = cells)) {
            reads <- reads[reads$cell %in% cells, ]
        }
        if (nrow(reads) == 0) {
            reads$ident <- integer()
            reads$length <- numeric()
            return(reads)
        }
        reads$length <- reads$end - reads$start
    } else {
        reads <- data.frame(
            "chr" = character(),
            "start" = numeric(),
            "end" = numeric(),
            "cell" = character(),
            "count" = numeric(),
            "ident" = integer(),
            "length" = numeric()
        )
    }
    return(reads)
}

#' @importFrom Matrix sparseMatrix
SingleFileCutMatrix <- function(
        cellmap,
        region,
        cells = NULL,
        tabix.file,
        verbose = TRUE) {
    # if multiple regions supplied, must be the same width
    cells <- SetIfNull(x = cells, y = names(x = cellmap))
    if (length(x = region) == 0) {
        return(NULL)
    }
    fragments <- GetReadsInRegion(
        region = region,
        cellmap = cellmap,
        cells = cells,
        tabix.file = tabix.file,
        verbose = verbose
    )
    start.lookup <- start(x = region)
    names(start.lookup) <- seq_along(region)
    # if there are no reads in the region
    # create an empty matrix of the correct dimension
    if (nrow(x = fragments) == 0) {
        cut.matrix <- sparseMatrix(
            i = NULL,
            j = NULL,
            dims = c(length(x = cells), width(x = region)[[1]])
        )
    } else {
        fragstarts <- start.lookup[fragments$ident] + 1
        cut.df <- data.frame(
            position = c(fragments$start, fragments$end) - fragstarts,
            cell = c(fragments$cell, fragments$cell),
            stringsAsFactors = FALSE
        )
        cut.df <- cut.df[
            (cut.df$position > 0) & (cut.df$position <= width(x = region)[[1]]),
        ]
        cell.vector <- seq_along(along.with = cells)
        names(x = cell.vector) <- cells
        cell.matrix.info <- cell.vector[cut.df$cell]
        cut.matrix <- sparseMatrix(
            i = cell.matrix.info,
            j = cut.df$position,
            x = 1,
            dims = c(length(x = cells), width(x = region)[[1]])
        )
    }
    rownames(x = cut.matrix) <- cells
    colnames(x = cut.matrix) <- seq_len(width(x = region)[[1]])
    return(cut.matrix)
}

#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom Signac Fragments GetFragmentData
CutMatrix <- function(
        object,
        region,
        group.by = NULL,
        assay = NULL,
        cells = NULL,
        verbose = TRUE) {
    # run SingleFileCutMatrix for each fragment file and combine results
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    cells <- SetIfNull(x = cells, y = colnames(x = object))
    fragments <- Fragments(object = object[[assay]])
    if (length(x = fragments) == 0) {
        stop("No fragment information found for requested assay")
    }
    res <- list()
    for (i in seq_along(along.with = fragments)) {
        fragment.path <- GetFragmentData(object = fragments[[i]], slot = "path")
        cellmap <- GetFragmentData(object = fragments[[i]], slot = "cells")
        tabix.file <- TabixFile(
            file = fragment.path,
            index = GetIndexFile(fragment = fragment.path, verbose = FALSE)
        )
        open(con = tabix.file)
        # remove regions that aren't in the fragment file
        seqnames.in.both <- intersect(
            x = seqnames(x = region),
            y = seqnamesTabix(file = tabix.file)
        )
        region <- keepSeqlevels(
            x = region,
            value = seqnames.in.both,
            pruning.mode = "coarse"
        )
        if (length(x = region) != 0) {
            cm <- SingleFileCutMatrix(
                region = region,
                cellmap = cellmap,
                tabix.file = tabix.file,
                cells = cells,
                verbose = FALSE
            )
            res[[i]] <- cm
        }
        close(con = tabix.file)
    }
    res <- Reduce(f = `+`, x = res)
    return(res)
}


ApplyMatrixByGroup <- function(
        mat,
        groups,
        fun,
        normalize = TRUE,
        group.scale.factors = NULL,
        scale.factor = NULL) {
    if (normalize) {
        if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
            stop("If normalizing counts, supply group scale factors")
        }
    }
    all.groups <- as.character(x = unique(x = groups))
    if (any(is.na(x = groups))) {
        all.groups <- c(all.groups, NA)
    }
    ngroup <- length(x = all.groups)
    npos <- ncol(x = mat)

    group <- unlist(
        x = lapply(X = all.groups, FUN = function(x) rep(x, npos))
    )
    position <- rep(x = as.numeric(x = colnames(x = mat)), ngroup)
    count <- vector(mode = "numeric", length = npos * ngroup)

    for (i in seq_along(along.with = all.groups)) {
        grp <- all.groups[[i]]
        if (is.na(x = grp)) {
            pos.cells <- names(x = groups)[is.na(x = groups)]
        } else {
            pos.cells <- names(x = groups)[groups == all.groups[[i]]]
        }
        if (length(x = pos.cells) > 1) {
            totals <- fun(x = mat[pos.cells, ])
        } else {
            totals <- mat[pos.cells, ]
        }
        count[((i - 1) * npos + 1):((i * npos))] <- totals
    }

    # construct dataframe
    coverages <- data.frame(
        "group" = group, "position" = position, "count" = count,
        stringsAsFactors = FALSE
    )

    if (normalize) {
        scale.factor <- SetIfNull(
            x = scale.factor, y = median(x = group.scale.factors)
        )
        coverages$norm.value <- coverages$count /
            group.scale.factors[coverages$group] * scale.factor
    } else {
        coverages$norm.value <- coverages$count
    }
    return(coverages)
}

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
        region <- StringToGRanges(regions = region)
    }
    if (is.null(x = peaks)) {
        peaks <- granges(x = object[[assay]])
        md <- object[[assay]][[]]
        mcols(x = peaks) <- md
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
    peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
    peak.df <- as.data.frame(x = peak.intersect)

    start.pos <- start(x = region)
    end.pos <- end(x = region)
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
    window.size <- width(x = region)
    levels.use <- levels(x = obj.groups)
    chromosome <- as.character(x = GenomicRanges::seqnames(x = region))
    start.pos <- start(x = region)
    end.pos <- end(x = region)
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
            coverages <- group_by(.data = coverages, group)
            coverages <- mutate(.data = coverages, coverage = roll_sum(
                x = norm.value, n = window, fill = NA, align = "center"
            ))
            coverages <- ungroup(x = coverages)
        } else {
            coverages$coverage <- coverages$norm.value
        }

        coverages <- coverages[!is.na(x = coverages$coverage), ]
        coverages <- group_by(.data = coverages, group)
        sampling <- min(max.downsample, window.size * downsample.rate)
        set.seed(seed = 1234)
        coverages <- slice_sample(.data = coverages, n = as.integer(x = sampling))
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
        colors_all <- hue_pal()(length(x = levels.use))
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
        gr <- GRanges(
            seqnames = chromosome,
            IRanges(start = start.pos, end = end.pos)
        )
        if (multicov) {
            p <- ggplot(
                data = coverages,
                mapping = aes(x = position, y = coverage, fill = Assay)
            )
        } else {
            p <- ggplot(
                data = coverages,
                mapping = aes(x = position, y = coverage, fill = group)
            )
        }
        p <- p +
            geom_area(
                stat = "identity",
                alpha = ifelse(test = !split.assays & multicov, yes = 0.5, no = 1)) +
            geom_hline(yintercept = 0, size = 0.1)
        if (split.assays) {
            p <- p +
                facet_wrap(facets = ~assay_group, strip.position = "left", ncol = 1)
        } else {
            p <- p + facet_wrap(facets = ~group, strip.position = "left", ncol = 1)
        }
        p <- p +
            xlab(label = paste0(chromosome, " position (bp)")) +
            ylab(label = paste0("Normalized signal \n(range ",
                                as.character(x = ymin), " - ",
                                as.character(x = ymax), ")")) +
            ylim(c(ymin, ymax)) +
            theme_browser(legend = multicov) +
            theme(panel.spacing.y = unit(x = 0, units = "line"))
        if (!is.null(x = levels.use) & !multicov) {
            p <- p + scale_fill_manual(values = colors_all)
        }
        if (!is.null(x = region.highlight)) {
            if (!inherits(x = region.highlight, what = "GRanges")) {
                warning("region.highlight must be a GRanges object")
            } else {
                md <- mcols(x = region.highlight)
                if ("color" %in% colnames(x = md)) {
                    color.use <- md$color
                } else {
                    color.use <- rep(x = "grey", length(x = region.highlight))
                }
                df <- data.frame(
                    "start" = start(x = region.highlight),
                    "end" = end(x = region.highlight),
                    "color" = color.use
                )
                df$start <- ifelse(
                    test = df$start < start.pos,
                    yes = start.pos,
                    no = df$start
                )
                df$end <- ifelse(
                    test = df$end > end.pos,
                    yes = end.pos,
                    no = df$end
                )
                p <- p +
                    geom_rect(
                        data = df,
                        inherit.aes = FALSE,
                        aes_string(
                            xmin = "start",
                            xmax = "end",
                            ymin = 0,
                            ymax = ymax),
                        fill = rep(x = df$color, length(x = unique(x = coverages$group))),
                        color = "transparent",
                        alpha = 0.2
                    )
            }
        }
        return(p)
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
        Idents(object = object) <- group.by
    }
    if (!is.null(x = idents)) {
        ident.cells <- WhichCells(object = object, idents = idents)
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
        grouping.var <- Idents(object = object)
        combined.var <- paste0(object[[split.by]][, 1], "_", grouping.var)
        object$grouping_tmp <- combined.var
        Idents(object = object) <- "grouping_tmp"
        group.by <- "grouping_tmp"
        if (!is.null(x = idents)) {
            # adjust idents parameter with new split.by information
            idents.keep <- combined.var[grouping.var %in% idents]
            idents <- unique(x = idents.keep)
        }
    }
    cells.per.group <- CellsPerGroup(
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
        reads.per.group <- AverageCounts(
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
        colnames(cutmat) <- start(x = region):end(x = region)
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
