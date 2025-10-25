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






####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
    return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
    weights <- pos_res <- list()
    for(i in 1:nrow(dat)){
        m <- mu[-i,!is.na(dat[i,])]
        x <- dat[i,!is.na(dat[i,])]
        gamma_sub <- gamma[-i]
        phi_sub <- phi[-i]

        # take a subset of genes to do integration - save time
        if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
            if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
            mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
            m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
            G_sub <- gene.subset.n
        }else{
            if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
            G_sub <- nrow(dat)-1
        }

        #LH <- sapply(1:G_sub, function(j){sum(log2(dnbinom(x, mu=m[j,], size=1/phi_sub[j])+1))})
        LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
        LH[is.nan(LH)]=0;
        if(sum(LH)==0 | is.na(sum(LH))){
            pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i]))
        }else{
            pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
        }

        weights[[i]] <- as.matrix(LH/sum(LH))
    }
    pos_res <- do.call(rbind, pos_res)
    weights <- do.call(cbind, weights)
    res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"], weights=weights)
    return(res)
}


####  Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
    new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
    for(a in 1:nrow(counts_sub)){
        for(b in 1:ncol(counts_sub)){
            if(counts_sub[a, b] <= 1){
                new_counts_sub[a,b] <- counts_sub[a, b]
            }else{
                tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
                if(abs(tmp_p-1)<1e-4){
                    new_counts_sub[a,b] <- counts_sub[a, b]
                    # for outlier count, if p==1, will return Inf values -> use original count instead
                }else{
                    new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
                }
            }
        }
    }
    return(new_counts_sub)
}



mapDisp <- function(old_mu, new_mu, old_phi, divider){
    new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
    for(a in 1:nrow(old_mu)){
        for(b in 1:ncol(old_mu)){
            old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
            new_var <- old_var / (divider[a, b]^2)
            new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
        }
    }
    return(new_phi)
}


#' Adjust for batch effects using an empirical Bayes framework in RNA-seq raw counts
#'
#' ComBat_seq is an improved model from ComBat using negative binomial regression,
#' which specifically targets RNA-Seq count data.
#'
#' @param counts Raw count matrix from genomic studies (dimensions gene x sample)
#' @param batch Vector / factor for batch
#' @param group Vector / factor for biological condition of interest
#' @param covar_mod Model matrix for multiple covariates to include in linear model (signals from these variables are kept in data after adjustment)
#' @param full_mod Boolean, if TRUE include condition of interest in model
#' @param shrink Boolean, whether to apply shrinkage on parameter estimation
#' @param shrink.disp Boolean, whether to apply shrinkage on dispersion
#' @param gene.subset.n Number of genes to use in empirical Bayes estimation, only useful when shrink = TRUE
#'
#' @return data A gene x sample count matrix, adjusted for batch effects.
#'
#' @importFrom edgeR DGEList estimateGLMCommonDisp estimateGLMTagwiseDisp glmFit glmFit.default getOffset
#' @importFrom stats dnbinom lm pnbinom qnbinom
#' @importFrom utils capture.output
#'
#' @examples
#' \dontrun{
#' count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0,1), 4)
#'
#' # include condition (group variable)
#' adjusted_counts <- ComBat_seq_devel(count_matrix, batch=batch, group=group, full_mod=TRUE)
#'
#' # do not include condition
#' adjusted_counts <- ComBat_seq_devel(count_matrix, batch=batch, group=NULL, full_mod=FALSE)
#' }
#'
#' @export
ComBat_seq_devel <- function(counts, batch, group=NULL, covar_mod=NULL, full_mod=TRUE,
                             shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL){
    ########  Preparation  ########
    ## Does not support 1 sample per batch yet
    batch <- as.factor(batch)
    if(any(table(batch)<=1)){
        stop("ComBat-seq doesn't support 1 sample per batch yet")
    }

    ## Remove genes with only 0 counts in any batch
    keep_lst <- lapply(levels(batch), function(b){
        which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
    })
    keep <- Reduce(intersect, keep_lst)
    rm <- setdiff(1:nrow(counts), keep)
    countsOri <- counts
    counts <- counts[keep, ]

    # require bioconductor 3.7, edgeR 3.22.1
    dge_obj <- edgeR::DGEList(counts=counts)

    ## Prepare characteristics on batches
    n_batch <- nlevels(batch)  # number of batches
    batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch
    n_batches <- sapply(batches_ind, length)
    #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
    n_sample <- sum(n_batches)
    cat("Found",n_batch,'batches\n')

    ## Make design matrix
    # batch
    batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
    # covariate
    group <- as.factor(group)
    if(full_mod & nlevels(group)>1){
        cat("Using full model in ComBat-seq.\n")
        mod <- model.matrix(~group)
    }else{
        cat("Using null model in ComBat-seq.\n")
        mod <- model.matrix(~1, data=as.data.frame(t(counts)))
    }
    # drop intercept in covariate model
    if(!is.null(covar_mod)){
        if(is.data.frame(covar_mod)){
            covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i){model.matrix(~covar_mod[,i])}))
        }
        covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]
    }
    # bind with biological condition of interest
    mod <- cbind(mod, covar_mod)
    # combine
    design <- cbind(batchmod, mod)

    ## Check for intercept in covariates, and drop if present
    check <- apply(design, 2, function(x) all(x == 1))
    #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
    design <- as.matrix(design[,!check])
    cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')

    ## Check if the design is confounded
    if(qr(design)$rank<ncol(design)){
        #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
        if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")}
        if(ncol(design)>(n_batch+1)){
            if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
            }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")}}
    }

    ## Check for missing values in count matrix
    NAs = any(is.na(counts))
    if(NAs){cat(c('Found',sum(is.na(counts)),'Missing Data Values\n'),sep=' ')}


    ########  Estimate gene-wise dispersions within each batch  ########
    cat("Estimating dispersions\n")
    ## Estimate common dispersion within each batch as an initial value
    disp_common <- sapply(1:n_batch, function(i){
        if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){
            # not enough residual degree of freedom
            return(edgeR::estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
        }else{
            return(edgeR::estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
        }
    })

    ## Estimate gene-wise dispersion within each batch
    genewise_disp_lst <- lapply(1:n_batch, function(j){
        if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
            # not enough residual degrees of freedom - use the common dispersion
            return(rep(disp_common[j], nrow(counts)))
        }else{
            return(edgeR::estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ],
                                                 dispersion=disp_common[j], prior.df=0))
        }
    })
    names(genewise_disp_lst) <- paste0('batch', levels(batch))

    ## construct dispersion matrix
    phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
    for(k in 1:n_batch){
        phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k])
    }


    ########  Estimate parameters from NB GLM  ########
    cat("Fitting the GLM model\n")
    glm_f <- edgeR::glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
    alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
    new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) +   # original offset - sample (library) size
        vec2mat(alpha_g, ncol(counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
    glm_f2 <- edgeR::glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4)

    gamma_hat <- glm_f2$coefficients[, 1:n_batch]
    mu_hat <- glm_f2$fitted.values
    phi_hat <- do.call(cbind, genewise_disp_lst)


    ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########
    if(shrink){
        cat("Apply shrinkage - computing posterior estimates for parameters\n")
        mcint_fun <- monte_carlo_int_NB
        monte_carlo_res <- lapply(1:n_batch, function(ii){
            if(ii==1){
                mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]],
                                   gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)
            }else{
                invisible(capture.output(mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]],
                                                            gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)))
            }
            return(mcres)
        })
        names(monte_carlo_res) <- paste0('batch', levels(batch))

        gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
        gamma_star_mat <- do.call(cbind, gamma_star_mat)
        phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
        phi_star_mat <- do.call(cbind, phi_star_mat)

        if(!shrink.disp){
            cat("Apply shrinkage to mean only\n")
            phi_star_mat <- phi_hat
        }
    }else{
        cat("Shrinkage off - using GLM estimates for parameters\n")
        gamma_star_mat <- gamma_hat
        phi_star_mat <- phi_hat
    }


    ########  Obtain adjusted batch-free distribution  ########
    mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
    for(jj in 1:n_batch){
        mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]])-vec2mat(gamma_star_mat[, jj], n_batches[jj]))
    }
    phi_star <- rowMeans(phi_star_mat)


    ########  Adjust the data  ########
    cat("Adjusting the data\n")
    adjust_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
    for(kk in 1:n_batch){
        counts_sub <- counts[, batches_ind[[kk]]]
        old_mu <- mu_hat[, batches_ind[[kk]]]
        old_phi <- phi_hat[, kk]
        new_mu <- mu_star[, batches_ind[[kk]]]
        new_phi <- phi_star
        adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=counts_sub,
                                                              old_mu=old_mu, old_phi=old_phi,
                                                              new_mu=new_mu, new_phi=new_phi)
    }

    #dimnames(adjust_counts) <- dimnames(counts)
    #return(adjust_counts)

    ## Add back genes with only 0 counts in any batch (so that dimensions won't change)
    adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
    dimnames(adjust_counts_whole) <- dimnames(countsOri)
    adjust_counts_whole[keep, ] <- adjust_counts
    adjust_counts_whole[rm, ] <- countsOri[rm, ]
    return(adjust_counts_whole)
}





# codes from https://github.com/kassambara/survminer/blob/e23fa1229c4e7624739d54d7b76885950ca61f74/R/surv_cutpoint.R
surv_cutpoint <- function(data, time = "time", event = "event", variables,
                          minprop = 0.1, progressbar = TRUE){
    if(!inherits(data, "data.frame"))
        stop("data should be an object of class data.frame")
    data <- as.data.frame(data)
    if(!all(c(time, event) %in% colnames(data)))
        stop("Specify correct column names containing time and event values.")
    if(!all(variables %in% colnames(data)))
        stop("Some variables are not found in the data: ",
             paste(setdiff(variables, colnames(data)), collapse =", "))

    not_numeric_vars <- .get_not_numeric_vars(data[, variables, drop = FALSE])
    variables <- setdiff(variables, not_numeric_vars) # keep only numeric variables
    if(length(variables)==0) stop("At least, one numeric variables required.")

    nvar <- length(variables)
    if(nvar <= 5) progressbar <- FALSE
    if(progressbar) pb <- utils::txtProgressBar(min = 0, max = nvar, style = 3)
    surv_data <- data.frame(time = data[, time], event = data[, event])
    res <- list()
    for (i in 1:nvar){
        var_i <- variables[i]
        surv_data$var <- data[, var_i]
        max_stat_i <- maxstat::maxstat.test(survival::Surv(time, event) ~ var, data = surv_data,
                                            smethod = "LogRank", pmethod="none",
                                            minprop = minprop, maxprop = 1-minprop,
                                            alpha = alpha)
        res[[var_i]] <- max_stat_i
        if(progressbar) utils::setTxtProgressBar(pb, i)
    }
    colnames(surv_data) <- c(time, event)
    res$data <- cbind.data.frame(surv_data[, 1:2, drop = FALSE], data[, variables, drop = FALSE])
    res$minprop <- minprop
    if(!is.null(not_numeric_vars)) res$not_numeric <- data[, not_numeric_vars, drop = FALSE]
    res <- structure(res, class = c("list", "surv_cutpoint"))
    res$cutpoint <- summary(res)
    res
}



# codes from https://github.com/kassambara/survminer/blob/e23fa1229c4e7624739d54d7b76885950ca61f74/R/surv_cutpoint.R
# Get not numeric columns in a data.frame
.get_not_numeric_vars <- function(data_frame){
    is_numeric <- sapply(data_frame, is.numeric)
    if(sum(!is_numeric) == 0) res = NULL
    else res <- colnames(data_frame[, !is_numeric, drop = FALSE])
    res
}




get_trans_pos <- function(strc = NULL){
    # Calculate gene-level genomic spans
    gene_spans <- strc |>
        dplyr::group_by(gene_name) |>
        dplyr::summarise(
            gene_start = min(start),
            gene_end = max(end),
            gene_span = max(end) - min(start)) |>
        arrange(gene_start)

    # Calculate transcript-level spans
    transcript_spans <- strc %>%
        dplyr::group_by(gene_name, transcript_id) |>
        dplyr::summarise(trans_start = min(start),
                         trans_end = max(end),
                         trans_span = max(end) - min(start),
                         .groups = "drop")  |>
        dplyr::left_join(gene_spans, by = "gene_name") |>
        dplyr::arrange(dplyr::desc(gene_span), gene_name, dplyr::desc(trans_span))

    # Assign tracks to avoid overlaps
    assign_tracks <- function(spans_df) {
        spans_df$track <- NA

        for (i in 1:nrow(spans_df)) {
            current_start <- spans_df$trans_start[i]
            current_end <- spans_df$trans_end[i]

            # Find available track
            track_num <- max(c(spans_df$track[!is.na(spans_df$track)], 0)) + 1

            # Check if this track is free for current transcript
            for (t in track_num:1) {
                overlaps <- spans_df %>%
                    filter(!is.na(track), track == t) %>%
                    filter(trans_end >= current_start & trans_start <= current_end)

                if (nrow(overlaps) == 0) {
                    track_num <- t
                } else {
                    break
                }
            }

            spans_df$track[i] <- track_num
        }
        return(spans_df)
    }

    # Assign tracks to transcripts
    transcript_spans <- assign_tracks(transcript_spans)

    max_track <- max(transcript_spans$track)
    transcript_spans <- transcript_spans |>
        dplyr::mutate(transcript_rank = max_track - track + 1) |>
        dplyr::select(gene_name, transcript_id, transcript_rank) |>
        dplyr::mutate(y_max = max(transcript_rank))

    return(transcript_spans)
}


