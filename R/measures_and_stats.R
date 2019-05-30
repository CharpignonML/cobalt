#Balance measures

balance.table <- function(C, weights, treat, continuous, binary, s.d.denom, un = FALSE,  
                          s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, 
                          types = NULL, pooled.sds = NULL, quick = TRUE, treat.type = "binary",
                          to.compute, to.disp, to.threshold, ...) {
    #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
    if (no.adj) weight.names <- NULL
    else weight.names <- names(weights)
    names(s.d.denom) <- weight.names
    
    args <- list(...)
    
    if (treat.type == "binary") {
        t.levels <- c("0", "1")
        t.levels_ <- list("0", "1")
        default.stat <- "mean.diffs"
        available.measures <- available.measures.bintreat()
    }
    else {
        t.levels_ <- NULL
        t.levels_ <- list(t.levels)
        default.stat <- "corrs"
        available.measures <- available.measures.conttreat()
    }
    
    measures <- sapply(to.compute[to.compute %in% available.measures], get0, simplify = FALSE)
    stats <- sapply(to.compute[to.compute %in% available.stats()], get0, simplify = FALSE)
    
    measures.to.threshold <- measures[to.threshold]
    
    n.measures <- length(measures)
    n.stats <- length(stats)
    
    to.disp.as.stat <- sapply(to.disp, function(x) tolower(get0(x)$names["stat"]))
    default.stat.as.stat <- tolower(get0(default.stat)$names["stat"])
    
    #B=Balance frame
    Bnames <- c("Type", 
                expand.grid_string(c(expand.grid_string(unlist(lapply(stats, function(S) S$names["stat"])), unlist(t.levels), collapse = "."), 
                                     unlist(lapply(measures, function(M) c(M$names["stat"], 
                                                                           paste.(toupper(M$names["abbrev"]), "Threshold"))))
                ),
                c("Un", weight.names), collapse = "."))
    B <- as.data.frame(matrix(nrow = NCOL(C), ncol = length(Bnames)))
    colnames(B) <- Bnames
    rownames(B) <- colnames(C)
    
    #Set var type (binary/continuous)
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- get.types(C)
    
    
    thresholds <- setNames(rep(NA_real_, length(measures.to.threshold)), 
                                      sapply(measures.to.threshold, function(M) paste.(M$names["abbrev"], "threshold")))
    disp <- setNames(rep(FALSE, n.measures + n.stats), 
                                c(sapply(measures, function(M) tolower(M$names["stat"])),
                                  sapply(stats, function(S) tolower(S$names["stat"]))))
    disp[names(disp) %in% c(default.stat.as.stat, to.disp.as.stat)] <- TRUE
    disp.cols <- setNames(rep(TRUE, ncol(B)), names(B))
    
    #Stats
    for (S in stats) {
        stat <- S$names["stat"]
        
        for (t in t.levels_) {
            B[[paste.(stat, t, "Un")]] <- S$fun(treat = treat, 
                                                C = C, 
                                                weights = NULL, 
                                                types = B[["Type"]], 
                                                continuous = continuous, 
                                                binary = binary, 
                                                s.weights = s.weights, 
                                                pooled.sds = pooled.sds,
                                                treat.level = t,
                                                ...)
            if (!un && !no.adj) {
                disp.cols[paste.(stat, t, "Un")] <- FALSE
            }
            
        }
        
        if (!no.adj) {
            for (j in weight.names) {
                for (t in t.levels_) {
                    B[[paste.(stat, t, j)]] <- S$fun(treat = treat, 
                                                     C = C,
                                                     weights = weights[[j]], 
                                                     types = B[["Type"]], 
                                                     continuous = continuous, 
                                                     binary = binary, 
                                                     s.weights = s.weights, 
                                                     pooled.sds = pooled.sds,
                                                     treat.level = t,
                                                     ...)
                }
            }
        }
        if (!any(sapply(B[expand.grid_string(stat, t.levels, c("Un", weight.names), collapse = ".")], is.finite))) {
            disp[tolower(stat)] <- FALSE
        }
        
        if (!disp[tolower(stat)]) {
            disp.cols[names(disp.cols) %in% expand.grid_string(
                stat, t.levels,
                c("Un", weight.names), collapse = "."
            )] <- FALSE
        }
    }
    
    #Measures
    for (M in measures) {
        stat <- M$names["stat"]
        abbrev <- M$names["abbrev"]
        ABBREV <- toupper(M$names["abbrev"])
        threshold <- args[[paste.(abbrev, "threshold")]]
        
        if (abs) a0 <- M$abs
        else a0 <- base::identity
        
        B[[paste.(stat, "Un")]] <- a0(M$fun(treat = treat, 
                                            C = C, 
                                            weights = NULL, 
                                            types = B[["Type"]], 
                                            continuous = continuous, 
                                            binary = binary, 
                                            s.d.denom=s.d.denom[1], 
                                            s.weights = s.weights, 
                                            pooled.sds = pooled.sds,
                                            ...))
        if (!un && !no.adj) {
            disp.cols[paste.(stat, "Un")] <- FALSE
        }
        
        if (!no.adj) {
            for (j in weight.names) {
                B[[paste.(stat, j)]] <- a0(M$fun(treat = treat, 
                                                 C = C,
                                                 weights = weights[[j]], 
                                                 types = B[["Type"]], 
                                                 continuous = continuous, 
                                                 binary = binary, 
                                                 s.d.denom = s.d.denom[j], 
                                                 s.weights = s.weights, 
                                                 pooled.sds = pooled.sds,
                                                 ...))
            }
        }
        
        if (!any(sapply(B[paste.(stat, c("Un", weight.names))], is.finite))) {
            disp[tolower(stat)] <- FALSE
            threshold <- NULL
        }
        
        if (is_not_null(threshold)) {
            
            assert.threshold.out.of.limits(threshold, M)
            threshold <- M$abs(threshold)
            
            if (no.adj) {
                B[[paste.(ABBREV, "Threshold", "Un")]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.(stat, "Un")]]), paste0(ifelse(M$abs(B[[paste.(stat, "Un")]]) < threshold, "Balanced, <", "Not Balanced, >"), round(threshold, 3)), "")
            }
            else {
                for (i in weight.names) {
                    B[[paste.(ABBREV, "Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.(stat, i)]]), paste0(ifelse(m$abs(B[[paste.(stat, i)]]) < threshold, "Balanced, <", "Not Balanced, >"), round(threshold, 3)), "")
                }
            }
            thresholds[paste.(abbrev, "threshold")] <- threshold
        }
        else {
            B[names(B) %in% paste.(ABBREV, "Threshold", c("Un", weight.names))] <- NULL
            disp.cols <- disp.cols[names(disp.cols) %nin% paste.(ABBREV, "Threshold", c("Un", weight.names))]
        }
        
        if (!disp[tolower(stat)]) {
            disp.cols[names(disp.cols) %in% expand.grid_string(
                c(stat, paste.(ABBREV, "Threshold")),
                c("Un", weight.names), collapse = "."
            )] <- FALSE
        }
        
        if (no.adj || NCOL(weights) <= 1) {
            names(B)[names(B) == paste.(ABBREV, "Threshold", "Adj")] <- paste.(ABBREV, "Threshold")
            names(disp.cols)[names(disp.cols) == paste.(ABBREV, "Threshold", "Adj")] <- paste.(ABBREV, "Threshold")
        }
        
    }
    
    attr(B, "disp") <- disp
    attr(B, "thresholds") <- thresholds
    attr(B, "disp.cols") <- disp.cols
    
    return(B)
    
}

#Sample Stats
available.stats <- function() c("means", "sds", "pop.diffs")

means <- setNames(vector("list", 2),
                  c("fun", "names"))
means$fun <- function(treat, C, weights, s.weights, treat.level = NULL, ...) {
    if (is_null(weights)) weights <- rep(1, length(treat))
    weights <- weights*s.weights
    
    if (is_null(treat.level)) in.treat <- rep(TRUE, length(treat))
    else in.treat <- treat == treat.level
    
    m <- col.w.m(C[in.treat, , drop = FALSE], w = weights[in.treat])
    return(m)
}
means$names <- c(stat = "Mean")

sds <- setNames(vector("list", 2),
                c("fun", "names"))
sds$fun <- function(treat, C, weights, s.weights, binary, types, treat.level = NULL, ...) {
    if (is_null(weights)) weights <- rep(1, length(treat))
    weights <- weights*s.weights
    
    if (missing(binary) || is_null(binary)) {
        binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
    }
    else binary <- match_arg(binary, c("raw", "std"))
    
    if (is_null(treat.level)) in.treat <- rep(TRUE, length(treat))
    else in.treat <- treat == treat.level
    
    sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else types != "Binary"
    s <- rep(NA_real_, NCOL(C))
    s[sd.computable] <- sqrt(col.w.v(C[in.treat, sd.computable, drop = FALSE], w = weights[in.treat]))
    
    return(s)
}
sds$names <- c(stat = "SD")

pop.diffs <- setNames(vector("list", 2),
                      c("fun", "names"))
pop.diffs$fun <- function(treat, C, types, weights, continuous, binary, pooled.sds = NULL, s.weights, treat.level = NULL, ...) {
    if (is_null(weights)) weights <- rep(1, length(treat))
    w <- weights*s.weights
    sw <- s.weights
    
    if (is_null(treat.level)) in.treat <- rep(TRUE, length(treat))
    else in.treat <- treat == treat.level
    
    #Check continuous and binary
    if (missing(continuous) || is_null(continuous)) {
        continuous <- match_arg(getOption("cobalt_continuous", "std"), c("std", "raw"))
    }
    else continuous <- match_arg(continuous, c("std", "raw"))
    if (missing(binary) || is_null(binary)) {
        binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
    }
    else binary <- match_arg(binary, c("raw", "std"))
    
    
    diffs <- col.w.m(C[in.treat, , drop = FALSE], w[in.treat]) - 
        col.w.m(C, w)
    diffs[check_if_zero(diffs)] <- 0
    denoms <- rep(1, NCOL(C))
    denoms.to.std <- ifelse(diffs == 0, FALSE, ifelse(types == "Binary", binary == "std", continuous == "std"))
    
    if (any(denoms.to.std)) {
        
        if (is_not_null(pooled.sds)) {
            denoms[denoms.to.std] <- pooled.sds[denoms.to.std]
        }
        else {
            denoms[denoms.to.std] <-  sqrt(.5*(col.w.v(C[in.treat, denoms.to.std, drop = FALSE], sw[in.treat]) +
                                                   col.w.v(C[!in.treat, denoms.to.std, drop = FALSE], sw[!in.treat])))
        }
        
    }
    
    std.diffs <- diffs/denoms
    if (any(!is.finite(std.diffs))) {
        warning("Some standardized mean differences were not finite. This can result from no variation in one of the treatment groups.", call. = FALSE)
        std.diffs[!is.finite(std.diffs)] <- NA_real_
    }
    
    return(std.diffs)
}
pop.diffs$names <- c(stat = "Pop.Diff")

#Balance Measures
available.measures.bintreat <- function() c("mean.diffs", "ks.statistics", "v.ratios", "log.sd.diffs")
available.measures.conttreat <- function() c("corr")

mean.diffs <- setNames(vector("list", 3),
                       c("fun", "abs", "names"))
mean.diffs$fun <- function(treat, C, types, weights, continuous, binary, s.d.denom, pooled.sds = NULL, s.weights, ...) {
    if (is_null(weights)) weights <- rep(1, length(treat))
    w <- weights*s.weights
    sw <- s.weights
    
    #Check continuous and binary
    if (missing(continuous) || is_null(continuous)) {
        continuous <- match_arg(getOption("cobalt_continuous", "std"), c("std", "raw"))
    }
    else continuous <- match_arg(continuous, c("std", "raw"))
    if (missing(binary) || is_null(binary)) {
        binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
    }
    else binary <- match_arg(binary, c("raw", "std"))
    
    
    diffs <- col.w.m(C[treat == 1, , drop = FALSE], w[treat == 1]) - 
        col.w.m(C[treat == 0, , drop = FALSE], w[treat == 0])
    diffs[check_if_zero(diffs)] <- 0
    denoms <- rep(1, NCOL(C))
    denoms.to.std <- ifelse(diffs == 0, FALSE, ifelse(types == "Binary", binary == "std", continuous == "std"))
    
    if (any(denoms.to.std)) {
        if (s.d.denom == "control") {
            denoms[denoms.to.std] <- sqrt(col.w.v(C[treat == 0, denoms.to.std, drop = FALSE], sw[treat == 0]))
        }
        else if (s.d.denom == "treated") {
            denoms[denoms.to.std] <- sqrt(col.w.v(C[treat == 1, denoms.to.std, drop = FALSE], sw[treat == 1]))
        }
        else if (s.d.denom == "pooled") {
            if (is_not_null(pooled.sds)) {
                denoms[denoms.to.std] <- pooled.sds[denoms.to.std]
            }
            else {
                denoms[denoms.to.std] <-  sqrt(.5*(col.w.v(C[treat == 0, denoms.to.std, drop = FALSE], sw[treat == 0]) +
                                                       col.w.v(C[treat == 1, denoms.to.std, drop = FALSE], sw[treat == 1])))
            }
        }
    }
    
    std.diffs <- diffs/denoms
    if (any(!is.finite(std.diffs))) {
        warning("Some standardized mean differences were not finite. This can result from no variation in one of the treatment groups.", call. = FALSE)
        std.diffs[!is.finite(std.diffs)] <- NA_real_
    }
    
    return(std.diffs)
}
mean.diffs$abs <- function(x) abs(x)
mean.diffs$names <- c(abbrev = "m", #m.threshold
                      what = "Means", #Balanced Means, Max Imbalance Means
                      stat = "Diff", #Diff.Adj
                      variable.with.greatest = "mean difference"
)
mean.diffs$threshold.limits <- c(min = 0, max = Inf)

ks.statistics <- setNames(vector("list", 3),
                          c("fun", "abs", "names"))
ks.statistics$fun <- function(treat, C, types, weights, s.weights, ...) {
    if (is_null(weights)) weights <- rep(1, length(treat))
    weights <- weights*s.weights
    
    weights[treat == 1] <- weights[treat==1]/sum(weights[treat==1])
    weights[treat == 0] <- -weights[treat==0]/sum(weights[treat==0])
    
    non.binary <- types != "Binary"
    ks[non.binary] <- apply(C[, non.binary, drop = FALSE], 2, function(x_) {
        x <- x_[!is.na(x_)]
        ordered.index <- order(x)
        cumv <- abs(cumsum(weights[ordered.index]))[diff(x[ordered.index]) != 0]
        return(if (is_null(cumv)) 0 else max(cumv))
    })
    return(ks)
}
ks.statistics$abs <- function(x) identity(x)
ks.statistics$names <- c(abbrev = "ks", #ks.threshold
                         what = "KS Statistics", #Balanced Means, Max Imbalance Means
                         stat = "KS", #KS.Adj
                         variable.with.greatest = "KS statistic"
) 
ks.statistics$threshold.limits <- c(min = 0, max = 1)

v.ratios <- setNames(vector("list", 3),
                     c("fun", "abs", "names"))
v.ratios$fun <- function(treat, C, types, weights, s.weights, ...) {
    if (is_null(weights)) weights <- rep(1, length(treat))
    weights <- weights*s.weights
    
    ratios <- rep(NA_real_, NCOL(C))
    non.binary <- types != "Binary"
    ratios[non.binary] <- col.w.v(C[treat == 1, non.binary, drop = FALSE], w = weights[treat == 1]) / col.w.v(C[treat == 0, non.binary, drop = FALSE], w = weights[treat == 0])
    return(ratios)
}
v.ratios$abs <- function(x) pmax(x, 1/x)
v.ratios$names <- c(abbrev = "v", #ks.threshold
                    what = "Variances", #Balanced Means, Max Imbalance Means
                    stat = "V.Ratio", #KS.Adj
                    variable.with.greatest = "variance ratio"
) 
v.ratios$threshold.limits <- c(min = 1, max = Inf)


log.sd.diffs <- setNames(vector("list", 3),
                         c("fun", "abs", "names"))
log.sd.diffs$fun <- function(treat, C, types, weights, s.weights, ...) {
    if (is_null(weights)) weights <- rep(1, length(treat))
    weights <- weights*s.weights
    
    log.sd.diffs <- rep(NA_real_, NCOL(C))
    non.binary <- types != "Binary"
    log.sd.diffs[non.binary] <- log(sqrt(col.w.v(C[treat == 1, non.binary, drop = FALSE], w = weights[treat == 1]) / col.w.v(C[treat == 0, non.binary, drop = FALSE], w = weights[treat == 0])))
    return(log.sd.diffs)
}
log.sd.diffs$abs <- function(x) abs(x)
log.sd.diffs$names <- c(abbrev = "sd", #ks.threshold
                        what = "Log Standard Deviations", #Balanced Means, Max Imbalance Means
                        stat = "L.SD.Diff", #KS.Adj
                        variable.with.greatest = "log standard deviation difference"
) 
log.sd.diffs$threshold.limits <- c(min = 0, max = Inf)

corrs <- setNames(vector("list", 3),
                       c("fun", "abs", "names"))
corrs$fun <- function(treat, C, types, weights, s.weights, ...) {
    if (is_null(weights)) weights <- rep(1, length(treat))
    w <- weights*s.weights

    corr <- apply(C, 2, w.r, y = treat, w = weights, s.weights = s.weights)
    
    return(corr)
}
corrs$abs <- function(x) abs(x)
corrs$names <- c(abbrev = "r", #m.threshold
                      what = "Correlations", #Balanced Means, Max Imbalance Means
                      stat = "Corr", #Diff.Adj
                 variable.with.greatest = "treatment correlation"
) 
corrs$threshold.limits <- c(min = 0, max = 1)


print.bal.tab <- function(x, disp.m.threshold = "as.is", disp.v.threshold = "as.is", disp.ks.threshold = "as.is", disp.r.threshold = "as.is", imbalanced.only = "as.is", un = "as.is", disp.bal.tab = "as.is", disp.means = "as.is", disp.sds = "as.is", disp.v.ratio = "as.is", disp.ks = "as.is", digits = max(3, getOption("digits") - 3), ...) {
    
    call <- x$call
    balance <- x$Balance
    nn <- x$Observations
    p.ops <- x$print.options
    
    #Prevent exponential notation printing
    op <- options(scipen=getOption("scipen"))
    options(scipen = 999)
    on.exit(options(op))
    
    #Adjustments to print options
    if (!identical(un, "as.is") && p.ops$disp.adj) {
        if (!is.logical(un)) stop("un must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$un == FALSE && un == TRUE) {
            warning("un cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$un <- un
    }
    if (!identical(disp.means, "as.is")) {
        if (!is.logical(disp.means)) stop("disp.means must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.means == FALSE && disp.means == TRUE) {
            warning("disp.means cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.means <- disp.means
    }
    if (!identical(disp.sds, "as.is")) {
        if (!is.logical(disp.sds)) stop("disp.sds must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.sds == FALSE && disp.sds == TRUE) {
            warning("disp.sds cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.sds <- disp.sds
    }
    if (!identical(disp.v.ratio, "as.is")) {
        if (!is.logical(disp.v.ratio)) stop("disp.v.ratio must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.v.ratio == FALSE && disp.v.ratio == TRUE) {
            warning("disp.v.ratio cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.v.ratio <- disp.v.ratio
    }
    if (!identical(disp.ks, "as.is")) {
        if (!is.logical(disp.ks)) stop("disp.ks must be TRUE, FALSE, or \"as.is\"")
        if (p.ops$quick && p.ops$disp.ks == FALSE && disp.ks == TRUE) {
            warning("disp.ks cannot be set to TRUE if quick = TRUE in the original object.", call. = FALSE)
        }
        else p.ops$disp.ks <- disp.ks
    }
    if (!identical(disp.r.threshold, "as.is")) {
        if (!is.logical(disp.r.threshold)) stop("disp.r.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$r.threshold) && !disp.r.threshold) {
            p.ops$r.threshold <- NULL
            baltal.r <- NULL
            maximbal.r <- NULL
        }
    }
    if (!identical(disp.m.threshold, "as.is")) {
        if (!is.logical(disp.m.threshold)) stop("disp.m.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$m.threshold) && !disp.m.threshold) {
            p.ops$m.threshold <- NULL
            baltal.m <- NULL
            maximbal.m <- NULL
        }
    }
    if (!identical(disp.v.threshold, "as.is")) {
        if (!is.logical(disp.v.threshold)) stop("disp.v.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$v.threshold) && !disp.v.threshold) {
            p.ops$v.threshold <- NULL
            baltal.v <- NULL
            maximbal.v <- NULL
        }
    }
    if (is_null(p.ops$disp.v.ratio)|| !p.ops$disp.v.ratio) {
        p.ops$v.threshold <- NULL
        baltal.v <- NULL
        maximbal.v <- NULL
    }
    if (!identical(disp.ks.threshold, "as.is")) {
        if (!is.logical(disp.ks.threshold)) stop("disp.ks.threshold must be FALSE or \"as.is\"")
        if (is_not_null(p.ops$ks.threshold) && !disp.ks.threshold) {
            p.ops$ks.threshold <- NULL
            baltal.ks <- NULL
            maximbal.ks <- NULL
        }
    }
    if (is_null(p.ops$disp.ks) || !p.ops$disp.ks) {
        p.ops$ks.threshold <- NULL
        baltal.ks <- NULL
        maximbal.ks <- NULL
    }
    if (!identical(disp.bal.tab, "as.is")) {
        if (!is.logical(disp.bal.tab)) stop("disp.bal.tab must be TRUE, FALSE, or \"as.is\"")
        p.ops$disp.bal.tab <- disp.bal.tab
    }
    if (p.ops$disp.bal.tab) {
        if (!identical(imbalanced.only, "as.is")) {
            if (!is.logical(imbalanced.only)) stop("imbalanced.only must be TRUE, FALSE, or \"as.is\"")
            p.ops$imbalanced.only <- imbalanced.only
        }
        if (p.ops$imbalanced.only) {
            if (all(sapply(c(p.ops$m.threshold, 
                             p.ops$v.threshold, 
                             p.ops$ks.threshold, 
                             p.ops$r.threshold), is_null))) {
                warning("A threshold must be specified if imbalanced.only = TRUE. Displaying all covariates.", call. = FALSE)
                p.ops$imbalanced.only <- FALSE
            }
        }
    }
    else p.ops$imbalanced.only <- FALSE
    
    keep <- names(attr(balance, "disp.cols"))[attr(balance, "disp.cols")]
    
    S <- setNames(lapply(names(p.ops$thresholds), function(th) {
        measure <- get0(available.measures.bintreat()[sapply(available.measures.bintreat(), function(m) {
            paste.(get0(m)$names["abbrev"], "threshold") == th
        })])
        list(threshold = p.ops$thresholds[th],
             Names = measure$names["what"],
             Baltal = paste0(measure$names["variable.with.greatest"], "s"),
             Maximbal = measure$names["variable.with.greatest"])
        
    }), names(p.ops$thresholds))
    
    baltal. <- setNames(lapply(S, function(s) x[[paste.("Balanced", s[["Names"]])]]), names(p.ops$thresholds))
    maximbal. <- setNames(lapply(S, function(s) x[[paste.("Max.Imbalance", s[["Names"]])]]), names(p.ops$thresholds))

    if (is_not_null(call)) {
        cat(underline("Call") %+% "\n " %+% paste(deparse(call), collapse = "\n") %+% "\n\n")
    }
    
    if (p.ops$disp.bal.tab) {
        if (p.ops$imbalanced.only) {
            keep.row <- rowSums(apply(balance[grepl(".Threshold", names(balance), fixed = TRUE)], 2, function(x) !is.na(x) & startsWith(x, "Not Balanced"))) > 0
        }
        else keep.row <- rep(TRUE, nrow(balance))
        
        cat(underline("Balance Measures") %+% "\n")
        if (all(!keep.row)) cat(italic("All covariates are balanced.") %+% "\n")
        else print.data.frame_(round_df_char(balance[keep.row, keep], digits))
        cat("\n")
    }
    
    for (s_ in names(p.ops$thresholds)) {
        if (is_not_null(baltal.[[s_]])) {
            cat(underline("Balance tally for " %+% S[[s_]][["Baltal"]]) %+% "\n")
            print.data.frame_(baltal.[[s_]])
            cat("\n")
        }
        if (is_not_null(maximbal.[[s_]])) {
            cat(underline("Variable with the greatest " %+% S[[s_]][["Maximbal"]]) %+% "\n")
            print.data.frame_(round_df_char(maximbal.[[s_]], digits), row.names = FALSE)
            cat("\n")
        }
    }

    if (is_not_null(nn)) {
        for (i in rownames(x$Observations)) {
            if (all(x$Observations[i,] == 0)) x$Observations <- x$Observations[rownames(x$Observations)!=i, , drop = FALSE]
        }
        if ("Matched (Unweighted)" %in% rownames(x$Observations) && all(check_if_zero(x$Observations["Matched",] - x$Observations["Matched (Unweighted)",]))) x$Observations <- x$Observations[rownames(x$Observations)!="Matched (Unweighted)", , drop = FALSE]
        cat(underline(attr(x$Observations, "tag")) %+% "\n")
        print.warning <- FALSE
        if (length(attr(x$Observations, "ss.type")) > 1 && nunique.gt(attr(x$Observations, "ss.type")[-1], 1)) {
            ess <- ifelse(attr(x$Observations, "ss.type") == "ess", "*", "")
            x$Observations <- setNames(cbind(x$Observations, ess), c(names(x$Observations), ""))
            print.warning <- TRUE
        }
        print.data.frame_(round_df_char(x$Observations, digits = max(0, digits-1)))
        if (print.warning) cat(italic("* indicates effective sample size"))
    }
    invisible(x)
}

assert.threshold.out.of.limits <- function(given.threshold, measure) {
    limits <- measure$threshold.limits
    if (!between(given.threshold, limits, inclusive = FALSE)) {
        
        threshname <- paste.(measure$names["abbrev"], "threshold")
        
        if (limits["min"] == -Inf && limits["max"] == Inf) {
            stop(paste0(threshname, " must be a real number."), call. = FALSE)
        }
        else if (limits["min"] == -Inf && limits["max"] != Inf) {
            stop(paste0(threshname, " must be less than ", limits["max"], "."), call. = FALSE)
        }
        else if (limits["min"] != -Inf && limits["max"] == Inf) {
            stop(paste0(threshname, " must be greater than ", limits["min"], "."), call. = FALSE)
            
        }
        else if (limits["min"] != -Inf && limits["max"] != Inf) {
            stop(paste0(threshname, " must be between ", limits["min"], " and ", limits["max"], "."), call. = FALSE)
        }
    }
}