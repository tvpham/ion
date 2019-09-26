# Author: Thang V Pham
#
# This code contains a modified version of the heatmap.2 function
# in the R gplots package (version 3.0.1)
#
# Copyright Thang Pham, 2018-2019

# NOTE: sourcing this will replace any existing 'ion' object

ion <- list()

# Normalization ----

ion$normalize_global <- function(d, total_count = NULL) {

    total <- if (is.null(total_count)) apply(d, 2, sum) else total_count

    m <- mean(total)
    factor <- total / m;

    tmp <-  d / (matrix(1, nrow(d), 1) %*% factor)

    return(tmp)

}

ion$normalize_median <- function(d_log, fixed_median = NULL) {

    m <- apply(d_log, 2, median, na.rm=TRUE)
    f <- if (is.null(fixed_median)) (mean(m)-m) else (fixed_median-m)
    return(d_log + matrix(1, nrow = nrow(d_log), ncol = 1) %*% as.numeric(f))

}


# Utilities ----

ion$row_ok <- function(d) {
    return(rowSums(is.na(d)) < ncol(d))
}

ion$dev_off <- function() {
    while (dev.off()>1) {}
}

ion$impute <- function(d, method = "constant", value = 0, seed = 1203) {

    n <- sum(is.na(d))

    if (n == 0) {
        return(d)
    } else {

        dd <- d

        if (method == "constant") {
            dd[is.na(dd)] <- value
        } else if (method == "normal") {

            set.seed(seed)
            global_min <- min(d, na.rm=TRUE)
            global_sd <- mean(apply(d, 1, sd, na.rm = TRUE), na.rm = TRUE)
            dd[is.na(dd)] <- rnorm(n,
                                   mean = global_min,
                                   sd = global_sd)
        } else {
            stop("Unknown imputation method.")
        }
        return(dd)
    }
}

ion$load <- function(filename,
                     # a tab-deliminated text file with a header line
                     ...
                     ) {
    return(read.delim(filename, quote = "", ...))
}

ion$save <- function(d,
                    # a data table

                    filename,
                    # e.g. "output.txt", "d.txt"
                    ...
                    ) {
    write.table(d, file = filename, sep = "\t", row.names = FALSE, quote = FALSE, qmethod = "double", ...)
}

ion$boxplot_column <- function(d, ...) {
    dl <- list()
    for (i in 1:ncol(d)) {
        dl[i] <- list(d[!is.na(d[, i]), i])
    }
    boxplot(dl,...)
}

# Statistical testing ----

## limma

ion$logFC_to_fc <- function(logFC, dat, group1, group2, log_base = 10, BIG = 1e4) {

    N <- nrow(dat)

    fc <- rep(1, N)

    dd <- dat[, c(group1,  group2)]

    for (r in (1:N)) {

        if (sum(dd[r, ] == -Inf, na.rm = TRUE) > 0) {
            dd[r, dd[r, ] == -Inf] <- NA
        }

        x <- as.numeric(dd[r, 1:length(group1)])
        y <- as.numeric(dd[r, (length(group1)+1):(length(group1)+length(group2))])

        # re-calculate fc

        mx <- mean(x, na.rm = TRUE)

        my <- mean(y, na.rm = TRUE)

        if (is.nan(mx)) {
            if (is.nan(my)) {
                fc[r] <- 1
            } else {
                fc[r] <- BIG
            }
        } else {
            if (is.nan(my)) {
                fc[r] <- -BIG
            } else {
                if (is.na(logFC[r])) {
                    fc[r] <- 1
                } else {
                    real_fc <- log_base^logFC[r]
                    if (real_fc > 1) {
                        fc[r] <- real_fc
                    } else {
                        fc[r] <- -1.0/real_fc
                    }
                }
            }
        }

    }

    return(fc)
}

ion$limma_2g <- function(dat, group1, group2) {

    require(limma)
    require(Biobase)

    N <- nrow(dat)

    myeset <- ExpressionSet(assayData = as.matrix(dat[, c(group1, group2)]))

    groups <- factor(c(rep("g1", length(group1)),
                       rep("g2", length(group2))))

    design <- model.matrix(~ 0 + groups)
    colnames(design) <- c("g1", "g2")

    contrast.matrix <- makeContrasts("g2-g1", levels = design)

    fit <- lmFit(myeset, design)

    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)

    a <- topTable(fit2, sort="none", n=Inf)

    return (list(dd = dat[, c(group1,  group2)],
                 logFC = a[, "logFC"],
                 pval = ion$impute(a[, "P.Value"], value = 1.0),
                 pval.BH = ion$impute(a[, "adj.P.Val"], value = 1.0)))
}

ion$limma_2g_paired <- function(dat, group1, group2) {

    if (length(group1) != length(group2)) {
        stop("Unequal group for a paired test.")
    }

    N <- nrow(dat)

    pval <- rep(1, N)

    # begin Limma
    require(limma)
    require(Biobase)

    myeset <- ExpressionSet(assayData = as.matrix(dat[, c(group1, group2)]))

    groups <- factor(c(rep("g1", length(group1)),
                       rep("g2", length(group2))))
    patients <- factor(c(1:length(group1), 1:length(group1)))

    design <- model.matrix(~ patients + groups)
    colnames(design)[1] <- "Intercept" # to get rid of the warning

    fit <- lmFit(myeset, design)

    fit2 <- eBayes(fit)

    a <- topTable(fit2, coef="groupsg2", sort="none", n=Inf)

    return (list(dd = dat[, c(group1,  group2)],
                 logFC = a[, "logFC"],
                 pval = ion$impute(a[, "P.Value"], value = 1.0),
                 pval.BH = ion$impute(a[, "adj.P.Val"], value = 1.0)))

}

ion$limma_3g <- function(dat, group1, group2, group3) {

    require(limma)
    require(Biobase)

    N <- nrow(dat)


    myeset <- ExpressionSet(assayData = as.matrix(dat[, c(group1, group2, group3)]))

    groups <- factor(c(rep("g1", length(group1)),
                       rep("g2", length(group2)),
                       rep("g3", length(group3))))

    design <- model.matrix(~ 0 + groups)

    colnames(design) <- c("g1", "g2", "g3")

    contrast.matrix <- makeContrasts(g2-g1, g3-g1, g3-g2, levels = design)

    fit <- lmFit(myeset, design)

    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)

    pval  <- fit2$F.p.value

    return (list(dd = dat[, c(group1,  group2, group3)],
                 pval = ion$impute(pval, value = 1.0),
                 pval.BH = ion$impute(p.adjust(pval, method = "BH"), value = 1.0)))
}

ion$limma_F <- function(dat, groups) {

    require(limma)
    require(Biobase)

    N <- nrow(dat)

    myeset <- ExpressionSet(assayData = as.matrix(dat))

    design <- model.matrix(~ 0 + factor(groups))

    fit <- lmFit(myeset, design)
    fit2 <- eBayes(fit)

    pval  <- fit2$F.p.value

    return (list(dd = dat,
                 pval = ion$impute(pval, value = 1.0),
                 pval.BH = ion$impute(p.adjust(pval, method = "BH"), value = 1.0)))
}

## t-test

ion$t_test <- function(dat, group1, group2, paired = FALSE) {

    if (paired) {
        if (length(group1) != length(group2)) {
            stop("Unequal group for a paired test.")
        }
    }

    N <- nrow(dat)

    pval <- rep(1, N)
    logFC <- rep(0, N)

    for (r in (1:N)) {
        x <- as.numeric(dat[r, group1])
        y <- as.numeric(dat[r, group2])

        # calculate FC in case t-test fails
        logFC[r] <- mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE)
        if (is.nan(logFC[r])) {
            logFC[r] <- 0
        }

        try({
            res <- t.test(x, y = y, paired = paired, na.action=na.omit)
            pval[r] <- res$p.value
            if (paired) {
                logFC[r] <- res$estimate[1]
            }
        }, silent = TRUE)
    }

    pval.BH <- pval

    ind <- pval < 1

    pval.BH[ind] <- p.adjust(pval[ind], method = "BH")

    return (list(dd = dat[, c(group1,  group2)],
                 logFC = logFC,
                 pval = pval,
                 pval.BH = pval.BH))

}

ion$wilcox_test <- function(dat, group1, group2, paired = FALSE) {

    if (paired) {
        if (length(group1) != length(group2)) {
            stop("Unequal group for a paired test.")
        }
    }

    N <- nrow(dat)

    pval <- rep(1, N)
    logFC <- rep(0, N)

    for (r in (1:N)) {
        x <- as.numeric(dat[r, group1])
        y <- as.numeric(dat[r, group2])

        # calculate FC in case the test fails
        logFC[r] <- mean(y - x, na.rm = TRUE)
        if (is.nan(logFC[r])) {
            logFC[r] <- 0
        }

        try({
            res <- wilcox.test(x, y = y, paired = paired, na.action=na.omit, conf.int = TRUE)
            pval[r] <- res$p.value
            if (paired) {
                logFC[r] <- res$estimate[1]
            }
        }, silent = TRUE)
    }

    pval.BH <- pval

    ind <- pval < 1

    pval.BH[ind] <- p.adjust(pval[ind], method = "BH")

    return (list(dd = dat[, c(group1,  group2)],
                 logFC = logFC,
                 pval = pval,
                 pval.BH = pval.BH))

}

## beta-binomial test

ion$fold_change <- function(c1, c2, BIG = 1e4) {

    if (is.vector(c1)) {
        c1 <- matrix(c1, nrow = length(c1))
    }

    if (is.vector(c2)) {
        c2 <- matrix(c2, nrow = length(c2))
    }

    val <- matrix(0, nrow(c1), 1)

    for (i in 1:nrow(c1)) {

        v1 <- sum(c1[i,]) / ncol(c1)

        v2 <- sum(c2[i,]) / ncol(c2)

        if (v1 == 0) {

            if (v2 == 0) {
                val[i] <- 1
            } else {
                val[i] <- BIG
            }
        } else {
            if (v2==0) {
                val[i] <- -BIG
            } else {
                if (v1 > v2) {
                    val[i] <- -v1/v2
                } else {
                    if (v1 < v2) {
                        val[i] <- v2/v1
                    } else {
                        val[i] <- 1
                    }
                }
            }
        }
    }

    return(val)
}

ion$beta_binomial_2g <- function(dat, group1, group2, total_count = NULL) {

    require(ibb)

    d <- dat[, c(group1, group2)]

    out <- bb.test(d,
                   if (is.null(total_count)) colSums(d) else total_count,
                   c(rep("a", length(group1)), rep("b", length(group2))),
                   n.threads = -1)

    d.norm <- ion$normalize_global(d, total_count)

    return (list(fc = ion$fold_change(d.norm[, group1], d.norm[, group2]),
                 pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}

ion$beta_binomial_2g_paired <- function(dat, group1, group2, total_count = NULL, BIG = 1e4) {

    require(ibb)

    d <- dat[, c(group1, group2)]

    out <- ibb.test(d,
                    if (is.null(total_count)) colSums(d) else total_count,
                    c(rep("a", length(group1)), rep("b", length(group2))),
                    n.threads = -1)

    fc <- out$fc

    total_g1 <- if (length(group1) > 1) rowSums(d[, group1]) else d[, group1]
    total_g2 <- if (length(group2) > 1) rowSums(d[, group2]) else d[, group2]

    fc[total_g1 == 0] <- fc[total_g1 == 0]*0 + BIG
    fc[total_g2 == 0] <- fc[total_g2 == 0]*0 - BIG
    fc[total_g1 == 0 & total_g2 == 0] <- 1

    return (list(fc = out$fc,
                 pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}

ion$beta_binomial_3g <- function(dat, group1, group2, group3, total_count = NULL) {

    require(ibb)

    d <- dat[, c(group1, group2, group3)]

    out <- bb.test(d,
                  if (is.null(total_count)) colSums(d) else total_count,
                  c(rep("a", length(group1)), rep("b", length(group2)), rep("c", length(group3))),
                  n.threads = -1)

    return (list(pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}

ion$beta_binomial_4g <- function(dat, group1, group2, group3, group4, total_count = NULL) {

    require(ibb)

    d <- dat[, c(group1, group2, group3, group4)]

    out <- bb.test(d,
                   if (is.null(total_count)) colSums(d) else total_count,
                   c(rep("a", length(group1)), rep("b", length(group2)), rep("c", length(group3)), rep("d", length(group4))),
                   n.threads = -1)

    return (list(pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}

# Text processing ----

ion$str_split <- function(comma_separated_text, sep = ",") {
    return(gsub(" ", "", unlist(strsplit(comma_separated_text, sep))))
}



# Heatmap ----

ion$heatmap <- function(d,
                       # a numeric matrix to show in the heatmap

                       color = c("#4292C6", "#519CCC", "#61A7D2", "#72B1D7", "#85BCDB",
                                 "#98C7DF", "#A9CEE4", "#B8D5EA", "#C6DBEF", "#CFE1F2",
                                 "#D9E7F5", "#E2EDF8", "#EBF3FB", "#F5F9FE", "#F9F9F8",
                                 "#FCF6F1", "#FEF3E8", "#FEEEDE", "#FEE8D3", "#FDE1C4",
                                 "#FDD9B4", "#FDD0A3", "#FDC48F", "#FDB77A", "#FDAA66",
                                 "#FD9E54", "#FD9142", "#FA8432", "#F57622", "#F16913"),
                       # heatmap color, e.g. colorRampPalette(c("white", "Orange"))(32), rev(heat.colors(32)), or bluered(59), or c(colorRampPalette(c("blue", "white"))(30), colorRampPalette(c("white", "red"))(30)[-1])

                       color_min = NULL,
                       # values below this value will get the first color value, default minimum of d_heatmap

                       color_max = NULL,
                       # values above this value will get the last color value, default maximum of d_heatmap

                       z_transform = "row",
                       # set this to "col", "row", or "none".

                       col_data = NA,
                       # data use for column clustering

                       col_distance = "euclidean",
                       # distance for column clustering, e.g. "pearson" = 1-Pearson correlation, "spearman" = 1-Spearman correlation, "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".

                       col_linkage = "complete",
                       # linkage for column clustering, e.g. "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

                       col_reorder = TRUE,
                       # reorder column as input while maintaining the tree structure

                       col_labels = NULL,
                       # column labels, e.g. colnames(d_heatmap),

                       col_label_colors = NULL,
                       # column labels colors

                       col_label_rotated = 90,
                       # degree of labels rotation

                       col_color_bar = NULL,
                       # a list of color vectors, each vector for a bar
                       
                       col_colors = NULL,
                       # a vector of colors

                       row_colors = NULL,
                       # a vector of colors
                       
                       col_margin = 0.5,
                       # margin for column labels

                       row_data = NA,
                       # data for row clustering

                       row_distance = "euclidean",
                       # "pearson", "spearman", "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".

                       row_linkage = "complete",
                       # "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

                       row_reorder = TRUE,
                       # see column

                       row_labels = NULL,
                       # see column, rownames(d_heatmap),

                       row_label_colors = NULL,
                       # see column

                       row_margin = 0.5,
                       # see column

                       sep = FALSE,
                       # a separator between cells

                       main = NULL,
                       cexCol = 1,
                       cexRow = 1,
                       lwid = NULL,
                       lhei = NULL,
                       key = TRUE, ...) {

    zscore <- z_transform != "none"
    
    if (z_transform == "row") {
        cat("\nRows are tranformed into z-values\n\n")
        d_heatmap <- t(scale(t(d)))
    } else {
        if (z_transform == "col") {
            cat("\nColumns are tranformed into z-values\n\n")
            d_heatmap <- scale(d)
        } else {
            d_heatmap <- d
        }
    }
    
    if (!is.null(col_data)) {

        if (identical(col_data, NA)) {
            col_data <- d_heatmap
        }
        
        if (col_distance == "pearson") {
            col_dist <- as.dist(1 - cor(col_data, method = "pearson", use = "pairwise.complete.obs"))
            cat("Using (1 - Pearson correlation) as distance for columns.\n")
        } else if (col_distance == "spearman") {
            col_dist <- as.dist(1 - cor(col_data, method = "spearman", use = "pairwise.complete.obs"))
            cat("Using (1 - Spearman correlation) as distance for columns.\n")
        } else {
            col_dist <- dist(t(col_data), method = col_distance)
            cat("Using", col_distance, "distance for columns.\n")
        }

        if (sum(is.na(col_dist)) > 0) {
            cat("\n=== NA values in col distance matrix ===\n")
            col_dist[is.na(col_dist)]  <-  max(col_dist, na.rm =  TRUE)
        }

        col_tree <- hclust(col_dist, method = col_linkage)
        cat("Using", col_linkage, "linkage for columns.\n")

        if (col_reorder) {
            col_clustering <- reorder(as.dendrogram(col_tree), 1:ncol(as.matrix(col_dist)), agglo.FUN = mean)
        } else {
            col_clustering <- as.dendrogram(col_tree)
        }
    }

    cat("\n")
    
    if (!is.null(row_data)) {
        
        if (identical(row_data, NA)) {
            row_data <- d_heatmap
        }
        
        if (row_distance == "pearson") {
            row_dist <- as.dist(1 - cor(t(row_data), method = "pearson", use = "pairwise.complete.obs"))
            cat("Using (1 - Pearson correlation) as distance for rows.\n")
        } else if (row_distance == "spearman") {
            row_dist <- as.dist(1 - cor(t(row_data), method = "spearman", use = "pairwise.complete.obs"))
            cat("Using (1 - Spearman correlation) as distance for rows.\n")
        } else {
            row_dist <- dist(row_data, method = row_distance)
            cat("Using", row_distance, "distance for rows.\n")
        }

        if (sum(is.na(row_dist)) > 0) {
            cat("\n=== NA values in row distance matrix ===\n")
            row_dist[is.na(row_dist)]  <-  max(row_dist, na.rm =  TRUE)
        }

        row_tree <- hclust(row_dist, method = row_linkage)
        cat("Using", row_linkage, "linkage for rows.\n")

        if (row_reorder) {
            row_clustering <- reorder(as.dendrogram(row_tree), 1:ncol(as.matrix(row_dist)), agglo.FUN = mean)
        }
        else {
            row_clustering <- as.dendrogram(row_tree)
        }
    }
    
    cat("\n")

    colsep  <- seq(1, ncol(d_heatmap))
    rowsep  <- seq(1, nrow(d_heatmap))

    sepcolor  <-  "grey100"
    sepwidth  <-  c(0.02, 0.02)


    if (!sep) {
        colsep <- NULL
        rowsep <- NULL
        sepcolor <- NULL
        sepwidth <- NULL
    }

    col_color_bar_matrix  <-  NULL

    if (!is.null(col_color_bar)) {
        col_color_bar_matrix <- as.matrix(as.data.frame(col_color_bar))
        colnames(col_color_bar_matrix) <- names(col_color_bar)
        if (ncol(col_color_bar_matrix) > 1) {
            col_color_bar_matrix <- col_color_bar_matrix[, ncol(col_color_bar_matrix):1]
        }
    }

    if (is.null(col_data)) {
        if (is.null(row_data)) {
            dendrogram = "none"
        } else {
            dendrogram = "row"
        }
    } else {
        if (is.null(row_data)) {
            dendrogram = "column"
        } else {
            dendrogram = "both"
        }
    }


    ion$.heatmap.2_gplots.3.0.1_modified(as.matrix(d_heatmap),
                                        dendrogram = dendrogram,
                  Colv = if (is.null(col_data)) FALSE else col_clustering,
                  Rowv = if (is.null(row_data)) FALSE else row_clustering,
                  col = color,
                  breaks = if (zscore) seq(-2, 2, length.out = (length(color) + 1)) else seq(ifelse(is.null(color_min), min(d_heatmap, na.rm = TRUE), color_min), ifelse(is.null(color_max), max(d_heatmap, na.rm = TRUE), color_max), length.out = (length(color) + 1)),
                  key = key,
                  keysize = 1,

                  cexRow = cexRow,
                  cexCol = cexCol,

                  colsep = colsep,
                  rowsep = rowsep,
                  sepcolor = sepcolor,
                  sepwidth = sepwidth,

                  scale = "none",
                  symkey = zscore,
                  density.info = "none",
                  trace = "none",
                  margins = c(col_margin, row_margin),

                  labRow = if (is.null(row_labels)) "" else row_labels,
                  labCol = if (is.null(col_labels)) "" else col_labels,

                  ColSideColors = if (is.null(col_colors)) col_color_bar_matrix else col_colors,
                  
                  RowSideColors = row_colors,
                  
                  colRow = row_label_colors,
                  colCol = col_label_colors,

                  srtCol = col_label_rotated,

                  lwid = lwid,
                  lhei = lhei,
                  main = main, ...)
}

ion$.heatmap.2_gplots.3.0.1_modified <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
    distfun = dist, hclustfun = hclust, dendrogram = c("both",
        "row", "column", "none"), reorderfun = function(d, w) reorder(d,
        w), symm = FALSE, scale = c("none", "row", "column"),
    na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks,
    symbreaks = any(x < 0, na.rm = TRUE) || scale != "none",
    col = "heat.colors", colsep, rowsep, sepcolor = "white",
    sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan",
    na.color = par("bg"), trace = c("column", "row", "both",
        "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks),
    linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors,
    cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
    labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0,
        NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5,
    colRow = NULL, colCol = NULL, key = TRUE, keysize = 1.5,
    density.info = c("histogram", "density", "none"), denscol = tracecol,
    symkey = any(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25,
    key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL,
    key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL,
    ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL,
    # Thang
    key_margins = NULL, ...)
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && any(duplicated(breaks)))
        stop("breaks may not contian duplicate values")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (all(Colv == "Rowv"))
        Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) &&
            (dendrogram %in% c("both", "row"))) {
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
            if (dendrogram == "both")
                dendrogram <- "column"
            else dendrogram <- "none"
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) &&
            (dendrogram %in% c("both", "column"))) {
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
            if (dendrogram == "both")
                dendrogram <- "row"
            else dendrogram <- "none"
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd >
            nr))
            stop("Rowv dendrogram doesn't match size of x")
        if (length(rowInd) < nr)
            nr <- length(rowInd)
    }
    else if (is.integer(Rowv)) {
        distr <- distfun(x)
        hcr <- hclustfun(distr)
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        distr <- distfun(x)
        hcr <- hclustfun(distr)
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (!isTRUE(Rowv)) {
        rowInd <- nr:1
        ddr <- as.dendrogram(hclust(dist(diag(nr))))
    }
    else {
        rowInd <- nr:1
        ddr <- as.dendrogram(Rowv)
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd >
            nc))
            stop("Colv dendrogram doesn't match size of x")
        if (length(colInd) < nc)
            nc <- length(colInd)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        distc <- distfun(if (symm)
            x
        else t(x))
        hcc <- hclustfun(distc)
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        distc <- distfun(if (symm)
            x
        else t(x))
        hcc <- hclustfun(distc)
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (!isTRUE(Colv)) {
        colInd <- 1:nc
        ddc <- as.dendrogram(hclust(dist(diag(nc))))
    }
    else {
        colInd <- 1:nc
        ddc <- as.dendrogram(Colv)
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (!is.null(colRow))
        colRow <- colRow[rowInd]
    if (!is.null(colCol))
        colCol <- colCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
        1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        #Thang BEGIN
        if (!is.null(ColSideColors)) {

            if (!is.matrix(ColSideColors)){
                if (!is.character(ColSideColors) || length(ColSideColors) != nc)
                    stop("'ColSideColors' must be a character vector of length ncol(x)")
                lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
                              1)
                lhei <- c(lhei[1], 0.2, lhei[2])
            } else {
                if (!is.character(ColSideColors) || dim(ColSideColors)[1] != nc)
                #if (dim(ColSideColors)[1] != nc)
                    stop("'ColSideColors' dim()[2] must be of length ncol(x)")
                lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
                lhei <- c(lhei[1], 0.2 * dim(ColSideColors)[2], lhei[2])

            }
        }

        #if (!missing(ColSideColors)) {
        #    if (!is.character(ColSideColors) || length(ColSideColors) !=
        #        nc)
        #        stop("'ColSideColors' must be a character vector of length ncol(x)")
        #    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
        #        1)
        #    lhei <- c(lhei[1], 0.2, lhei[2])
        #}
        
        # Thang END
        
        # Thang
        #if (!missing(RowSideColors)) {
        if (!is.null(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) !=
                nr)
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    plot.index <- 1

    # Thang
    #if (!missing(RowSideColors)) {
    if (!is.null(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        plot.index <- plot.index + 1
    }
    #if (!missing(ColSideColors)) {
    #    par(mar = c(0.5, 0, 0, margins[2]))
    #    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    #    plot.index <- plot.index + 1
    #}

# Thang BEGIN
if (!is.null(ColSideColors)) {
if (!is.matrix(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    #image(cbind(1:nc), col = ColSideColors[colInd, 1], axes = FALSE)
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    plot.index <- plot.index + 1
}
else{
    par(mar = c(0.5, 0, 0, margins[2]))
    csc = matrix(ColSideColors[colInd, ], nrow = length(colInd))
    csc.colors = matrix()
    csc.names = names(table(csc))
    csc.i = 1
    for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
    }
    csc = matrix(as.numeric(csc), nrow = dim(csc)[1])

a <- as.vector(csc.colors)
a[nchar(a) <= 1] <- "white"
a[startsWith(a, "`")] <- "white"

image(csc, col = a, axes = F)

a <- as.vector(csc.colors)
for (i in 1:nrow(csc)) {
    for (j in 1:ncol(csc)) {
        if (nchar(a[csc[i,j]]) <= 1) {
            text((i-1)/(nrow(csc)-1), (j-1)/(ncol(csc)-1),a[csc[i,j]], cex = 1.5)
        }
        if (startsWith(a[csc[i,j]],"`")) {
            eval(parse(text=paste0("points(",
                                   (i-1)/(nrow(csc)-1),
                                   ",",
                                   (j-1)/(ncol(csc)-1), 
                                   ",",
                                   substring(a[csc[i,j]],2,1000),
                                   ")")))
        }
    }
}
#text(0.1, 0.1, "+")
#text(1, 1, "++")

    abline(h=(0:(dim(csc)[2]-2))/(dim(csc)[2] - 1) + (0.5/(dim(csc)[2] - 1)), col="white", lwd=5);

    if (length(colnames(ColSideColors)) > 0) {
        if (dim(csc)[2] > 1) {
            axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1),
                 colnames(ColSideColors),
                 las = 2, tick = FALSE)
        } else {
            axis(2, 0,
                 colnames(ColSideColors),
                 las = 2, tick = FALSE)
        }
    }
    plot.index <- plot.index + 1
}
}
    # Thang END

    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col

# Thang
#    if (!invalid(na.color) & any(is.na(x))) {

    if (any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol) && is.null(colCol))
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 +
            offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1],
            padj = adjCol[2])
    else {
        if (is.null(srtCol) || is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol))
                adjCol = c(1, NA)
            if (is.null(srtCol))
                srtCol <- 90
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2,
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) *
                strheight("M"), labels = labCol, adj = adjCol,
                cex = cexCol, srt = srtCol, col = colCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow) && is.null(colRow)) {
        axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow,
            tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
    }
    else {
        if (is.null(srtRow) || is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2,
                line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
                y = ypos, labels = labRow, adj = adjRow, cex = cexRow,
                srt = srtRow, col = colRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0,
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) +
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in 1:length(colInd)) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in 1:length(rowInd)) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol,
                  lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    plot.index <- plot.index + 1
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
# Thang
#        flag <- try(plot.dendrogram(ddr, horiz = TRUE, axes = FALSE,
        #                                    yaxs = "i", leaflab = "none"))
        flag <- try(plot(ddr, horiz = TRUE, axes = FALSE,
                         yaxs = "i", leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
                stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
# Thang
#        flag <- try(plot.dendrogram(ddc, axes = FALSE, xaxs = "i",
        #            leaflab = "none"))
        flag <- try(plot(ddc, axes = FALSE, xaxs = "i",
                         leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
                stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        mar <- c(5, 4, 2, 1)
        if (!is.null(key.xlab) && is.na(key.xlab))
            mar[1] <- 2
        if (!is.null(key.ylab) && is.na(key.ylab))
            mar[2] <- 2
        if (!is.null(key.title) && is.na(key.title))
            mar[3] <- 1
        par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
        if (length(key.par) > 0)
            do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            
            #max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            #min.raw <- -max.raw
            #tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            #tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)

            # Thang 180628 - unsorted break complain
            min.raw <- tmpbreaks[1] <- -max(abs(c(x, breaks)), na.rm = TRUE)
            max.raw <- tmpbreaks[length(tmpbreaks)] <- -min.raw
        }
        else {
            min.raw <- min.breaks
            max.raw <- max.breaks
        }
        z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
        
        
        # Thang BEGIN
        if (is.null(key_margins)) {
            tmp <- dev.size("in")
            key_margins <- c(tmp[1]/2, 1, 0.5, 2)
        }
        
        par(mar = key_margins)
        
        do_stop <- FALSE
        tryCatch({
            image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
                  xaxt = "n", yaxt = "n")},
            warning = function(w) {
                do_stop <<- TRUE
            },
            error = function (e) {
                do_stop <<- TRUE
            })
        if (do_stop) {
            cat("\nCannot make legend key.\nRun dev.off(), then make the margin bigger or set key = FALSE.\n")
            return(invisible(0))
        }
        # Thang END
        
        par(usr = c(0, 1, 0, 1))
        
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row")
                key.xlab <- "Row Z-Score"
            else if (scale == "column")
                key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5,
                cex = par("cex") * par("cex.lab"))
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE,
                from = min.scale, to = max.scale)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[!omit]
            dens$y <- dens$y[!omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(dens$y)/max(dens$y) *
                  0.95, labels = pretty(dens$y))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title))
                key.title <- "Color Key\nand Density Plot"
            if (!is.na(key.title))
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab))
                key.ylab <- "Density"
            if (!is.na(key.ylab))
                mtext(side = 2, key.ylab, line = par("mgp")[1],
                  padj = 0.5, cex = par("cex") * par("cex.lab"))
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(hy)/max(hy) * 0.95,
                  labels = pretty(hy))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title))
                key.title <- "Color Key\nand Histogram"
            if (!is.na(key.title))
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab))
                key.ylab <- "Count"
            if (!is.na(key.ylab))
                mtext(side = 2, key.ylab, line = par("mgp")[1],
                  padj = 0.5, cex = par("cex") * par("cex.lab"))
        }
        else if (is.null(key.title)) {
            # Thang
            #title("Color Key")
        }
        if (trace %in% c("both", "column")) {
            vline.vals <- scale01(vline, min.raw, max.raw)
            if (!is.null(vline)) {
                abline(v = vline.vals, col = linecol, lty = 2)
            }
        }
        if (trace %in% c("both", "row")) {
            hline.vals <- scale01(hline, min.raw, max.raw)
            if (!is.null(hline)) {
                abline(v = hline.vals, col = linecol, lty = 2)
            }
        }
    } else {
        par(mar = c(0, 0, 0, 0))
        plot.new()
    }
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    retval$layout <- list(lmat = lmat, lhei = lhei, lwid = lwid)
    if (!is.null(extrafun))
        extrafun()
    invisible(retval)
}