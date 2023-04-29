# Author: Thang V Pham
#
# Copyright Thang V Pham, 2018-2022
#
# Version: 22.02.04
#
# NOTE: sourcing this will replace any existing 'ion' object

ion <- list()

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


                        #--- COLUMNS ---

                        col_data = NA,
                        # data for column clustering, use d if NA

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

                        #--- ROW ---

                        row_data = NA,
                        # data for row clustering, use d if NA

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


                        #--- OTHERS ---

                        color_key_margins = NULL,
                        # margins of the color key, default c(5, 1, 0.5, 2) if NULL
                        
                        key_margins,
                        # not used

                        separator = FALSE,
                        # a separator between cells

                        ...) {

    if (!missing(key_margins)) {
        message("Please use 'color_key_margins' instead of 'key_margins'.")
        return(invisible(NULL))
    }
    
    zscore <- z_transform != "none"

    if (z_transform == "row") {
        message("Rows are tranformed into z-values\n")
        d_heatmap <- t(scale(t(d)))
    } else {
        if (z_transform == "col") {
            message("Columns are tranformed into z-values\n")
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
            message("Using (1 - Pearson correlation) as distance for columns.")
        } else if (col_distance == "spearman") {
            col_dist <- as.dist(1 - cor(col_data, method = "spearman", use = "pairwise.complete.obs"))
            message("Using (1 - Spearman correlation) as distance for columns.")
        } else {
            col_dist <- dist(t(col_data), method = col_distance)
            message("Using ", col_distance, " distance for columns.")
        }

        if (sum(is.na(col_dist)) > 0) {
            message("\n=== NA values in col distance matrix ===\n")
            col_dist[is.na(col_dist)]  <-  max(col_dist, na.rm =  TRUE)
        }

        col_tree <- hclust(col_dist, method = col_linkage)
        message("Using ", col_linkage, " linkage for columns.\n")

        if (col_reorder) {
            col_clustering <- reorder(as.dendrogram(col_tree), 1:ncol(as.matrix(col_dist)), agglo.FUN = mean)
        } else {
            col_clustering <- as.dendrogram(col_tree)
        }
    }


    if (!is.null(row_data)) {

        if (identical(row_data, NA)) {
            row_data <- d_heatmap
        }

        if (row_distance == "pearson") {
            row_dist <- as.dist(1 - cor(t(row_data), method = "pearson", use = "pairwise.complete.obs"))
            message("Using (1 - Pearson correlation) as distance for rows.")
        } else if (row_distance == "spearman") {
            row_dist <- as.dist(1 - cor(t(row_data), method = "spearman", use = "pairwise.complete.obs"))
            message("Using (1 - Spearman correlation) as distance for rows.")
        } else {
            row_dist <- dist(row_data, method = row_distance)
            message("Using ", row_distance, " distance for rows.")
        }

        if (sum(is.na(row_dist)) > 0) {
            message("\n=== NA values in row distance matrix ===\n")
            row_dist[is.na(row_dist)]  <-  max(row_dist, na.rm =  TRUE)
        }

        row_tree <- hclust(row_dist, method = row_linkage)
        message("Using ", row_linkage, " linkage for rows.")

        if (row_reorder) {
            row_clustering <- reorder(as.dendrogram(row_tree), 1:ncol(as.matrix(row_dist)), agglo.FUN = mean)
        }
        else {
            row_clustering <- as.dendrogram(row_tree)
        }
    }

    colsep  <- seq(1, ncol(d_heatmap))
    rowsep  <- seq(1, nrow(d_heatmap))

    sepcolor  <-  "grey100"
    sepwidth  <-  c(0.02, 0.02)


    if (!separator) {
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

                                         key_margins = if (is.null(color_key_margins)) c(dev.size("in")[2]+0.5, 1, 0.5, 2) else color_key_margins,

                                         ...)
}

ion$plot_pca <- function(d, 
                         pch = 19, 
                         main = "PCA plot", 
                         xlab = "PC 1", 
                         ylab = "PC 2", ...) {
    
    pca <- prcomp(t(d), center = TRUE, scale = TRUE)
    
    projection <- predict(pca, newdata = t(d))
    
    plot(projection[, c("PC1", "PC2")], pch = pch, main = main, 
         xlab = xlab, ylab = ylab, ...)
}


ion$plot_volcano <- function (logFC, pval, 
                              main = "", 
                              pval_thresh = 0.05, 
                              fc_thres = log2(2.0), 
                          xlab = "Log2 fold change") {
    x <- logFC
    x[is.na(x)] <- 0
    
    y <- -log10(pval)
    log_pval_thresh <- -log10(pval_thresh)
    
    point_color <- rep("gray10",length(x))
    
    r_significant <- y > log_pval_thresh
    r_up   <- x > fc_thres
    r_down <- x < -fc_thres
    
    up_color <- "red"
    down_color <- "blue"
    
    ret <- r_significant & (r_up | r_down)
    ret[is.na(ret)] <- FALSE
    
    point_color[ r_significant & r_up ] <- up_color
    point_color[ r_significant & r_down ] <- down_color
    
    x_width <- max(abs(x)) * 1.2
    
    plot(x, y, main = main, xlab = xlab,
         pch = 20,
         ylim = c(0, max(y) * 1.2),
         xlim = c(-x_width, x_width),
         ylab = "-log10 p-value",
         col = point_color)
    
    abline(h = log_pval_thresh, lty=2, col="gray")
    
    text(1 - x_width, log_pval_thresh,
         labels = paste0("pval=", sprintf("%.2f",pval_thresh)), pos=3)
    
    abline(v = fc_thres, lty = 2, col = "gray")
    abline(v = -fc_thres, lty=2, col = "gray")
    
    legend("topright", border = up_color, col = up_color,
           legend= paste0("UP: ", sum(r_significant & r_up )), bty = "n", pch = 20, pt.cex = 3)
    
    legend("topleft", border = down_color, col = down_color,
           legend= paste0("DOWN: ", sum(r_significant & r_down )), bty = "n", pch = 20, pt.cex = 3)
    
    return(ret)
}


# Normalization ----

ion$normalize_global <- function(d, total_count = NULL) {

    total <- if (is.null(total_count)) apply(d, 2, sum) else total_count

    m <- mean(total)
    factor <- total / m

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
        } else if (method == "min") {
            global_min <- min(d, na.rm=TRUE)
            dd[is.na(dd)] <- global_min
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
                     col_names = NULL,
                     stringsAsFactors = FALSE,
                     ...
                     ) {
    if (is.null(col_names)) {
        return(read.delim(filename, quote = "", stringsAsFactors = stringsAsFactors, ...))
    } else {
        h <- read.delim(filename, nrow = 1)

        a <- setdiff(col_names, colnames(h))

        if (length(a) > 0) {
            message(paste("Column not available:", a, "\n"))
            return(NULL)
        }

        v <- rep("NULL", ncol(h))
        names(v) <- colnames(h)
        v[col_names] <- NA

        return(read.delim(filename, colClasses = v, stringsAsFactors = stringsAsFactors))
    }
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

ion$limma_4g <- function(dat, group1, group2, group3, group4) {
    
    require(limma)
    require(Biobase)
    
    N <- nrow(dat)
    
    
    myeset <- ExpressionSet(assayData = as.matrix(dat[, c(group1, group2, group3, group4)]))
    
    groups <- factor(c(rep("g1", length(group1)),
                       rep("g2", length(group2)),
                       rep("g3", length(group3)),
                       rep("g4", length(group4))))
    
    design <- model.matrix(~ 0 + groups)
    
    colnames(design) <- c("g1", "g2", "g3", "g4")
    
    contrast.matrix <- makeContrasts(g2-g1, g3-g1, g4-g1, g3-g2, g4-g2, g4-g3, levels = design)
    
    fit <- lmFit(myeset, design)
    
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    pval  <- fit2$F.p.value
    
    return (list(dd = dat[, c(group1,  group2, group3, group4)],
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

## anova
ion$anova_subject_time <- function(mat, subjects, time_points, random_effect = FALSE) {

    pval <- rep(NA, nrow(mat))

    if (random_effect) {

        library(lme4)

        for (i in 1:nrow(mat)) {
            try({
                ds = list(y = as.numeric(mat[i, ]), subjects = subjects , time_points = time_points)

                m1 <- lmer(y ~ time_points + (1|subjects), ds, REML = FALSE)
                m2 <- lmer(y ~ (1|subjects), ds, REML = FALSE)

                pval[i] <- anova(m2, m1, test = "Chisq")[2, "Pr(>Chisq)"]

            }, silent = TRUE)
        }
    } else {
        for (i in 1:nrow(mat)) {
            try({
                ds = list(y = as.numeric(mat[i, ]), subjects = subjects , time_points = time_points)

                m1 <- lm(y ~ subjects + time_points, data = ds)

                m2 <- lm(y ~ subjects, data = ds)

                pval[i] <- anova(m2, m1, test = "Chisq")[2, "Pr(>Chi)"]
            }, silent = TRUE)
        }
    }

    p.BH <- pval
    p.BH[!is.na(pval)] <- p.adjust(pval[!is.na(pval)], method = "BH")

    return (list(dd = mat,
                 pval = ion$impute(pval, value = 1.0),
                 pval.BH = ion$impute(p.BH, value = 1.0)))
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
        logFC[r] <- median(y, na.rm = TRUE) - median(x, na.rm = TRUE)
        if (is.na(logFC[r]) || is.nan(logFC[r])) {
            logFC[r] <- 0
        }

        try({
            if (paired) {
                res <- wilcox.test(x, y = y, paired = TRUE, na.action=na.omit, conf.int = TRUE)
                pval[r] <- res$p.value
                logFC[r] <- res$estimate[1]
            } else {
                res <- wilcox.test(x, y = y, na.action=na.omit)
                pval[r] <- res$p.value    
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

    require(countdata)

    d <- dat[, c(group1, group2)]

    out <- countdata::bb.test(d,
                              if (is.null(total_count)) colSums(d) else total_count,
                              c(rep("a", length(group1)), rep("b", length(group2))),
                              n.threads = -1)

    d.norm <- ion$normalize_global(d, total_count)

    return (list(fc = ion$fold_change(d.norm[, group1], d.norm[, group2]),
                 pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}

ion$beta_binomial_2g_paired <- function(dat, group1, group2, total_count = NULL, BIG = 1e4) {

    require(countdata)

    d <- dat[, c(group1, group2)]

    out <- countdata::ibb.test(d,
                               if (is.null(total_count)) colSums(d) else total_count,
                               c(rep("a", length(group1)), rep("b", length(group2))),
                               n.threads = -1, BIG = BIG)

    return (list(fc = out$fc,
                 pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}

ion$beta_binomial_3g <- function(dat, group1, group2, group3, total_count = NULL) {

    require(countdata)

    d <- dat[, c(group1, group2, group3)]

    out <- countdata::bb.test(d,
                              if (is.null(total_count)) colSums(d) else total_count,
                              c(rep("a", length(group1)), rep("b", length(group2)), rep("c", length(group3))),
                              n.threads = -1)

    return (list(pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}

ion$beta_binomial_4g <- function(dat, group1, group2, group3, group4, total_count = NULL) {

    require(countdata)

    d <- dat[, c(group1, group2, group3, group4)]

    out <- countdata::bb.test(d,
                              if (is.null(total_count)) colSums(d) else total_count,
                              c(rep("a", length(group1)), rep("b", length(group2)), rep("c", length(group3)), rep("d", length(group4))),
                              n.threads = -1)

    return (list(pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}

ion$beta_binomial_5g <- function(dat, group1, group2, group3, group4, group5, total_count = NULL) {
    
    require(countdata)
    
    d <- dat[, c(group1, group2, group3, group4, group5)]
    
    out <- countdata::bb.test(d,
                              if (is.null(total_count)) colSums(d) else total_count,
                              c(rep("a", length(group1)), rep("b", length(group2)), rep("c", length(group3)), 
                                rep("d", length(group4)), rep("e", length(group5))),
                              n.threads = -1)
    
    return (list(pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}


ion$beta_binomial_6g <- function(dat, group1, group2, group3, group4, group5, group6, total_count = NULL) {
    
    require(countdata)
    
    d <- dat[, c(group1, group2, group3, group4, group5, group6)]
    
    out <- countdata::bb.test(d,
                              if (is.null(total_count)) colSums(d) else total_count,
                              c(rep("a", length(group1)), rep("b", length(group2)), rep("c", length(group3)), 
                                rep("d", length(group4)), rep("e", length(group5)), rep("f", length(group6))),
                              n.threads = -1)
    
    return (list(pval = out$p.value,
                 pval.BH = p.adjust(out$p.value, method = "BH")))
}


ion$cuzick_test <- function(x, z, test.type=c("two.sided", "upper", "lower")) {
    
    # https://r.789695.n4.nabble.com/Cuzick-s-test-for-trend-td807101.html
    
    N = length(z)
    n = unique(z)
    
    ranks=rank(x)
    
    T = sum(ranks*z)
    
    p = (table(z)/N)
    E_Z = sum(unique(z)*p)
    E_T = 0.5*N*(N+1)*E_Z
    
    Var_Z = sum(unique(z)^2*p) - E_Z^2
    Var_T = N^2*(N+1)/12*Var_Z
    
    Zscore = (T-E_T)/sqrt(Var_T)
    
    if(test.type == "two.sided") {
        pval = 2*pnorm(-abs(Zscore))
    } else if(test.type == "upper") {
        pval = pnorm(Zscore,lower.tail=F)
    } else pval = pnorm(Zscore,lower.tail=T)
    
    out = data.frame(cbind(Zscore,pval,test.type))
    colnames(out) = c("Z","p","testType")
    #return(out)
    return(pval)
}

# Text processing ----

ion$str_split <- function(comma_separated_text, sep = ",") {
    return(unlist(strsplit(comma_separated_text, sep)))
}

ion$str_split_trim <- function(comma_separated_text, sep = ",") {
    return(gsub(" ", "", unlist(strsplit(comma_separated_text, sep))))
}

# Machine learning ----

ion$auc <- function(y, yhat, make_plot = TRUE, ...) {
    require(pROC)
    r <- roc(controls = yhat[y == 1], 
             cases = yhat[y != 1], quiet = TRUE)
    if (make_plot) {
        plot(r, ...)
    }
    return(r$auc)
}


ion$find_max_auc <- function(y, X) {
    require(pROC)
    # return column with maximum AUC
    a <- 0
    col_max <- -1
    for (i in 1:ncol(X)) {
        b <- ion$auc(y, X[,i], make_plot = FALSE)
        if (b > a) {
            a <- b
            col_max <- i
        }
    }
    return(col_max)
}


ion$get_folds <- function(y, nfold = 10) {
    
    g <- unique(y)
    
    a <- which(y == g[1])
    b <- which(y == g[2])
    
    a_random <- sample(a)
    b_random <- sample(b)
    
    folds <- vector(mode="list", length = nfold)
    
    for (i in 1:length(a_random)) {
        j <- (i %% nfold) + 1
        folds[[j]] <- c(folds[[j]], a_random[i])
    }
    
    for (i in 1:length(b_random)) {
        j <- (i %% nfold) + 1
        folds[[j]] <- c(folds[[j]], b_random[i])
    }
    
    return(folds)
}

ion$regression_form <- function(vars, is_quadratic = FALSE, env = environment()) {
    
    if (is_quadratic) {
        ret <- paste0(vars[1], " ~ ", vars[2], " + I(", vars[2], "^2)")
        
        if (length(vars) > 2) {
            for (i in 3:length(vars)) {
                ret <- paste0(ret, " + ", vars[i] , " + I(", vars[i], "^2)")
            }
            
            for (i in 2:(length(vars)-1)) {
                for (j in (i+1):length(vars)) {
                    ret <- paste0(ret, " + ", vars[i], ":", vars[j])
                }
            } 
        }
    } else {
        ret <- paste0(vars[1], " ~ ", vars[2])
        if (length(vars) > 2) {
            for (i in 3:length(vars)) {
                ret <- paste0(ret, " + ", vars[i])
            }
        }
    }
    
    return(as.formula(ret, env = env))
}

ion$logistic_regression_LOOCV <- function(y, X, is_quadratic = FALSE) {
    
    if (nrow(X) != length(y)) {
        stop("The number of samples (X, y) does not match.\n")
    }
    
    v <- rep(0, nrow(X))
    
    f <- colnames(X)
    
    form <- ion$regression_form(c("y", f), is_quadratic = is_quadratic)
    
    for (r in 1:nrow(X)) {
        
        d1 <-  data.frame(y = y[-r], X[-r, , drop = FALSE])
        
        logistic_model <- glm(form,
                              data = d1,
                              family = "binomial")
        
        v[r] <- predict(logistic_model,
                        newdata = split(X[r, ], f),
                        type = "response")
    }
    
    return(v)
}


ion$logistic_regression_bruteforce_4 <- function(y, X, alpha_top = 0.01, n_threads = 1) {
    
    require(pROC)    
    
    message("Checking ", choose(ncol(X), 4), " combinations.")
    # brute-force
    res <- list()
    cc <- 0
    
    N <- ncol(X)
    
    for (i1 in 1:(N-3)) { 
        for (i2 in (i1+1):(N-2)) {
            for (i3 in (i2+1):(N-1)) {
                for (i4 in (i3+1):N) {
                    cc <- cc + 1
                    res[[cc]] <- c(i1, i2, i3, i4)
                }
            }
        }
    }
    
    message(length(res), " combinations.")
    
    vec_auc <- rep(NA, length(res))
    
    if (n_threads == 1) {
        for (i in 1:length(res)) {
            if (i %% 10000 == 0) {
                message(i, "/", length(res), " = ", floor(i * 100/length(res)), "%")
            }
            v <- ion$logistic_regression_LOOCV(y, X[, res[[i]]], is_quadratic = FALSE)
            r <- roc(controls = v[y == 1], 
                     cases = v[y != 1],
                     auc = TRUE, quiet = TRUE)
            vec_auc[i] <- r$auc
        }
    } else {
        require(parallel)
        
        nt <- ifelse(n_threads > 1, n_threads,
                     parallel::detectCores(logical = FALSE) + n_threads)
        
        message("Using ", nt, " core(s)")
        
        cl <- makeCluster(nt)
        clusterExport(cl=cl, varlist=c("y", "X", "roc", "ion"),  envir=environment())
        
        vec_auc <- parSapply(cl, res, function(combi) {
            v <- ion$logistic_regression_LOOCV(y, X[, combi], is_quadratic = FALSE)
            r <- roc(controls = v[y == 1], 
                     cases = v[y != 1],
                     auc = TRUE, quiet = TRUE)
            r$auc
        })
        stopCluster(cl)
    }
    
    ind <- order(vec_auc, decreasing = TRUE)
    
    N <- floor(length(res) * alpha_top)
    
    message("Accumulate variables in top ", N, " out of ", length(res), " combinations.")
    
    h <- rep(0, ncol(X))
    names(h) <- colnames(X)
    for (i in 1:N) {
        h[res[[ind[i]]]] <- h[res[[ind[i]]]] + 1 
    }
    ih <- order(h, decreasing = TRUE)
    
    return(list(all_combinations = res, vec_auc = vec_auc, frequency_top = h[ih]))
    
}

ion$logistic_regression_bruteforce_5 <- function(y, X, alpha_top = 0.01, n_threads = 1) {
    
    require(pROC)    
    
    message("Checking ", choose(ncol(X), 5), " combinations.")
    # brute-force
    res <- list()
    cc <- 0
    
    N <- ncol(X)
    
    for (i1 in 1:(N-4)) { 
        for (i2 in (i1+1):(N-3)) {
            for (i3 in (i2+1):(N-2)) {
                for (i4 in (i3+1):(N-1)) {
                    for (i5 in (i4+1):N) {
                        cc <- cc + 1
                        res[[cc]] <- c(i1, i2, i3, i4, i5)
                    }
                }
            }
        }
    }
    
    message(length(res), " combinations.")
    
    vec_auc <- rep(NA, length(res))
    
    if (n_threads == 1) {
        for (i in 1:length(res)) {
            if (i %% 10000 == 0) {
                message(i, "/", length(res), " = ", floor(i * 100/length(res)), "%")
            }
            v <- ion$logistic_regression_LOOCV(y, X[, res[[i]]], is_quadratic = FALSE)
            r <- roc(controls = v[y == 1], 
                     cases = v[y != 1],
                     auc = TRUE, quiet = TRUE)
            vec_auc[i] <- r$auc
        }
    } else {
        require(parallel)
        
        nt <- ifelse(n_threads > 1, n_threads,
                     parallel::detectCores(logical = FALSE) + n_threads)
        
        message("Using ", nt, " core(s)")
        
        cl <- makeCluster(nt)
        clusterExport(cl=cl, varlist=c("y", "X", "roc", "ion"),  envir=environment())
        
        vec_auc <- parSapply(cl, res, function(combi) {
            v <- ion$logistic_regression_LOOCV(y, X[, combi], is_quadratic = FALSE)
            r <- roc(controls = v[y == 1], 
                     cases = v[y != 1],
                     auc = TRUE, quiet = TRUE)
            r$auc
        })
        stopCluster(cl)
    }
    
    ind <- order(vec_auc, decreasing = TRUE)
    
    N <- floor(length(res) * alpha_top)
    
    message("Accumulate variables in top ", N, " out of ", length(res), " combinations.")
    
    h <- rep(0, ncol(X))
    names(h) <- colnames(X)
    for (i in 1:N) {
        h[res[[ind[i]]]] <- h[res[[ind[i]]]] + 1 
    }
    ih <- order(h, decreasing = TRUE)
    
    return(list(all_combinations = res, vec_auc = vec_auc, frequency_top = h[ih]))
    
}

ion$logistic_regression_random_combination <- function(y, X, k, n_samples, alpha_top = 0.01, n_threads = 1, seed = 1) {
    
    require(pROC)    
    
    message("Checking ", n_samples, " out of ", choose(ncol(X), k), " combinations.")
    
    set.seed(seed = seed)
    
    N <- ncol(X)
    
    res <- list()
    for (i in 1:n_samples) { 
        res[[i]] <- sample(N, k)
    }
    
    vec_auc <- rep(NA, length(res))
    
    if (n_threads == 1) {
        for (i in 1:length(res)) {
            if (i %% 10000 == 0) {
                message(i, "/", length(res), " = ", floor(i * 100/length(res)), "%")
            }
            v <- ion$logistic_regression_LOOCV(y, X[, res[[i]]], is_quadratic = FALSE)
            r <- roc(controls = v[y == 1], 
                     cases = v[y != 1],
                     auc = TRUE, quiet = TRUE)
            vec_auc[i] <- r$auc
        }
    } else {
        require(parallel)
        
        nt <- ifelse(n_threads > 1, n_threads,
                     parallel::detectCores(logical = FALSE) + n_threads)
        
        message("Using ", nt, " core(s)")
        
        cl <- makeCluster(nt)
        clusterExport(cl=cl, varlist=c("y", "X", "roc", "ion"),  envir=environment())
        
        vec_auc <- parSapply(cl, res, function(combi) {
            v <- ion$logistic_regression_LOOCV(y, X[, combi], is_quadratic = FALSE)
            r <- roc(controls = v[y == 1], 
                     cases = v[y != 1],
                     auc = TRUE, quiet = TRUE)
            r$auc
        })
        stopCluster(cl)
    }
    
    ind <- order(vec_auc, decreasing = TRUE)
    
    N <- floor(length(res) * alpha_top)
    
    message("Accumulate variables in top ", N, " out of ", length(res), " combinations.")
    
    h <- rep(0, ncol(X))
    names(h) <- colnames(X)
    for (i in 1:N) {
        h[res[[ind[i]]]] <- h[res[[ind[i]]]] + 1 
    }
    ih <- order(h, decreasing = TRUE)
    
    return(list(all_combinations = res, vec_auc = vec_auc, frequency_top = h[ih]))
    
}

ion$logistic_regression_thresholded_predict <- function(str_model, X) {
    
    Xnew <- data.frame(X)
    if (!identical(colnames(Xnew), colnames(str_model$X))) {
        stop("Variables in data and model do not match.")
    }
    
    for (i in 1:ncol(Xnew)) {
        Xnew[, i] <- as.numeric(Xnew[, i] > str_model$threshold_values[i])
    }
    
    return(predict(str_model$logistic_model,
                   newdata = Xnew,
                   type = "response"))
}

ion$logistic_regression_thresholded_LOO_AUC <- function(y, X, max_refine_round = 100, 
                                                n_quantile = 4,
                                                pdf_out = "output-str") {
    
    quantile_cutoff <- seq(0, 1, 1 / n_quantile)
    
    require(pROC)
    
    Xnew <- as.matrix(X)
    for (i in 1:ncol(Xnew)) {
        Xnew[, i] <- Xnew[, i] > median(Xnew[, i]) 
    }
    
    best_threshold <- NULL
        
    for (k in 1:max_refine_round) {
            
        best_thres <- rep(0, ncol(Xnew))
            
        for (i in 1:ncol(Xnew)) {
                
            thres <- quantile(X[, colnames(Xnew)[i]], probs = quantile_cutoff)
                
            Xnew2 <- Xnew
            max_val <- -Inf
                
            for (j in 2:(length(thres)-1)) {
                    
                Xnew2[, i] <- as.numeric(X[, colnames(Xnew)[i]] > thres[j])
                    
                v <- ion$logistic_regression_LOOCV(y, Xnew2)
                    
                r <- roc(controls = v[y == 1], 
                         cases = v[y != 1],
                         auc = TRUE, quiet = TRUE)
                    
                if (r$auc > max_val) {
                    max_val <- r$auc
                    Xnew <- Xnew2
                    best_thres[i] <- j
                }
            }
        }
        
        message("AUC round ", k, " = ", max_val)
            
        if (identical(best_threshold, best_thres)) {
            break;
        }
        best_threshold <- best_thres
            
    }
    
    threshold_values <- best_threshold
    for (i in 1:ncol(X)) {
        threshold_values[i] <- quantile(X[, i], probs = quantile_cutoff)[best_threshold[i]]
    }
    
    d1 <-  data.frame(y, Xnew)
    form <- ion$regression_form(c("y", colnames(Xnew)), is_quadratic = FALSE)
    
    logistic_model <- glm(form,
                          data = d1,
                          family = "binomial")
    
    
    ret <- list(X = Xnew, thresholds = best_threshold, auc = max_val, 
                threshold_values = threshold_values, logistic_model = logistic_model, variables = colnames(Xnew))
        
    message("AUC = ", max_val)
        
    if (!is.null(pdf_out)) {
        pdf(paste0(pdf_out,".pdf"), width = 6, height = 6)
        v <- ion$logistic_regression_LOOCV(y, ret$X)
        r <- roc(controls = v[y == 1], 
                 cases = v[y != 1], quiet = TRUE)
        plot(r, print.auc = TRUE)
        
        for (i in 1:ncol(ret$X)) {
            str <- colnames(ret$X)[i]
            thres <- quantile(X[, str], probs = quantile_cutoff)
            
            dat <- list(v = X[, str], 
                        groups = y)
            boxplot(v ~ groups, 
                    data = dat,
                    ylab="",
                    main = str,
                    whisklty = 1,
                    staplelty = 0, range = 0)
            
            stripchart(v ~ groups, 
                       vertical = TRUE, 
                       data = dat, 
                       method = "jitter", 
                       add = TRUE, 
                       pch = 20, col = "blue", cex = 2)
            
            abline(h = ret$threshold_values[i], lwd = 2, col = "red")
        }
        dev.off()
    }    
    
    return(ret)
    
}

ion$logistic_regression_stepwise_thresholded_LOO_AUC <- function(y, X, initial_columns = NULL, 
                                                         n = ncol(X) - (if (is.null(initial_columns)) 1 else length(initial_columns)), 
                                                         max_optimization_round = 10, 
                                                         n_quantile = 4,
                                                         pdf_out = "output") {
    
    quantile_cutoff <- seq(0, 1, 1/n_quantile)
    
    if (is.null(initial_columns)) {
        initial_columns <- colnames(X)[ion$find_max_auc(y, X)]    
        message("Starts with maximal AUC: ", initial_columns)
    }
    
    res <- list()
    
#    if (length(initial_columns) == 1) {
#        res[[1]] <- list(X = X[, initial_columns, drop = FALSE], auc = 0)
#        for (i in 1:ncol(res[[1]]$X)) {
#            res[[1]]$X[, i] <- res[[1]]$X[, i] > median(res[[1]]$X[, i]) 
#        }
#    } else {
        res[[1]] <- ion$logistic_regression_thresholded_LOO_AUC(y, X[, initial_columns, drop = FALSE], 
                                                        n_quantile = n_quantile, pdf_out = NULL)
#        #res[[1]] <- list(X = out$X, auc = out$auc)
#    }
    
    # build up stepwise
    for (s in (length(initial_columns)+1):n) {
        
        XX <- res[[s-1]]$X
        
        ff <- setdiff(colnames(X), colnames(XX))
        
        max_val <- -Inf 
        max_f <- -1
        max_thres <- -Inf 
        
        for (i in 1:length(ff)) {
            #cat(i, ":", max_val, "\n")
            thres <- quantile(X[, ff[i]], probs = quantile_cutoff)
            for (j in 2:(length(thres)-1)) {
                Xnew <- cbind(XX, "x" = as.numeric(X[, ff[i]] > thres[j]))
                v <- ion$logistic_regression_LOOCV(y, Xnew)
                
                r <- roc(controls = v[y == 1], 
                         cases = v[y != 1],
                         auc = TRUE, quiet = TRUE)
                
                if (r$auc > max_val) {
                    max_val <- r$auc
                    max_f <- i
                    max_thres <- thres[j]
                }
            }
        }
        
        message("adding ", ff[max_f], ". Optimizing combination ...")
        Xnew <- cbind(XX, "x" = as.numeric(X[, ff[max_f]] > max_thres))
        colnames(Xnew)[ncol(Xnew)] <- ff[max_f]
        
        # Optimize existing combination
        
        best_threshold <- NULL
        
        for (k in 1:max_optimization_round) {
            #cat(k, ":", max_val, "\n")
            
            best_thres <- rep(0, ncol(Xnew))
            
            for (i in 1:ncol(Xnew)) {
                
                thres <- quantile(X[, colnames(Xnew)[i]], probs = quantile_cutoff)
                
                Xnew2 <- Xnew
                max_val <- -Inf
                
                for (j in 2:(length(thres)-1)) {
                    
                    Xnew2[, i] <- as.numeric(X[, colnames(Xnew)[i]] > thres[j])
                    
                    v <- ion$logistic_regression_LOOCV(y, Xnew2)
                    
                    r <- roc(controls = v[y == 1], 
                             cases = v[y != 1],
                             auc = TRUE, quiet = TRUE)
                    
                    if (r$auc > max_val) {
                        max_val <- r$auc
                        Xnew <- Xnew2
                        best_thres[i] <- j
                    }
                }
            }
            
            if (identical(best_threshold, best_thres)) {
                break;
            }
            best_threshold <- best_thres
            
        }
        
        #res[[s]] <- list(X = Xnew, auc = max_val)
        res[[s]] <- ion$logistic_regression_thresholded_LOO_AUC(y, X[, colnames(Xnew), drop = FALSE], 
                                                        n_quantile = n_quantile, pdf_out = NULL)
        
        #message("AUC = ", max_val)
        
        if (!is.null(pdf_out)) {
            pdf(paste0(pdf_out,"-",s,".pdf"), width = 6, height = 6)
            v <- ion$logistic_regression_LOOCV(y, res[[s]]$X)
            r <- roc(controls = v[y == 1], 
                     cases = v[y != 1], quiet = TRUE)
            plot(r, print.auc = TRUE)
            
            for (i in 1:ncol(res[[s]]$X)) {
                str <- colnames(res[[s]]$X)[i]
                thres <- quantile(X[, str], probs = quantile_cutoff)
                
                dat <- list(v = X[, str], 
                            groups = y)
                boxplot(v ~ groups, 
                        data = dat,
                        ylab="",
                        main = str,
                        whisklty = 1,
                        staplelty = 0, range = 0)
                
                stripchart(v ~ groups, 
                           vertical = TRUE, 
                           data = dat, 
                           method = "jitter", 
                           add = TRUE, 
                           pch = 20, col = "blue", cex = 2)
            
                abline(h = thres[best_threshold[i]], lwd = 2, col = "red")
            }
        
            dev.off()
        }
    }
    
    return(res)
}


ion$logistic_regression_thresholded_add <- function(y, X, Xi, quantile_cutoff) {

    d0 <- data.frame(y = y, X)
    m0 <- glm(ion$regression_form(colnames(d0)),
              data = d0,
              family = "binomial")
    thres <- quantile(Xi, probs = quantile_cutoff)
    
    min_p <- 1
    min_i <- floor(length(quantile_cutoff) / 2) + 1 # just take the middle one
    
    for (i in 2:(length(thres)-1)) {
        d1 <- cbind(d0, yy = as.numeric(Xi > thres[i])) 
        m1 <- glm(ion$regression_form(colnames(d1)),
                  data = d1,
                  family = "binomial")
        
        pval <- anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2]
        if (!is.na(pval) && pval < min_p) {
            min_p <- pval
            min_i <- i
        }
    }
    return(list(i = min_i, p = min_p, thres = thres[min_i]))
}


ion$logistic_regression_thresholded_add_0 <- function(y, X, Xi, quantile_cutoff) {
    
    d0 <- data.frame(y = y, X)
    
    thres <- quantile(Xi, probs = quantile_cutoff) 
    
    d1 <- d0
    
    for (i in 2:(length(thres)-1)) {
        d1 <- cbind(d1, as.numeric(Xi > thres[i])) 
        colnames(d1)[ncol(d1)] <- paste0("yy", i)
    }
    
    m0 <- glm(ion$regression_form(colnames(d0), env = environment()),
              data = d1,
              family = "binomial")
    
    scope <- ion$regression_form(colnames(d1), env = environment())
    
    a <- add1(m0, scope = scope, test = "LRT", data = d1)

    return(a)
    
    #return(list(i = min_i, p = min_p, thres = thres[min_i]))
}


ion$logistic_regression_thresholded <- function(y, X, max_optimization_round = 100, 
                                                n_quantile = 4,
                                                pdf_out = "output-str") {
    
    quantile_cutoff <- seq(0, 1, 1 / n_quantile)
    p <- quantile_cutoff[floor(n_quantile/2) + 1]
    
    Xnew <- as.matrix(X)
    for (i in 1:ncol(Xnew)) {
        thres <- quantile(X[, i], probs = p)
        Xnew[, i] <- as.numeric(Xnew[, i] > thres)
    }
    
    best_threshold <- NULL
    
    for (k in 1:max_optimization_round) {
        
        if (ncol(Xnew) > 1) {
        best_thres <- rep(0, ncol(Xnew))
        
        for (i in 1:ncol(Xnew)) {
            
            a <- ion$logistic_regression_thresholded_add(y, Xnew[, -i, drop = FALSE],
                                                         X[, i, drop = FALSE], 
                                                         quantile_cutoff)
            Xnew[, i] <- as.numeric(X[, i, drop = FALSE] > a$thres)
            best_thres[i] <- a$i
        }
        } else {
            best_thres <- floor(length(quantile_cutoff)/2)+1
        }
        
        if (identical(best_threshold, best_thres)) {
            break;
        }
        best_threshold <- best_thres
        
    }
    
    threshold_values <- best_threshold
    for (i in 1:ncol(X)) {
        threshold_values[i] <- quantile(X[, i], probs = quantile_cutoff[best_threshold[i]])
    }
    
    d1 <-  data.frame(y, Xnew)
    form <- ion$regression_form(colnames(d1), is_quadratic = FALSE)
    
    logistic_model <- glm(form,
                          data = d1,
                          family = "binomial")
    
    
    ret <- list(X = Xnew, thresholds = best_threshold,
                threshold_values = threshold_values, 
                logistic_model = logistic_model, variables = colnames(Xnew))
    
    if (!is.null(pdf_out)) {
        require(pROC)
        pdf(paste0(pdf_out,".pdf"), width = 6, height = 6)
        v <- ion$logistic_regression_LOOCV(y, ret$X)
        r <- roc(controls = v[y == 1], 
                 cases = v[y != 1], quiet = TRUE)
        plot(r, print.auc = TRUE)
        
        for (i in 1:ncol(ret$X)) {
            str <- colnames(ret$X)[i]
            thres <- quantile(X[, str], probs = quantile_cutoff)
            
            dat <- list(v = X[, str], 
                        groups = y)
            boxplot(v ~ groups, 
                    data = dat,
                    ylab="",
                    main = str,
                    whisklty = 1,
                    staplelty = 0, range = 0)
            
            stripchart(v ~ groups, 
                       vertical = TRUE, 
                       data = dat, 
                       method = "jitter", 
                       add = TRUE, 
                       pch = 20, col = "blue", cex = 2)
            
            abline(h = ret$threshold_values[i], lwd = 2, col = "red")
        }
        dev.off()
    }    
    
    return(ret)
    
}

ion$logistic_regression_stepwise_thresholded <- function(y, X, initial_columns = NULL, 
                                                           n = ncol(X) - (if (is.null(initial_columns)) 1 else length(initial_columns)), 
                                                           max_optimization_round = 10, 
                                                           n_quantile = 4,
                                                           pdf_out = "output") {
    
    quantile_cutoff <- seq(0, 1, 1/n_quantile)
    
    if (is.null(initial_columns)) {
        initial_columns <- colnames(X)[ion$find_max_auc(y, X)]    
        message("Starts with maximal AUC: ", initial_columns)
    }
    
    res <- list()
    
    res[[1]] <- ion$logistic_regression_thresholded(y, X[, initial_columns, drop = FALSE], 
                                                    n_quantile = n_quantile, pdf_out = NULL)

    # build up stepwise
    for (s in (length(initial_columns)+1):n) {
        
        ff <- setdiff(colnames(X), colnames(res[[s-1]]$X))
        
        min_p <- 2
        min_f <- -1
        min_t <- -1
        
        for (i in 1:length(ff)) {

            a <- ion$logistic_regression_thresholded_add(y, res[[s-1]]$X,
                                                         X[, ff[i], drop = FALSE], 
                                                         quantile_cutoff)
            
            if (a$p < min_p) {
                min_p <- a$p;
                min_f <- i
                min_t <- a$thres
            }
        }
        
        message("adding ", ff[min_f], "; p = ", min_p, "\nOptimizing combination ...")
        Xnew <- cbind(res[[s-1]]$X, "yy" = as.numeric(X[, ff[min_f]] > min_t))
        colnames(Xnew)[ncol(Xnew)] <- ff[min_f]
        
        res[[s]] <- ion$logistic_regression_thresholded(y, X[, colnames(Xnew), drop = FALSE], 
                                                        max_optimization_round = max_optimization_round,
                                                        n_quantile = n_quantile, pdf_out = NULL)
        
        if (!is.null(pdf_out)) {
            pdf(paste0(pdf_out, "-", s, ".pdf"), width = 6, height = 6)
            v <- ion$logistic_regression_LOOCV(y, res[[s]]$X)
            r <- roc(controls = v[y == 1], 
                     cases = v[y != 1], quiet = TRUE)
            plot(r, print.auc = TRUE)
            
            for (i in 1:ncol(res[[s]]$X)) {
                str <- colnames(res[[s]]$X)[i]
                thres <- quantile(X[, str], probs = quantile_cutoff)
                
                dat <- list(v = X[, str], 
                            groups = y)
                boxplot(v ~ groups, 
                        data = dat,
                        ylab="",
                        main = str,
                        whisklty = 1,
                        staplelty = 0, range = 0)
                
                stripchart(v ~ groups, 
                           vertical = TRUE, 
                           data = dat, 
                           method = "jitter", 
                           add = TRUE, 
                           pch = 20, col = "blue", cex = 2)
                
                abline(h = res[[s]]$threshold_values[i], lwd = 2, col = "red")
            }
            
            dev.off()
        }
    }
    
    return(res)
}


ion$logistic_regression_stepwise_thresholded_LOOCV <- function(y, X, initial_columns = NULL, 
                                                               n = ncol(X) - (if (is.null(initial_columns)) 1 else length(initial_columns)), 
                                                               n_quantile = seq(2, 20, 2),
                                                               max_optimization_round = 10) {
    if (nrow(X) != length(y)) {
        stop("The number of samples (X, y) does not match.\n")
    }
    
    v <- matrix(0, nrow = length(n_quantile), ncol = n)
    
    for (nq in 1:length(n_quantile)) {    

        z <- matrix(0, nrow = nrow(X), ncol = n)
        
        for (r in 1:nrow(X)) {
            message("n_quantile = ", nq, "/", length(n_quantile), "; leaving out ", r, "/", nrow(X))
            out <- ion$logistic_regression_stepwise_thresholded(y[-r], X[-r, , drop = FALSE], 
                                                                n = n, n_quantile = n_quantile[nq], 
                                                                max_optimization_round = max_optimization_round,
                                                                pdf_out = NULL)
            for (i in 1:n) {
                z[r, i] <- ion$logistic_regression_thresholded_predict(out[[i]], X[r, out[[i]]$variables, drop = FALSE])
            }
        }
        
        for (i in 1:n) {
            v[nq, i] <- ion$auc(y, z[, i], make_plot = FALSE)
        }
    }
    
    return(v)
    
}

# Others ----

# This code contains a modified version of the heatmap.2 function in the R gplots
# package (version 3.0.1).

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
    
    #TP 220430
    #if (is.null(Rowv) || is.na(Rowv))
    if (is.null(Rowv))
        Rowv <- FALSE
    #if (is.null(Colv) || is.na(Colv))
    if (is.null(Colv))
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

    csc_Thang <- matrix(0, nrow = dim(csc)[1], ncol = dim(csc)[2])

    csc.colors = matrix()
    csc.names = names(table(csc))
    csc.i = 1
    for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc_Thang[csc == csc.name] = csc.i
        csc.i = csc.i + 1
    }
    #csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
    csc <- csc_Thang

a <- as.vector(csc.colors)
a[nchar(a) <= 1] <- "white"
a[startsWith(a, "`")] <- "white"

image(csc, col = a, axes = F)

a <- as.vector(csc.colors)

for (i in 1:nrow(csc)) {
    for (j in 1:ncol(csc)) {
        if (nchar(a[csc[i,j]]) <= 1) {
            if (ncol(csc) > 1) {
                text((i-1)/(nrow(csc)-1), (j-1)/(ncol(csc)-1), a[csc[i,j]], cex = 1.5)
            }
            else {
                text((i-1)/(nrow(csc)-1), 0, a[csc[i,j]], cex = 1.5)
            }
        }
        if (startsWith(a[csc[i,j]],"`")) {
            if (ncol(csc) > 1) {
                eval(parse(text=paste0("points(",
                                       (i-1)/(nrow(csc)-1),
                                       ",",
                                       (j-1)/(ncol(csc)-1),
                                       ",",
                                       substring(a[csc[i,j]],2,1000),
                                       ")")))
            } else {
                eval(parse(text=paste0("points(",
                                       (i-1)/(nrow(csc)-1),
                                       ",",
                                       0,
                                       ",",
                                       substring(a[csc[i,j]],2,1000),
                                       ")")))
            }
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
            #tmp <- dev.size("in")
            #key_margins <- c(tmp[1]/2, 1, 0.5, 2)

            #if (dev.size("in")[2] < 2) {
            #    message("Why such a short figure?")
            #}
            #key_margins <- c(dev.size("in")[2] * 1.5 , 1, 0.5, 2) # 6/4
            key_margins <- c(dev.size("in")[2]+0.5, 1, 0.5, 2)
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
            message("\nCannot make legend key.\nRun dev.off(), then make the margin bigger, or color_key_margins smaller, or set key = FALSE.\n")
            cat("Current color_key_margins", key_margins, "\n")
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
