# function copied from `lmPerm::anova.lmp` but modified to silently return the
# ANOVA table instead of returning `NULL`

anova.lmp <- function (object, ...) {
    if (length(list(object, ...)) > 1)
        stop("The first argument cannot be an array or list.")
    mo <- NCOL(object$resid)
    if (mo > 1)
        respNames <- colnames(object$residuals)
    else respNames <- deparse(formula(object)[[2]])
    for (id in 1:mo) {
        resid <- as.matrix(object$resid)[, id]
        effects <- as.matrix(object$effects)[, id]
        w <- object$weights
        ssr <- sum(if (is.null(w)) resid^2 else w * resid^2)
        dfr <- df.residual(object)
        p <- object$rank
        if (p > 0) {
            p1 <- 1:p
            comp <- effects[p1]
            asgn <- object$assign[object$qr$pivot][p1]
            nmeffects <- c("(Intercept)", attr(object$terms,
                                               "term.labels"))
            tlabels <- nmeffects[1 + unique(asgn)]
            ss <- c(unlist(lapply(split(comp^2, asgn), sum)),
                    ssr)
            df <- c(unlist(lapply(split(asgn, asgn), length)),
                    dfr)
        }
        else {
            ss <- ssr
            df <- dfr
            tlabels <- character(0)
        }
        ms <- ss/df
        if (is.null(object$perm))
            doPerm <- FALSE
        else {
            perm = object$perm$perm
            doPerm <- (perm == "Prob" || perm == "SPR" || perm ==
                           "Exact")
        }
        if (doPerm) {
            table <- data.frame(df, ss, ms)
            const <- attr(object$terms, "intercept")
            P <- c(as.vector(object$perm$P[, id]), NA)
            perm <- object$perm$perm
            if (perm == "Prob") {
                table$Mn <- c(as.vector(object$perm$Mn[, id]),
                              NA)
                table$P <- P
                dimnames(table) <- list(c(tlabels, "Residuals"),
                                        c("Df", "R Sum Sq", "R Mean Sq", "Iter",
                                          "Pr(Prob)"))
            }
            else if (perm == "SPR") {
                table$Mn <- c(as.vector(object$perm$Mn[, id]),
                              NA)
                table$P <- P
                table$accept <- c(as.vector(object$perm$accept),
                                  NA)
                dimnames(table) <- list(c(tlabels, "Residuals"),
                                        c("Df", "R Sum Sq", "R Mean Sq", "Iter",
                                          "Pr(SPR)", "Accept"))
            }
            else {
                table$P <- P
                dimnames(table) <- list(c(tlabels, "Residuals"),
                                        c("Df", "R Sum Sq", "R Mean Sq",
                                          "Pr(Exact)"))
            }
            if (const)
                table <- table[-1, , drop = FALSE]
            print(structure(table, heading = c("Analysis of Variance Table\n",
                                               paste("Response:",
                                                     respNames[id])),
                            class = c("anova", "data.frame")))
        }
        else {
            f <- ms/(ssr/dfr)
            P <- pf(f, df, dfr, lower.tail = FALSE)
            table <- data.frame(df, ss, ms, f, P)
            table[length(P), 4:5] <- NA
            dimnames(table) <- list(c(tlabels, "Residuals"),
                                    c("Df", "R Sum Sq", "R Mean Sq", "F value",
                                      "Pr(>F)"))
            if (attr(object$terms, "intercept"))
                table <- table[-1, ]
            print(structure(table, heading = c("Analysis of Variance Table\n",
                                               paste("Response:", respNames[id],
                                                     object$info)),
                            class = c("anova", "data.frame")))
        }
    }

    invisible(table)
}
