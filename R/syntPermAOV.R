#' Permutational ANOVA for specific analysis of relationship between
#' synteny and phylogenetic distance across different order
#'
#' @param treeD vector of phylogenetic distances
#' @param syntD vector of synteny distances
#' @param ord vector of order names
#' @param B number of permutations to preform
#' @return returns a data.frame similar to a standard ANOVA table

syntPermAOV <- function(treeD, syntD, ord, B = 999) {
    dat <- data.frame(y = syntD, x = treeD, g = factor(ord))
    nord <- nlevels(dat$g)

    df <- c(1, rep(nord - 1, 2))
    df <- c(df, nrow(dat) - sum(df) - 1)

    # make models of reducing complexity from full (X3) to slope only (X1)
    X3 <- model.matrix(y ~ x * g, data = dat)
    X2 <- X3[, 1:(2 + nord - 1)]
    X1 <- X3[, 1:2]

    # do matrix calculations on different models
    pHat1 <- partialHat(X1)
    pHat2 <- partialHat(X2)
    pHat3 <- partialHat(X3)

    # calculate F statistics
    fobs <- fval(list(pHat1, pHat2, pHat3),
                 list(X1, X2, X3),
                 df, dat$y)

    # calculate null F distributions
    fnull <- lapply(1:B, function(i) {
        newy <- sample(dat$y)
        f <- fval(list(pHat1, pHat2, pHat3),
                  list(X1, X2, X3),
                  df, newy)
        return(f)
    })

    fnull <- do.call(rbind, c(fnull, fobs))

    # calculate pvals from fnull
    pval <- numeric(length(fobs))
    for(i in 1:length(pval)) {
        pval[i] <- mean(fnull[, i] >= fobs[i])
    }

    tab <- data.frame(DF = paste(c(df[1:3], sum(df[1:3])), df[4], sep = ', '),
                      FVal = fobs, Pval = pval)

    return(tab)
}


#' preforms part of the matrix calculation needed to solve for the parameters of
#' a linear regression; the full solution is
#' `bHat <- solve(t(X) %*% X) %*% t(X) %*% y`. We leave out the part about `y`
#' because we are going to permute `y`, but all calculations not involving `y`
#' can be done once to save computation
#' @param X the model matrix
#' @note only intended for use inside `syntPermAOV`

partialHat <- function(X) {
    return(solve(t(X) %*% X) %*% t(X))
}


#' calculate predicted y using the output of `partialHat`
#' @param pHat the output of `partialHat`
#' @param X the model matrix
#' @param y the observed (possibly permuted) response variable
#' @note only intended for use inside `syntPermAOV`

yhat <- function(pHat, X, y) {
    b <- pHat %*% y
    return(X %*% b)
}


#' function to calculate F values
#' @param pHat the output of `partialHat`
#' @param X the model matrix
#' @param df vector of degrees of freedom
#' @param y the observed (possibly permuted) response variable
#' @note only intended for use inside `syntPermAOV`

fval <- function(pHat, X, df, y) {
    y1 <- yhat(pHat[[1]], X[[1]], y)
    y2 <- yhat(pHat[[2]], X[[2]], y)
    y3 <- yhat(pHat[[3]], X[[3]], y)

    ssq <- c(sum((y1 - mean(y))^2),
             sum((y2 - y1)^2),
             sum((y3 - y2)^2),
             sum((y3 - y)^2))

    msq <- ssq / df

    msqTot <- sum(ssq[1:3]) / sum(df[1:3])

    f <- c(msq[1:3], msqTot) / msq[4]

    return(f)
}

