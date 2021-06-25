library(lmPerm) # package to perform permutational significance tests

# phylogenetic distances
treeD <- read.csv("data/rescaled_tree_insecta6.csv",
                  header = TRUE, row.names = 1)

# synteny distances
syntD <- read.table("data/Insecta_matrix_matched_to_phylo_mod3.txt",
                   header = TRUE, row.names = 1)

# map of location of orders in distance matrices
taxaii <- data.frame(ord = c('Hemiptera', 'Hymenoptera', 'Coleoptera',
                             'Diptera', 'Lepidoptera'),
                     i0 = c(2, 16, 42, 49, 65),
                     i1 = c(14, 40, 48, 64, 143))

# function to etract lower triangles for a given order
getLowerTri <- function(ord) {
    o <- which(taxaii$ord == ord)
    ii <- taxaii$i0[o]:taxaii$i1[o]

    mt <- treeD[ii, ii]
    ms <- syntD[ii, ii]

    return(data.frame(ord = ord,
                      treeD = mt[lower.tri(mt)],
                      syntD = ms[lower.tri(ms)]))
}


# make data.frame with all data on distances and grouops
ordDat <- lapply(taxaii$ord, getLowerTri)
ordDat <- do.call(rbind, ordDat)


# make three different models: linear, exponential, power law

# linear
yLin <- ordDat$syntD
xLin <- ordDat$treeD
ord <- ordDat$ord
modLin <- lmp(yLin ~ xLin * ord, perm = 'Prob')
summary(modLin)[['r.squared']]

# exponential
yExp <- log(ordDat$syntD)
xExp <- ordDat$treeD
modExp <- lmp(yExp ~ xExp * ord, perm = 'Prob', center = FALSE)
summary(modExp)[['r.squared']]

foo <- lm(yExp ~ xExp * ord)

cbind(modExp$coefficients, foo$coefficients)


# power law
yPow <- log(ordDat$syntD)
xPow <- log(ordDat$treeD)
modPow <- lmp(yPow ~ xPow * ordDat$ord, perm = 'Prob')
summary(modPow)[['r.squared']]


predict(modExp)

modExp$coefficients
plot(ordDat$treeD, ordDat$syntD, col = factor(ordDat$ord))
legend('topright', legend = levels(factor(ordDat$ord)),
       col = 1:5, pch = 1)
boo <- data.frame(xExp = seq(0, 1, length.out = 100), ord = 'Lepidoptera')
lines(boo$xExp, exp(predict(foo, boo)))


modExpCurve <- function(i) {
    b <- modExp$coefficients
    b <- c(b, 'ord0' = 0, 'xExp:ord0' = 0)

    curve(exp(b[1] + b[paste0('ord', i - 1)]) *
              exp(x * (b['xExp'] + b[paste0('xExp:ord', i - 1)])),
          col = i, add = TRUE)
}

modExpCurve(5)


