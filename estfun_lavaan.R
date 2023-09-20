
setwd("C:/Users/classe/Desktop/Diss/Paper3/R")
source("support.R")




# shortcuts
lavdata        <- object@Data
lavmodel       <- object@Model
lavsamplestats <- object@SampleStats
lavoptions     <- object@Options


ntab <- unlist(lavdata@norig)
ntot <- sum(ntab)
npar <- lav_object_inspect_npar(object)
nvar <- ncol(lavsamplestats@cov[[1]])


moments <- lavaan::fitted(object)
N1 <- 1
X <- lavdata@X[[1]]

Score.mat <- matrix(NA, ntot, npar) #empty matrix


J <- matrix(1, 1L, ntab[1]) ## FIXME: needed? better maybe rowSums/colSums?
J2 <- matrix(1, nvar, nvar)
diag(J2) <- 0.5
Mu.hat <- moments$mean
Sigma.hat <- moments$cov
Sigma.inv <- inv.chol(Sigma.hat, logdet=FALSE)
group.w <- (unlist(lavsamplestats@nobs)/lavsamplestats@ntotal)

################################################################################
################################## ML ##########################################
################################################################################

# scores.H1 (H1 = saturated model)
#Equation in Appendix B 
mean.diff <- t(t(X) - Mu.hat %*% J) 
dx.Mu <- -1 * mean.diff %*% Sigma.inv 
dx.Sigma <- t(matrix(apply(mean.diff, 1L, 
                           function(x) lav_matrix_vech(- J2 * (Sigma.inv %*% (tcrossprod(x)*N1 - Sigma.hat) %*% Sigma.inv))), ncol=nrow(mean.diff))) #fuer WLS anders!

scores.H1 <- cbind(dx.Mu, dx.Sigma)          ########### !!! ncol =  (vectorized non-redundant elements of pxp matrix) + (n means)

Delta <- computeDelta(lavmodel = lavmodel)[[1]] #should also work for WLS...

wi <- lavdata@case.idx[[1]]
Score.mat[wi,] <- -scores.H1 %*% Delta


Score.mat[wi,] <- (-1/ntot) * Score.mat[wi,] #scaling


# provide column names
colnames(Score.mat) <- names(lav_object_inspect_coef(object,
                                                     type = "free", add.labels = TRUE))
