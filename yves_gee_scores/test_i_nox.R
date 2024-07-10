 library(lavaan)
 library(pbivnorm)
Data <- read.csv("yves_gee_scores/test_data.csv")

model <- '
    f1 =~ u1 + u2 + u3
    f2 =~ u4 + u5 + u6
'
fit.wls <- sem(model, data= Data, estimator = "WLS",
           ordered = c("u1", "u2", "u3", "u4", "u5", "u6"))
coef(fit.wls)

fit.dwls <- sem(model, data= Data, estimator = "DWLS",
                ordered = c("u1", "u2", "u3", "u4", "u5", "u6"))
coef(fit.dwls)

# PML
fit.pml <- cfa(model, data=Data, estimator="PML", verbose=TRUE, ordered = TRUE)
coef(fit.pml)


# GEE based on Reboussin & Liang (1998)
# Reboussin, B. A., & Liang, K. Y. (1998). An estimating equations approach for
# the LISCOMP model. Psychometrika, 63, 165-182.

# product unfitted model
fit0 <- sem(model, data = Data, do.fit = FALSE, estimator = "WLS",
            ordered = c("u1", "u2", "u3", "u4", "u5", "u6"))
start.x <- coef(fit0)

get_y_yy_exp_nox.i <- function(x, lavmodel = NULL) {
  # insert parameters in model matrices
  this.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = x)
  Sigma.star <- lavaan:::computeSigmaHat(lavmodel = this.model)[[1]]
  MLIST <- this.model@GLIST # only 1 group
  LAMBDA <- MLIST$lambda
  ALPHA  <- MLIST$alpha
  TAU    <- MLIST$tau
  IB.inv <- lavaan:::.internal_get_IB.inv(MLIST = MLIST)

  mu.i <- drop(LAMBDA %*% IB.inv %*% ALPHA)
  upper.i <- mu.i - drop(TAU)
  y_exp.i <- pnorm(upper.i)

  nvar <- length(y_exp.i)
  pstar <- nvar * (nvar - 1 ) / 2

  Eyy_xi.i <- w_exp.i <- numeric(pstar)
  pstar.idx <- 0L
  for (k in seq_len(nvar - 1L)) {
    for (j in (k + 1L):nvar) {
      pstar.idx <- pstar.idx + 1L
      Eyy_xi.i[pstar.idx] <- pbivnorm(upper.i[j], upper.i[k], rho = Sigma.star[j, k])
      w_exp.i[pstar.idx] <- Eyy_xi.i[pstar.idx] - (y_exp.i[j] * y_exp.i[k])
    }
  }
  c(y_exp.i, w_exp.i)
}


get_scores <- function(x) {
  # insert parameters in model matrices
  this.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = x)
  Sigma.star <- lavaan:::computeSigmaHat(lavmodel = this.model)[[1]]

  # E(Y|x) = 1 * P(Y=1|x) + 0 * P(Y=0|x)
  #        = P(Y=1|x)
  #        = P(y_star > tau | x)
  #        = int_{tau}^{Inf} f(y; u_ij, sigma_jj) dy
  #        = int_{tau-u_ij}^{Inf} f(z) dz
  mu.star.i <- lavaan:::computeEY(lavmodel = this.model,
    lavsamplestats = fit0@SampleStats)[[1]] # all zeros

  # lower integration points:
  # tmp <- (matrix(this.model@GLIST$tau,
  #                nrow = nrow(mu.star.i), ncol = ncol(mu.star.i), byrow = TRUE)
  #          - mu.star.i)
  # diag.Sigma.star <- diag(Sigma.star) # all ones here
  # tmp <- t( t(tmp) / diag.Sigma.star ) #not needed here
  # Ey_xi <- pnorm(tmp, lower.tail = FALSE)

  # using upper integration point
  # = int_{-inf}^{-tau+u_ij} f(z) dz
  upper <- mu.star.i - drop(this.model@GLIST$tau)
  N <- nrow(fit0@Data@X[[1]])
  nvar <- ncol(fit0@Data@X[[1]])
  y_exp <- matrix(pnorm(upper), nrow = N, ncol = nvar, byrow = TRUE)
  # note: each row of y_exp is the same as colMeans(fit0@Data@X[[1]] - 1)

  # E(Y_j E_k|x) = 1 * P(Y_j=1, Y_k=1|x) + 0 
  #              = P(yj_star > tau_j, yk_star > tau_k | x)
  #              = int_{tau_j}^{inf} 
  #                int_{tau_k}^{inf} f(yj,yk;x,rho_jk) dyj dyk
  #              = int_{tau_j-mu_ij}^{inf} 
  #                int_{tau_k-mu_ik}^{inf} g(zj,zk;x,rho_jk) dzj dzk
  nvar <- this.model@nvar[1]
  pstar <- nvar * (nvar - 1) / 2
  w_obs <- w_exp <- Eyy_xi <- matrix(0, N, pstar)
  PSTAR <- matrix(0, nvar, nvar) # utility matrix, to get indices
  PSTAR[lav_matrix_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar

  y_obs <- fit0@Data@X[[1]] - 1
  y_obs.c <- y_obs - y_exp # model-based centering?

  # changing sign, to get 'upper' integration limits
  for (k in seq_len(nvar - 1L)) {
    for (j in (k + 1L):nvar) {
      pstar.idx <- PSTAR[j, k]
      Eyy_xi[, pstar.idx] <- pbivnorm(upper[j], upper[k],
                                      rho = Sigma.star[j, k])
      w_exp[, pstar.idx] <- Eyy_xi[, pstar.idx] - (y_exp[,j] * y_exp[,k])
      w_obs[, pstar.idx] <- y_obs.c[,j] * y_obs.c[,k]
    }
  }

  # error
  ey_error <- y_obs - y_exp
  ew_error <- w_obs - w_exp
  both_error <- cbind(ey_error, ew_error)

  # covy
  covy <- lav_matrix_vech_reverse(w_exp[1,], diagonal = FALSE)
  diag(covy) <- y_exp[1,] * (1 - y_exp[1,])

  # covw
  covw <- diag(pstar)
  #diag(covw) <- w_obs[i,]^2 - w_exp[i,]^2 # correct?
  diag(covw) <- colMeans(w_obs^2) - w_exp[1,]^2

  # B_i
  B_i <- lav_matrix_bdiag(covy, covw)
  B_i.inv <- solve(B_i)

  # Delta
  Delta.i <- numDeriv:::jacobian(func = get_y_yy_exp_nox.i, x = x,
                                 lavmodel = fit0@Model)

  SCORES <- both_error %*% B_i.inv %*% Delta.i

  SCORES
}

# estimate parameters using GEE
max_iter <- 200L
old.x <- start.x
# old.x[5:10] <- (1:6)/10 # just to see what happens
old.x[5:10] <- 0
for(iter in seq_len(max_iter)) {

  # verbose
  cat("iter = ", iter, "\n")
  print(old.x)
  cat("\n")

  # insert parameters in model matrices
  this.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = old.x)
  Sigma.star <- lavaan:::computeSigmaHat(lavmodel = this.model)[[1]]

  MLIST <- this.model@GLIST # only 1 group 
  LAMBDA <- MLIST$lambda
  ALPHA  <- MLIST$alpha
  TAU    <- MLIST$tau
  IB.inv <- lavaan:::.internal_get_IB.inv(MLIST = MLIST)
  
  mu.star.i <- drop(LAMBDA %*% IB.inv %*% ALPHA)

  upper <- mu.star.i - drop(TAU)
  N <- nrow(fit0@Data@X[[1]])
  nvar <- ncol(fit0@Data@X[[1]])
  y_exp <- matrix(pnorm(upper), nrow = N, ncol = nvar, byrow = TRUE)

  nvar <- this.model@nvar[1]
  pstar <- nvar * (nvar - 1) / 2
  w_obs <- w_exp <- Eyy_xi <- matrix(0, N, pstar)

  y_obs <- fit0@Data@X[[1]] - 1
  y_obs.c <- y_obs - y_exp # model-based centering?

  # changing sign, to get 'upper' integration limits
  pstar.idx <- 0L
  for (k in seq_len(nvar - 1L)) {
    for (j in (k + 1L):nvar) {
	  pstar.idx <- pstar.idx + 1L
      Eyy_xi[, pstar.idx] <- pbivnorm(upper[j], upper[k],
                                      rho = Sigma.star[j, k])
      w_exp[, pstar.idx] <- Eyy_xi[, pstar.idx] - (y_exp[,j] * y_exp[,k])
      w_obs[, pstar.idx] <- y_obs.c[,j] * y_obs.c[,k]
    }
  }

  # error
  ey_error <- y_obs - y_exp
  ew_error <- w_obs - w_exp
  both_error <- cbind(ey_error, ew_error)

  # covy 
  covy <- lav_matrix_vech_reverse(w_exp[1,], diagonal = FALSE)
  diag(covy) <- y_exp[1,] * (1 - y_exp[1,])
  
  # covw
  covw <- diag(pstar)
  #diag(covw) <- w_obs[i,]^2 - w_exp[i,]^2 # correct?
  diag(covw) <- colMeans(w_obs^2) - w_exp[1,]^2

  # B_i
  B_i <- lav_matrix_bdiag(covy, covw)
  B_i.inv <- solve(B_i)

  # Delta
  Delta.i <- numDeriv:::jacobian(func = get_y_yy_exp_nox.i, x = old.x,
                                 lavmodel = fit0@Model)

  # Scores
  SCORES <- both_error %*% B_i.inv %*% Delta.i
  sum.SCORES <- colSums(SCORES)

  DBiD <- N * (t(Delta.i) %*% B_i.inv %*% Delta.i)

  # update parameters
  step <- 0.2
  for(nstep in 1:200L) {
    if(nstep > 1L) {
      cat("\n")
      cat("step halving: stepsize = ", step / 2, "\n")
    }
    step <- step / 2
    new.x <- old.x + step * drop( solve(DBiD) %*% sum.SCORES )
    # check
    tmp.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = new.x)
    Sigma.star <- lavaan:::computeSigmaHat(lavmodel = tmp.model)[[1]]
    if(all(abs(lav_matrix_vech(Sigma.star, diagonal = FALSE)) < 0.99)) {
      # good, proceed
      break
    }
  }
  if(nstep == 200L) {
    cat("step halving failed; bailing out...\n")
    break
  }

  # check for convergence
  diff <- abs(new.x - old.x)
  max.diff <- max(diff)
  if(max.diff < 1e-08) {
    cat("\n")
    cat("CONVERGED!\n")
    cat("final estimates:\n")
    print(new.x)
    break
  } else {
    cat("max diff = ", max.diff, "\n")
  }

  old.x <- new.x
}

# interestingly: GEE == PML == DWLS in this case:
coef(fit.pml) - new.x
coef(fit.dwls) - new.x

# compute GEE scores
SC <- get_scores( new.x )
colnames(SC) <- names(coef(fit0))

# compute WLS based scores
fit1 <- sem(model, data = Data, start = new.x, 
            do.fit = FALSE, estimator = "WLS",
            ordered = c("u1", "u2", "u3", "u4", "u5", "u6"))

SC2 <- lavScores(fit1)

# Franz Classe version of GEE scores
source("estfun_WLS_complex.R")
SC3 <- estfun.WLS(fit1)

# round(diag(cor(SC, SC2)), 3)
# f1=~u2 f1=~u3 f2=~u5 f2=~u6  u1|t1  u2|t1  u3|t1  u4|t1  u5|t1  u6|t1 f1~~f1
#  0.869  0.887  0.875  0.868 -0.993 -0.989 -0.993 -0.992 -0.994 -0.994  0.957
# f2~~f2 f1~~f2
#  0.957  0.991

# round(diag(cor(SC, SC3)), 3)
# f1=~u2 f1=~u3 f2=~u5 f2=~u6  u1|t1  u2|t1  u3|t1  u4|t1  u5|t1  u6|t1 f1~~f1
#  1.000  1.000  1.000  1.000 -0.726 -0.681 -0.723 -0.760 -0.706 -0.745  1.000
# f2~~f2 f1~~f2
#  1.000  1.000

