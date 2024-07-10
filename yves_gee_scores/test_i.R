 library(lavaan)
 library(pbivnorm)
Data <- read.csv("test_data.csv")

model <- '
    f1 =~ u1 + u2 + u3
    f2 =~ u4 + u5 + u6
    f1 + f2 ~ start(0.1)*x1 + start(-0.1)*x2
'
fit <- sem(model, data= Data, estimator = "WLS",
           ordered = c("u1", "u2", "u3", "u4", "u5", "u6"))
# summary(fit)

# PML
# fit <- cfa(model, data=Data, estimator="PML", verbose=TRUE)
# summary(fit, fit.measures=TRUE, standardized=TRUE)


# GEE based on Reboussin & Liang (1998)
# Reboussin, B. A., & Liang, K. Y. (1998). An estimating equations approach for
# the LISCOMP model. Psychometrika, 63, 165-182.

# product unfitted model
fit0 <- sem(model, data = Data, do.fit = FALSE, estimator = "WLS",
            ordered = c("u1", "u2", "u3", "u4", "u5", "u6"))
start.x <- coef(fit0)

get_y_yy_exp.i <- function(x, lavmodel = NULL, x.i) {
  # insert parameters in model matrices
  this.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = x)
  Sigma.star <- lavaan:::computeSigmaHat(lavmodel = this.model)[[1]]
  MLIST <- this.model@GLIST # only 1 group
  LAMBDA <- MLIST$lambda
  BETA   <- MLIST$beta
  GAMMA  <- MLIST$gamma
  ALPHA  <- MLIST$alpha
  TAU    <- MLIST$tau
  IB.inv <- lavaan:::.internal_get_IB.inv(MLIST = MLIST)

  mu.i <- drop(LAMBDA %*% IB.inv %*% ALPHA +
               LAMBDA %*% IB.inv %*% GAMMA %*% x.i)
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
  mu.star.i <- lavaan:::computeEYx(lavmodel = this.model,
    lavsamplestats = fit0@SampleStats, eXo = fit0@Data@eXo)[[1]]

  # lower integration points:
  # tmp <- (matrix(this.model@GLIST$tau,
  #                nrow = nrow(mu.star.i), ncol = ncol(mu.star.i), byrow = TRUE)
  #          - mu.star.i)
  # diag.Sigma.star <- diag(Sigma.star) # all ones here
  # tmp <- t( t(tmp) / diag.Sigma.star ) #not needed here
  # Ey_xi <- pnorm(tmp, lower.tail = FALSE)

  # using upper integration point
  # = int_{-inf}^{-tau+u_ij} f(z) dz
  upper <- t( t(mu.star.i) - drop(this.model@GLIST$tau) )
  y_exp <- pnorm(upper)

  # E(Y_j E_k|x) = 1 * P(Y_j=1, Y_k=1|x) + 0 
  #              = P(yj_star > tau_j, yk_star > tau_k | x)
  #              = int_{tau_j}^{inf} 
  #                int_{tau_k}^{inf} f(yj,yk;x,rho_jk) dyj dyk
  #              = int_{tau_j-mu_ij}^{inf} 
  #                int_{tau_k-mu_ik}^{inf} g(zj,zk;x,rho_jk) dzj dzk
  nvar <- this.model@nvar[1]
  pstar <- nvar * (nvar - 1) / 2
  w_obs <- w_exp <- Eyy_xi <- matrix(0, nrow(mu.star.i), pstar)
  PSTAR <- matrix(0, nvar, nvar) # utility matrix, to get indices
  PSTAR[lav_matrix_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar

  y_obs <- fit@Data@X[[1]] - 1
  y_obs.c <- y_obs - y_exp # model-based centering?

  # changing sign, to get 'upper' integration limits
  for (k in seq_len(nvar - 1L)) {
    for (j in (k + 1L):nvar) {
      pstar.idx <- PSTAR[j, k]
      Eyy_xi[, pstar.idx] <- pbivnorm(upper[,j], upper[,k],
                                      rho = Sigma.star[j, k])
      w_exp[, pstar.idx] <- Eyy_xi[, pstar.idx] - (y_exp[,j] * y_exp[,k])
      w_obs[, pstar.idx] <- y_obs.c[,j] * y_obs.c[,k]
    }
  }

  # error
  ey_error <- y_obs - y_exp
  ew_error <- w_obs - w_exp
  both_error <- cbind(ey_error, ew_error)

  SCORES <- matrix(0, nrow = nrow(y_obs), ncol = length(x))
  for(i in 1:nrow(y_obs)) {
    cat("i = ", i, "\n")

    # covy
    covy <- lav_matrix_vech_reverse(w_exp[i,], diagonal = FALSE)
    diag(covy) <- y_exp[i,] * (1 - y_exp[i,])

    # covw
    covw <- diag(pstar)
    #diag(covw) <- w_obs[i,]^2 - w_exp[i,]^2 # correct?
	diag(covw) <- colMeans(w_obs^2) - w_exp[i,]^2

    # B_i
    B_i <- lav_matrix_bdiag(covy, covw)
    B_i.inv <- solve(B_i)

    # Delta
    Delta.i <- numDeriv:::jacobian(func = get_y_yy_exp.i, x = x,
                                   lavmodel = fit0@Model, 
								   x.i = fit0@Data@eXo[[1]][i,])

    SCORES[i, ] <- drop( t(Delta.i) %*% B_i.inv %*% as.matrix(both_error[i,]) )
  }

  SCORES
}

# SC <- get_scores( coef(fit) )
# colSums(SC) 
# these values should be close to zero, 
# although GEE estimates are not WLS estimates
