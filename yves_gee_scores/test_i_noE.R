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
summary(fit, fit.measures=TRUE, standardized=TRUE)

# PML
# fit <- cfa(model, data=Data, estimator="PML", verbose=TRUE)
# summary(fit, fit.measures=TRUE, standardized=TRUE)


# GEE based on Reboussin & Liang (1998)
# Reboussin, B. A., & Liang, K. Y. (1998). An estimating equations approach for
# the LISCOMP model. Psychometrika, 63, 165-182.

# product unfitted model
fit0 <- sem(model, data = Data, do.fit = FALSE,
            ordered = c("u1", "u2", "u3", "u4", "u5", "u6"))
start.x <- coef(fit0)


# helper functions (mostly for Delta.i)
get_w_obs.i <- function(y) {
  nvar <- length(y)
  pstar <- nvar * (nvar - 1 ) / 2
  w_obs.i <- numeric(pstar)
  pstar.idx <- 0L
  for (k in seq_len(nvar - 1L)) {
    for (j in (k + 1L):nvar) {
      pstar.idx <- pstar.idx + 1L
      w_obs.i[pstar.idx] <- y[j] * y[k]
    }
  }
  w_obs.i
}

get_mu_star.i <- function(x, lavmodel = NULL, x.i) {
  # insert parameters in model matrices
  this.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = x)
  MLIST <- this.model@GLIST # only 1 group
  LAMBDA <- MLIST$lambda
  BETA   <- MLIST$beta
  GAMMA  <- MLIST$gamma
  ALPHA  <- MLIST$alpha
  IB.inv <- lavaan:::.internal_get_IB.inv(MLIST = MLIST)
  
  mu.i <- drop(LAMBDA %*% IB.inv %*% ALPHA + 
               LAMBDA %*% IB.inv %*% GAMMA %*% x.i)
  mu.i
}

get_upper.i <- function(x, lavmodel = NULL, x.i) {
  # insert parameters in model matrices
  this.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = x)
  MLIST <- this.model@GLIST # only 1 group
  TAU    <- MLIST$tau
  mu.i <- get_mu_star.i(x = x, lavmodel = lavmodel, x.i = x.i)
  upper.i <- mu.i - drop(TAU)  
  upper.i
}

get_y_exp.i <- function(x, lavmodel = NULL, x.i) {
  # insert parameters in model matrices
  this.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = x)
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
  y_exp.i
}

grad_y_exp.i <- function(x, lavmodel = NULL, x.i) {
  this.model <- lav_model_set_parameters(lavmodel = fit0@Model, x = x)
  MLIST <- this.model@GLIST # only 1 group
  TAU    <- MLIST$tau
  tmp1 <- dmu_x.dtheta(x = start.x, lavmodel = lavmodel, x.i = x.i)[[1]]
  tau.idx <- which(names(lavmodel@GLIST) == "tau")
  tmp1[cbind(1:nrow(tmp), lavmodel@x.free.idx[[tau.idx]])] <- -1
  upper.i <- get_upper.i(x = x, lavmodel = lavmodel, x.i = x.i)
  tmp2 <- diag(dnorm(upper.i))
  out <- tmp2 %*% tmp1
  out
}

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

get_w_exp.i <- function(x, lavmodel = NULL, x.i) {
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

  nvar <- this.model@nvar[1]
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
  w_exp.i
}

grad_pbivnorm <- function(x) {
  a = x[1]
  b = x[2]
  rho = x[3]
  R <- 1 - rho*rho
  da <- dnorm(a) * pnorm((b - rho*a)/sqrt(R))
  db <- dnorm(b) * pnorm((a - rho*b)/sqrt(R))
  drho <- lavaan:::dbinorm(u = a, v = b, rho = rho)
  c(da, db, drho)
}


get_scores <- function(x) {
  # raw data
  y_obs <- fit0@Data@X[[1]] - 1
  y_mu  <- colMeans(y_obs)
  eXo   <- fit0@Data@eXo[[1]]
  nvar  <- ncol(y_obs)
  pstar <- nvar * (nvar - 1 ) / 2

  SCORES <- matrix(0, nrow = nrow(y_obs), ncol = length(x))
  for(i in 1:nrow(y_obs)) {

    cat("i = ", i, "\n")

    all_i <- get_y_yy_exp.i(x = x, lavmodel = fit0@Model, x.i = eXo[i,])
	y_exp.i <- all_i[1:nvar]
	w_exp.i <- all_i[-(1:nvar)]
	y_obs.i <- y_obs[i,]

    ### how to center y?
    # option 1:
	#y_obs_c.i <- y_obs.i - y_mu 
	# option 2:
	#mu_i <- get_mu_star.i(x = x, lavmodel = fit0@Model, x.i = eXo[i,])
	#y_obs_c.i <- y_obs.i - mu_i
	# option 3:
	y_obs_c.i <- y_obs.i - y_exp.i

    w_obs.i <- get_w_obs.i(y = y_obs_c.i)
  
    # covy
    covy <- lav_matrix_vech_reverse(w_exp.i, diagonal = FALSE)
    diag(covy) <- y_exp.i * (1 - y_exp.i)

    # covw
    covw <- diag(pstar)
	### double check: paper uses E(w_obs.i^2) !! (but what is it?)
    diag(covw) <- w_obs.i^2 - w_exp.i^2

    # B_i
    B_i <- lav_matrix_bdiag(covy, covw)
    B_i.inv <- solve(B_i)

    # Delta
    Delta.i <- numDeriv:::jacobian(func = get_y_yy_exp.i, x = x,
                                   lavmodel = fit0@Model, x.i = eXo[i,])

    error_y <- y_obs.i - y_exp.i
    error_w <- w_obs.i - w_exp.i
    both_error <- c(error_y, error_w)

    SCORES[i, ] <- drop( t(Delta.i) %*% B_i.inv %*% as.matrix(both_error) )
  }

  SCORES
}

