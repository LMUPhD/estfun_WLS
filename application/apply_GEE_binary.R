setwd("//dss/dsshome1/0B/ra35tik2/einleitung")
#source("univ_simu.R") #simulate data (univariate model)
source("multi_simu.R") #simulate data (multivariate model)

################################################################################
########################## GEE estimation non-binary ###########################
################################################################################
model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'


fits_random <- datagen(model=model_lav, schwellen = 1, ID=2000, times=1, items=3,latvar=3) 
Data = fits_random[["data"]][["data1"]]
#model = fits_random[["model"]][[1]] #only works for univariate model


################################################################ My Code
fit0 <- lavaan::sem(model_lav, data = Data, do.fit = FALSE, estimator = "WLS",ordered = T)
coef(fit0)
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

# estimate parameters using GEE
library(pbivnorm)
max_iter <- 200L
old.x <- start.x


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







#### test GEE scores
fit.wls <- lavaan::cfa(model_lav, data = Data, ordered = TRUE, estimator = "WLS", std.lv=F )
coef(fit.wls) - new.x 


################################################################################
############################### Evaluate #######################################
################################################################################
SC_GEE = SCORES; colSums(SC_GEE)
SC_WLS = lavaan::lavScores(fit.wls); colSums(SC_WLS)
round(diag(cor(SC_GEE, SC_WLS)), 5)   






