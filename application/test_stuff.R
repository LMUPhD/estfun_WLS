setwd("C:\\Users\\classe\\Desktop\\Diss\\Paper3\\estfun_WLS")
source("application\\univ_simu.R")
library(pbivnorm)


################################################################################
############################# GEE estimation ###################################
################################################################################
fits_random <- datagen(schwellen = 2, ID=1000, times=1, items=3) # categories
Data = fits_random[["data"]][["data1"]]
#Data = Data-1
model = fits_random[["model"]][[1]]



################################################################ My Code
fit0 <- lavaan::sem(model, data = Data, do.fit = FALSE, estimator = "WLS",ordered = T)
coef(fit0)
start.x <- coef(fit0)

compute.moments <- function(params, lavmodel = NULL) {
  GLIST <- lav_model_x2GLIST(lavmodel = lavmodel, x=params, type="free")
  Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, GLIST = GLIST)
  polychors = Sigma.hat[[1]]
  th = as.vector(GLIST[["tau"]])
  th.pr = pnorm(th*-1)
  mus = unlist(lapply(1:nvar, function(x) get_mus(x, th, lv, nvar, catvals)   ))
  combs = rbind(  combn(1:nvar,2), lavaan::lav_matrix_vech(polychors,diagonal=FALSE) )           
  joint_exps = apply(combs, 2L, function(x) get_joint_exp(x, th, lv, nvar, catvals)  ) #E(y1y2) 
  sigma =  joint_exps - t(  lavaan::lav_matrix_vech(tcrossprod(mus) ,diagonal=FALSE) )  #E(y1y2)-mu1mu2
  return(c(th.pr,sigma))
}
source("support.R")

# estimate parameters using GEE
max_iter <- 200L
old.x <- start.x
#old.x[5:9] <- 0

for(iter in seq_len(max_iter)) {
  
  # verbose
  cat("iter = ", iter, "\n")
  print(old.x)
  cat("\n")
  
  # insert parameters in model matrices
  object <- sem(model, data = Data, do.fit = FALSE, start = old.x, estimator = "WLS",ordered = T)
  
  
  # shortcuts
  lavdata        <- object@Data
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavoptions     <- object@Options
  
  
  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)
  npar <- lav_object_inspect_npar(object)
  nvar <- ncol(lavsamplestats@cov[[1]])
  
  
  #moments <- lavaan::fitted(object)
  X <- lavdata@X[[1]]

  #polychoric corr
  polychors = object@Fit@Sigma.hat[[1]] 
  
  th = object@Fit@TH[[1]]
  th.pr = pnorm( th*-1)                                       #add the diagonal of the model implied matrix!                   
  
  #dummies
  lv = lavdata@ov[["nlev"]]
  Xd = do.call(cbind, lapply(1:nvar, function(i) doDummySingleVar(X,lv,ntot,i)  )) 
  
  
  ###e1
  e1 = t( apply(Xd, 1L, function(x) x-th.pr) )                     
  
  
  ###e2 
  catvals = lapply(1:nvar, function(x)  as.numeric(names(table(X[,x]))) )
  mus = unlist(lapply(1:nvar, function(x) get_mus(x, th, lv, nvar, catvals)   ))
  y_minus_mu = t( apply(X, 1L, function(x) x - mus ) ) 
  
  combs = rbind(  combn(1:nvar,2), lavaan::lav_matrix_vech(polychors,diagonal=FALSE) ) 
  joint_exps = apply(combs, 2L, function(x) get_joint_exp(x, th, lv, nvar, catvals)  ) #E(y1y2)
  sigma =  joint_exps - t(  lavaan::lav_matrix_vech(tcrossprod(mus) ,diagonal=FALSE) )  #E(y1y2)-mu1mu2
  
  s_vech = t(apply(y_minus_mu, 1L, function(i){    lavaan::lav_matrix_vech(tcrossprod(i) ,diagonal=FALSE) })) #s=c( (y1-mu1)(y2-mu2)....
  
  
  e2 = t( apply(s_vech, 1L, function(x) x - sigma ) ) 
  
  
  
  ###e
  e = cbind(e1,e2)
  
  ###weigthing matrix
  
  #get sigma_indi
  seqnc=c()
  for(l in seq(nvar)){ seqnc = c(seqnc,sapply(head(seq(lv[l]),-1), function(y) paste(l,y )  ))}
  combs_indi = combn(seqnc,2)
  sigma_indi = apply(combs_indi, 2L, function(x) get_sigmas_indi(x, th, lv, nvar, polychors)  ) 
  
  mat1 = matrix(0, ncol=length(th),nrow=length(th))
  mat1[lower.tri(mat1)] = sigma_indi 
  mat1[upper.tri(mat1)] <- t(mat1)[upper.tri(t(mat1))] # 
  diag(mat1) = th.pr*(1-th.pr)        #mus for indicators in diagonal                                           
  
  diag2 = colMeans(s_vech^2) - sigma^2
  mat2 = matrix(diag(c(diag2)),ncol=length(diag2) )
  
  W = lav_matrix_bdiag(mat1,mat2)
  W.inv = solve(W)
  
  
  #Delta
  params <- lav_object_inspect_coef(object,type = "free", add.labels = F)
  Delta <- numDeriv::jacobian(func=compute.moments, x = params, lavmodel = lavmodel)
  
  
  # Scores
  SCORES <- t( t(Delta) %*% W.inv %*% t(e)  ) 
  sum.SCORES <- colSums(SCORES)
  
  DBiD <- ntot * (t(Delta) %*% W.inv %*% Delta)
  
  # update parameters
  step <- 0.6
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
fit.wls <- lavaan::cfa(model, data = Data, ordered = TRUE, estimator = "WLS", std.lv=F )
coef(fit.wls) - new.x 

fit.dwls <- lavaan::cfa(model, data = Data, ordered = TRUE, estimator = "DWLS", std.lv=F )
coef(fit.dwls) - new.x #Binary: DWLS estimates almost identical parameters





################################################################################
############################### Evaluate #######################################
################################################################################
source("estfun_GEE.R")


SC <- get_scores(new.x) 
fit1 <- sem(model, data = Data, start = new.x, do.fit = FALSE, estimator = "WLS", ordered = T) #egal ob WLS oder DWLS
fit2 <- sem(model, data = Data, start = coef(fit.dwls), do.fit = FALSE, estimator = "WLS", ordered = T)

SC1 <- estfun.GEE(fit1)
SC2 <- estfun.GEE(fit2) 
SC3 <- estfun.GEE(fit.dwls)
SC4 <- lavaan::lavScores(fit1)
SC5 <- lavaan::lavScores(fit2)
SC6 <- lavaan::lavScores(fit.wls)
#SC=SC1; SC2 = SC3

round(diag(cor(SC1, SC2)), 5) #perfekt!
round(diag(cor(SC2, SC3)), 5) #perfekt! 
round(diag(cor(SC1, SC3)), 5) #perfekt!  


round(diag(cor(SC1, SC4)), 5) #ziemich unterschiedlich...  
round(diag(cor(SC1, SC5)), 5) #ziemich unterschiedlich...  








################################################################### Yves code
fit0 <- sem(model, data = Data, do.fit = FALSE, estimator = "WLS",ordered = T)
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
#old.x[5:9] <- 0

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





################################################################################
################################ Sandbox #######################################
################################################################################

#different test results even though scores correlate perfectly?


#fit1@Fit@Sigma.hat[[1]] 
#fit.dwls@Fit@Sigma.hat[[1]] 
#fit.wls@Fit@Sigma.hat[[1]] 
