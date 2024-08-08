source("application\\univ_simu.R") #simulate data (univariate model)


################################################################################
########################## GEE estimation non-binary ###########################
################################################################################
fits_random <- datagen(schwellen = 3, ID=1000, times=1, items=3) # categories
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
  Delta <- numDeriv::jacobian(func=compute.moments, x = params, lavmodel = lavmodel) #takes a long time!
  #Delta <- computeDelta(lavmodel = lavmodel)[[1]] #also works for WLS but is not correct according to papers
  
  
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


SC1 = SCORES; colSums(SC1)
SC2 = estfun.GEE(fit.wls); colSums(SC2) #doesn't sum up to one!
SC3 = lavaan::lavScores(fit.wls); colSums(SC3)

round(diag(cor(SC1, SC2)), 5) #near-perfect correlation
round(diag(cor(SC2, SC3)), 5)  
round(diag(cor(SC1, SC3)), 5)   














