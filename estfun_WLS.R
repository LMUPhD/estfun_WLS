
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")
source("support.R")


estfun.WLS <- function(object){
  
  # shortcuts
  lavdata        <- object@Data
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavoptions     <- object@Options
  
  
  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)
  npar <- lav_object_inspect_npar(object)
  nvar <- ncol(lavsamplestats@cov[[1]])
  
  
  moments <- fitted(object)
  N1 <- 1
  X <- lavdata@X[[1]]
  
  Score.mat <- matrix(NA, ntot, npar) #empty matrix
  
  
  
  ################################################################################
  ################################# WLS ##########################################
  ################################################################################
  
  #polychoric corr
  polychors = object@SampleStats@cov[[1]]
  #s = lavaan::fitted(object)$cov  #doesnt work...
  
  th = lavsamplestats@th[[1]]
  th.pr = VGAM::probitlink( th*-1,inverse=T) 
  
  #dummies
  lv = lavdata@ov[["nlev"]]
  Xd = do.call(cbind, lapply(1:nvar, function(i) doDummySingleVar(X,lv,ntot,i)  )) #Problems with doDummySingleVar...
  
  
  ###e1
  e1 = t( apply(Xd, 1L, function(x) x-th.pr ) ) 
  
  
  ###e2 
  catvals = lapply(1:nvar, function(x)  as.numeric(names(table(X[,x]))) )
  mus = unlist(lapply(1:nvar, function(x) get_mus(x, th, lv, nvar, catvals)   ))
  y_minus_mu = t( apply(X, 1L, function(x) x - mus ) ) 
  
  combs = rbind(  combn(1:nvar,2), lav_matrix_vech(polychors,diagonal=FALSE) ) 
  joint_exps = apply(combs, 2L, function(x) get_joint_exp(x, X, th, lv, nvar, catvals)  ) #E(y1y2)
  sigma =  joint_exps - t(  lav_matrix_vech(tcrossprod(mus) ,diagonal=FALSE) )  #E(y1y2)-mu1mu2
  
  s_vech = t(apply(y_minus_mu, 1L, function(i){    lav_matrix_vech(tcrossprod(i) ,diagonal=FALSE) })) #s=c( (y1-mu1)(y2-mu2)....
  
  e2 = t( apply(s_vech, 1L, function(x) x - sigma ) ) 
  
  
  
  ###e
  e = cbind(e1,e2)
  
  #weigthing matrix
  W = lavsamplestats@WLS.V[[1]] 
  
  #Delta
  Delta <- computeDelta(lavmodel = lavmodel)[[1]] #should also work for WLS...
  
  ### combine matrices
  Score.mat = t( t(Delta) %*% W %*% t(e)  )  #works without mu*...
  
  
  
  ################################################################################
  
  
  #Score.mat <- (-1/ntot) * Score.mat #scaling
  
  # provide column names
  colnames(Score.mat) <- names(lav_object_inspect_coef(object,
                                                       type = "free", add.labels = TRUE))
  
  return(Score.mat)
  
}
