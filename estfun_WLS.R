
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
  lv = lavdata@ov[["nlev"]]
  
  #moments <- lavaan::fitted(object)
  N1 <- 1
  X <- lavdata@X[[1]]
  
  Score.mat <- matrix(NA, ntot, npar) #empty matrix
  
  
  
  ################################################################################
  ################################# WLS ##########################################
  ################################################################################
  

  #dummies
  Xd = do.call(cbind, lapply(1:nvar, function(i) doDummySingleVar(X,lv,ntot,i)  )) #Problems with doDummySingleVar...
  
  
  ###e1
  musd = colMeans(Xd)
  e1 = t( apply(Xd, 1L, function(x) x-musd ) )
  
  ###e2 
  mus = colMeans(X)
  y_minus_mu = t( apply(X, 1L, function(x) x - mus ) ) 
  s_vech = t(apply(y_minus_mu, 1L, function(i){    lavaan::lav_matrix_vech(tcrossprod(i) ,diagonal=FALSE) })) #s=c( (y1-mu1)(y2-mu2)....
  sigma = colMeans(s_vech)
  
  
  e2 = t( apply(s_vech, 1L, function(x) x - sigma ) ) 
  
  
  
  ###e
  e = cbind(e1,e2)
  
  #weigthing matrix
  if(object@call[["estimator"]]=="WLS"){
    W = lavsamplestats@WLS.V[[1]] #WLS.V is already inverted W!
  } else if(object@call[["estimator"]]=="DWLS"){
    W = matrix(0,ncol=length(c(th,sigma)), nrow=length(c(th,sigma)))
    diag(W) = lavsamplestats@WLS.VD[[1]] #WLS.V is already inverted W!
  }

  
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
