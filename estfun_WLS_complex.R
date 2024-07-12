
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")
source("support.R")


estfun.WLS <- function(object){
 
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
  
  Score.mat <- matrix(NA, ntot, npar) #empty matrix
  
  
  
  ################################################################################
  ################################# WLS ##########################################
  ################################################################################
  
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

  
  ### combine matrices
  Score.mat = t( t(Delta) %*% W.inv %*% t(e)  ) 
  
  
  
  ################################################################################
  
  
  #Score.mat <- (-1/ntot) * Score.mat #scaling
  
  # provide column names
  colnames(Score.mat) <- names(lav_object_inspect_coef(object,
                                                       type = "free", add.labels = TRUE))
  
  return(Score.mat)
  
}
