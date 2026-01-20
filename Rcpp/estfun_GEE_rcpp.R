
Rcpp::sourceCpp("Rcpp/gee_rcpp_support.cpp")


estfun.GEE <- function(object){
 
  compute.moments <- function(params, lavmodel = NULL) {
    
    m.free.idx <- lavmodel@m.free.idx
    x.free.idx <- lavmodel@x.free.idx
    GLIST.template <- lavmodel@GLIST
    isSymmetric <- lavmodel@isSymmetric
    
    GLIST <- x2GLIST_withtheta(x = params,
                               m_free_idx = m.free.idx,
                               x_free_idx = x.free.idx,
                               GLIST_template = GLIST.template,
                               isSymmetric = isSymmetric)
    sigma_hat = compute_SigmaHat(GLIST$lambda,GLIST$psi,GLIST$theta)
    th = as.vector(GLIST[["tau"]])
    
    th.pr = pbv_rcpp_pnorm(th*-1)
    mus = apply_get_mus(th, lv, nvar,  catvals)
    combs = create_combs(nvar,sigma_hat)         
    joint_exps = apply_get_joint_exp(combs, th, lv, nvar, catvals)  #E(y1y2) 
    sigma =  compute_sigma(joint_exps,mus)  #E(y1y2)-mu1mu2
    return(c(th.pr,sigma))
  }
   
  # shortcuts
  lavdata        <- object@Data
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavoptions     <- object@Options
  
  
  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)
  npar <- sum(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
  nvar <- ncol(lavsamplestats@cov[[1]])
  
  #params
  idx <- which(object@ParTable$free > 0L)
  est <- object@Fit@est
  params <- est[idx]
  
  #moments <- lavaan::fitted(object)
  X <- lavdata@X[[1]]
  Score.mat <- matrix(NA, ntot, npar) #empty matrix
  
  
  ################################################################################
  ################################# GEE ##########################################
  ################################################################################
  
  #polychoric corr
  polychors = object@Fit@Sigma.hat[[1]] 

  th = object@Fit@TH[[1]]
  th.pr = pbv_rcpp_pnorm( th*-1)                                       #add the diagonal of the model implied matrix!                   
  
  #dummies
  lv = lavdata@ov[["nlev"]]
  Xd = createDummies(X,lv)
  
  
  ###e1
  e1 = minus_mat(Xd,th.pr)                     

  ###e2 
  catvals = get_unique_values_per_column(X)                                     
  mus = apply_get_mus(th, lv, nvar,  catvals)
  y_minus_mu = minus_mat(X, mus)                         
  combs = create_combs(nvar,polychors) 
  joint_exps = apply_get_joint_exp(combs, th, lv, nvar, catvals) #E(y1y2)  
  sigma =  compute_sigma(joint_exps,mus)  #E(y1y2)-mu1mu2
  s_vech = compute_vech_by_row(y_minus_mu) #s=c( (y1-mu1)(y2-mu2)....
  
  e2 = minus_mat(s_vech,sigma)    
  
  
  ###e
  e = cbind(e1,e2)
  
  ###weigthing matrix
  
  #get sigma_indi
  mat1 = create_upper_weight_matrix(th, lv, nvar, polychors)                                           

  diag2 = colMeans(s_vech^2) - sigma^2
  mat2 = matrix(diag(c(diag2)),ncol=length(diag2) )

  W = lav_matrix_bdiag(mat1,mat2)
  W.inv = solve(W)

  
  #Delta
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
