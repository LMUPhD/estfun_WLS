setwd("C:\\Users\\classe\\Desktop\\Diss\\Paper3\\estfun_WLS")
#source("application/univ_simu.R") #simulate data (univariate model)
source("application/multi_simu.R") #simulate data (multivariate model)


################################################################################
########################## GEE estimation non-binary ###########################
################################################################################

model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'
model = model_lav
fits_random <- datagen(model=model_lav, schwellen = 6, ID=1000, times=1, items=3,latvar=3, rmsea_cutoff = 0.9) 
Data = fits_random[["data"]][["data1"]]



#fits_random <- datagen(schwellen = 4, ID=1000, times=1, items=5) # categories
#real_params = c(fits_random$betas$beta1[-1],as.vector(fits_random$kappas$kappa1),fits_random$vars$var1)
#Data = fits_random[["data"]][["data1"]]
##Data = Data-1
#model = fits_random[["model"]][[1]]



################################################################ My Code
fit0 <- lavaan::sem(model, data = Data, do.fit = FALSE, estimator = "WLS",ordered = T)
coef(fit0)
start.x <- coef(fit0)

Rcpp::sourceCpp("Rcpp/gee_rcpp_support.cpp")

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
  npar <- sum(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
  nvar <- ncol(lavsamplestats@cov[[1]])
  
  #params
  idx <- which(object@ParTable$free > 0L)
  est <- object@Fit@est
  params <- est[idx]
  
  #moments <- lavaan::fitted(object)
  X <- lavdata@X[[1]]
  Score.mat <- matrix(NA, ntot, npar) #empty matrix
  
  
  lv = lavdata@ov[["nlev"]]
  catvals = get_unique_values_per_column(X)                                     
  
  
  #model context for compute.moments
  model.context = list()
  model.context[["m.free.idx"]] <- lavmodel@m.free.idx
  model.context[["x.free.idx"]] <- lavmodel@x.free.idx
  model.context[["GLIST.template"]] <- lavmodel@GLIST
  model.context[["isSymmetric"]] <- lavmodel@isSymmetric
  model.context[["lv"]] = lv
  model.context[["nvar"]] = nvar
  model.context[["catvals"]] = catvals
  
  
  ################################################################################
  ################################# GEE ##########################################
  ################################################################################
  
  #polychoric corr
  polychors = object@Fit@Sigma.hat[[1]] 
  
  th = object@Fit@TH[[1]]
  th.pr = pbv_rcpp_pnorm( th*-1)                                       #add the diagonal of the model implied matrix!                   
  
  #dummies
  Xd = createDummies(X,lv)
  
  
  ###e1
  e1 = minus_mat(Xd,th.pr)                     
  
  ###e2 
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
  Delta <- compute_jacobian_fast(params,model.context)
  
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


#### Compute model fit test statistics
nsatpar = ncol(W)
#R = sum(W != 0)
#GF = 1 - sum( apply(e, 1, FUN = function(ei){t(ei) %*% W.inv %*% ei}) ) / ((ntot-1)*nsatpar)
#this values coult be bootstrapped... that too cost intensiv though...

# compute naive deviance 
Tl = sum( apply(e, 1, FUN = function(ei){t(ei) %*% W.inv %*% ei}) ) / ntot

#compute P0 and P1
P0 = DBiD / ntot
P1 = (t(Delta) %*% W.inv %*% t(e)%*%e %*% W.inv %*% Delta) / ntot
P = solve(P0)%*%P1

#compute d_bar
r = nsatpar - npar
d_bar = tr(P)/r
#d_bar = tr(P1)/tr(P0) #doesnt work...

#compute test stat
Ta = Tl/d_bar

#compute p-value
pvalue_gee = pchisq(q=Ta, df=r, lower.tail=FALSE)
rmsea_gee = sqrt(pmax((Ta / ntot) / r - 1 / ntot, 0)) 


#### Compute standard errors





################################################################################
############################### Evaluate #######################################
################################################################################

#### compare GEE scores
fit.wls <- lavaan::cfa(model_lav, data = Data, ordered = TRUE, estimator = "WLS", std.lv=F )
coef(fit.wls) - new.x 

model_mirt <- mirt::mirt.model('
  Eta1 = 1-3
  Eta2 = 4-6
  Eta3 = 7-9
  COV=Eta1*Eta2*Eta3')
fit.mirt <- mirt::mirt(Data, model=model_mirt, itemtype="graded",method="MHRM" )  



SC_GEE = SCORES; colSums(SC_GEE)
SC_WLS = lavaan::lavScores(fit.wls); colSums(SC_WLS)
round(diag(cor(SC_GEE, SC_WLS)), 5)   





