setwd("C:\\Users\\classe\\Desktop\\Diss\\Paper3\\estfun_WLS")
source("application/multi_simu.R") #simulate data (multivariate model)
model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'
model = model_lav
fits_random <- datagen(model=model_lav, schwellen = 4, ID=2000, times=1, items=3,latvar=3) 
Data = fits_random[["data"]][["data1"]]

#estimate parmeters
fit.wls <- lavaan::cfa(model_lav, data = Data, ordered = TRUE, estimator = "WLS", std.lv=F )
object = fit.wls



################################################################################
################################################################################
null_model <- lavaan::cfa(model_lav, data = Data, ordered = TRUE, estimator = "WLS", std.lv=F, do.fit=F )
#only starting values

source("Rcpp/support_rcpp.R")
params <- lav_object_inspect_coef(object,type = "free", add.labels = F) #fitted params
lavmodel       <- null_model@Model #lavmodel from null_model

#code in jeweiliger funtion ausführen...

#####  R
res_long <- lav_model_x2GLIST(lavmodel = lavmodel, x=params, type="free")
Sigma.hat_long <- computeSigmaHat(lavmodel = lavmodel, GLIST = res_long)



#simpler...
m.free.idx <- lavmodel@m.free.idx
x.free.idx <- lavmodel@x.free.idx
GLIST.template <- lavmodel@GLIST
isSymmetric <- lavmodel@isSymmetric

simple_x2GLIST <- function(x, 
                           m.free.idx,     # list of matrix indices for each matrix
                           x.free.idx,     # list of parameter vector indices for each matrix
                           GLIST.template, # template with correct dimensions
                           isSymmetric = NULL) {  # optional for symmetric matrices
  
  GLIST <- GLIST.template  # start with template
  
  for (mm in 1:length(GLIST)) {
    # Skip if no free parameters for this matrix
    if (length(m.free.idx[[mm]]) == 0) next
    
    # Assign values from parameter vector to matrix
    GLIST[[mm]][m.free.idx[[mm]]] <- x[x.free.idx[[mm]]]
    
    # Make symmetric if needed (for psi/theta matrices)
    if (!is.null(isSymmetric) && isSymmetric[mm]) {
      GLIST[[mm]][upper.tri(GLIST[[mm]])] <- t(GLIST[[mm]])[upper.tri(GLIST[[mm]])]
    }
  }
  
  return(GLIST)
}

res_short= simple_x2GLIST(
  x = params,
  m.free.idx = m.free.idx,
  x.free.idx = x.free.idx,
  GLIST.template = GLIST.template,
  isSymmetric = isSymmetric
)


identical(res_short, lavmodel@GLIST)
identical(res_long, lavmodel@GLIST) #theta not filled up with estimates because they're not in params...
View(null_model@Model@GLIST)
View(object@Model@GLIST)


##### C++ 
Rcpp::sourceCpp("Rcpp/test_stuff.cpp")
res_c = x2GLIST_withtheta(x = params,
                          m_free_idx = m.free.idx,
                          x_free_idx = x.free.idx,
                          GLIST_template = GLIST.template,
                          isSymmetric = isSymmetric)
identical(res_long,res_c)


## compute sigma_hat
lambda = res_c$lambda
psi = res_c$psi
theta = res_c$theta

sigma_hat = compute_SigmaHat(lambda,psi,theta)

identical(Sigma.hat_long[[1]],sigma_hat)

Scores_c = compute_scores(Delta,W.inv,e)
