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
lavdata        <- object@Data
lavmodel  <- object@Model
lavsamplestats <- object@SampleStats
nvar <- ncol(lavsamplestats@cov[[1]])
lv = lavdata@ov[["nlev"]]
X <- lavdata@X[[1]]
catvals = lapply(1:nvar, function(x)  as.numeric(names(table(X[,x]))) )
#model context for compute.moments
model.context = list()
model.context[["m.free.idx"]] <- lavmodel@m.free.idx
model.context[["x.free.idx"]] <- lavmodel@x.free.idx
model.context[["GLIST.template"]] <- lavmodel@GLIST
model.context[["isSymmetric"]] <- lavmodel@isSymmetric
model.context[["lv"]] = lv
model.context[["nvar"]] = nvar
model.context[["catvals"]] = catvals

#params
idx <- which(object@ParTable$free > 0L)
est <- object@Fit@est
params <- est[idx]


#R
source("support.R")
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
moments_r = compute.moments(params,lavmodel)

#C++
Rcpp::sourceCpp("Rcpp/test_stuff.cpp") #shit...
moments_c = compute_moments(params,model.context)
moments_c
