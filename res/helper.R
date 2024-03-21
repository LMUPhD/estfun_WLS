
#install latest lavaan package!
#detach('package:lavaan', unload=TRUE)
#devtools::install_github("yrosseel/lavaan",force=TRUE, build_opts = c())
#library(lavaan)

################################################################################
setwd("C:\\Users\\classe\\Desktop\\Diss\\Paper3\\estfun_WLS")
source("estfun_WLS_complex.R") #simple
source("application\\multi_simu.R")

################################################################################
model = '
  Eta3 =~ simuvar7 + simuvar8 + simuvar9
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6' 
fits_random <- datagen(model = model, schwellen = 4, ID=1000, times=1, items=3, latvar = 3) #multi
simu=fits_random[["data"]][["data1"]]

#mirt
#model_mirt <- mirt::mirt.model('
#  Eta1 = 1-3
#  Eta2 = 4-6
#  Eta3 = 7-9
#  COV=Eta1*Eta2*Eta3') 
#fit_mirt <- mirt::mirt(data=simu[,1:9], model=model_mirt, itemtype="graded",method="MHRM" )  


#lavaan
fit_ord <- lavaan::cfa(model, data = simu, ordered = TRUE, estimator = "WLS",std.lv=TRUE )
colSums(estfun.WLS(fit_ord))
object = fit_ord
#summary(fit_ord)
#lavaan::estfun.lavaan(fit_ord)






################################################################################
######################### Compute Delta new?? ##################################
################################################################################
params <- lav_object_inspect_coef(object,type = "free", add.labels = F)
params <- lav_object_inspect_coef(object,type = "user", add.labels = F)



compute.moments <- function(params) {
  GLIST <- lav_model_x2GLIST(lavmodel = lavmodel, x=params, type="free")
  Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, GLIST = GLIST)
  polychors = Sigma.hat[[1]]
  th = as.vector(GLIST[["tau"]])
  th.pr = VGAM::probitlink( th*-1,inverse=T) 
  mus = unlist(lapply(1:nvar, function(x) get_mus(x, th, lv, nvar, catvals)   ))
  combs = rbind(  combn(1:nvar,2), lavaan::lav_matrix_vech(polychors,diagonal=FALSE) ) 
  joint_exps = apply(combs, 2L, function(x) get_joint_exp(x, X, th, lv, nvar, catvals)  ) #E(y1y2)
  sigma =  joint_exps - t(  lavaan::lav_matrix_vech(tcrossprod(mus) ,diagonal=FALSE) )  #E(y1y2)-mu1mu2
  return(c(th,sigma))
}
 
Delta <- numDeriv::jacobian(func=compute.moments, x = params)
#???? compute moments --> irgendeine Funktion, die geschätze Parameter aufnimmt und c(th,sigma) ausgibt... oder umgekehrt?







################################################################################
################################################################################
simu$rand1 <- replicate(nrow(simu),sample(c(1,2,3,4,5),1,prob=c(0.1,0.2,0.4,0.2,0.1))) #ordinal
simu$rand2 <- round(runif(nrow(simu),min=100,max=400),2)                               #metric
simu$rand3 <- replicate(nrow(simu),sample(c(1,2,3,4),1,prob=c(0.15,0.35,0.35,0.15)))   #nominal
#covariates


#mirt
fluc_ord <- strucchange::gefp(fit_mirt, fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand1"]) #ordinal/categorical
fluc_num <- strucchange::gefp(fit_mirt, fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand2"]) #numeric
fluc_cat <- strucchange::gefp(fit_mirt, fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand3"]) #nominal


#lavaan
fluc_ord <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand1"]) #ordinal/categorical
fluc_num <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand2"]) #numeric
fluc_cat <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand3"]) #nominal


#both
p_maxBB    <- strucchange::sctest(fluc_num, functional = strucchange::maxBB)$p.value #dm
p_meanL2BB <- strucchange::sctest(fluc_num, functional = strucchange::meanL2BB)$p.value #cvm 
p_sumlm    <- strucchange::sctest(fluc_num, functional=  strucchange::supLM(from = 0.15, to = NULL))$p.value #suplm
p_ordL2BB  <- strucchange::sctest(fluc_ord, functional = strucchange::ordL2BB(fluc_ord))$p.value[1]  #LM #istead of fluc1 ... factor(simu[,"rand1"])
p_ordwmax  <- strucchange::sctest(fluc_ord, functional = strucchange::ordwmax(fluc_ord))$p.value[1]  #wdmo 
p_catl2bb  <- strucchange::sctest(fluc_cat, functional = strucchange::catL2BB(fluc_cat))$p.value  #LMuo



