
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")
source("estfun_WLS.R")
source("multi_simu.R")

################################################################################
######################### Ordinal Data ######################################### 
################################################################################



model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'


fits_random <- datagen(model = model_lav,rmsea_cutoff = .05,ID=1000,items=3,latvar=3,schwellen=4) 
simu=fits_random[["data"]][["data1"]]



################################################################################
############################### Models ######################################### 
################################################################################


#lavaan
fit_ord <- lavaan::cfa(model_lav, data = simu, ordered = TRUE, estimator = "WLS",std.lv=F ) #estimate variances and covariances of latent variables
colSums(estfun.WLS(fit_ord)) 



#mirt
model_mirt <- mirt::mirt.model('
Eta1 = 1-3
Eta2 = 4-6
Eta3 = 7-9
Beta2 = 2,5,8
Beta3 = 3,6,9
COV=Eta1*Eta2*Eta3*Beta2*Beta3
FIXED = (2-9, a1)
FIXED = (1-3,5-9, a2)
FIXED = (1-6,8,9, a3)
FIXED = (1,3-9, a4)
FIXED = (1,2,4-9, a5)') #????
fit_mirt <- mirt::mirt(data=simu, model=model_mirt, itemtype="graded",method='MHRM' ) #pars="values" #'EM',TOL = NaN 
mirt::M2(fit_mirt)


####select a simpler model!!

################################################################################
########################## try strucchange #####################################
################################################################################
simu$rand1 <- replicate(nrow(simu),sample(c(1,2,3,4,5),1,prob=c(0.1,0.2,0.4,0.2,0.1))) #ordninal
simu$rand2 <- round(runif(nrow(simu),min=0,max=200),2)                               #numeric (univ)


flucs1 <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand1"]) #ordinal
flucs2 <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand2"]) #numeric

strucchange::sctest(flucs2, functional = strucchange::maxBB)$p.value #dm    #--> passt... andere stats funktionieren mit >25 parametern nicht
strucchange::sctest(flucs1, functional = strucchange::ordwmax(flucs1))$p.value[1]  #wdmo
