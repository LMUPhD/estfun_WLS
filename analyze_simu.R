
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")
source("estfun_WLS.R")

################################################################################
######################### Ordinal Data ##################################### 
################################################################################

source("univ_simu.R")

fits_random <- datagen(schwellen = 2,                                           #je mehr Schwellen, desdo ungenauer wird es...
                       ID=1000,
                       times=1,
                       items=4) 

simu=fits_random[["data"]][["data1"]][[1]]
model= fits_random[["model"]][[1]]

fit_ord <- lavaan::cfa(model, data = simu, ordered = TRUE, estimator = "WLS" ) 
WLS_scores <- estfun.WLS(fit_ord)
colSums(WLS_scores)

#muss mu anders berechnet werden?? E(Y1) = sum y1*P(y1y2)

################################################################################
######################### Dichotomous Data ##################################### 
################################################################################

## mirt:
library(mirt)
dt <- mirt::expand.table(LSAT7)
fit_mirt <- mirt(dt[,1:5], 1, itemtype = "Rasch", SE = TRUE)
mirt_scores <- estfun.AllModelClass(fit_mirt)

####################
library(lavaan)


model = 'Eta1 =~ 1*Item.1 + 1*Item.2 + 1*Item.3 + 1*Item.4 + 1*Item.5'

fit_ord <- lavaan::cfa(model, data = dt, ordered = TRUE, estimator = "WLS" ) #parameterization = "theta"
WLS_scores <- estfun.WLS(fit_ord)

object = fit_ord



################################################################################
########################## try strucchange #####################################
################################################################################




#compare with ML
simu_ml=fits_random[["data"]][["data1"]][[2]]

fit_num <- lavaan::cfa(model, data = simu_ml, estimator = "ML" ) 
ML_scores <- lavaan::lavScores(fit_num)
colSums(ML_scores)







################################################################################
############################## Sandbox ######################################### 
################################################################################

