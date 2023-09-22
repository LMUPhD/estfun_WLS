
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")
source("estfun_WLS.R")


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
######################### Ordinal Data ##################################### 
################################################################################

source("univ_simu.R")

model = 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5'
fits_random <- datagen(model = model, 
                       schwellen = 2,   #je mehr Schwellen, desdo ungenauer wird es...
                       rmsea_cutoff = .05,
                       ID=500,
                       times=1,
                       items=5) 

simu=fits_random[["data"]][["data1"]][[1]]

fit_ord <- lavaan::cfa(model, data = simu, ordered = TRUE, estimator = "WLS" ) 
WLS_scores <- estfun.WLS(fit_ord)
colSums(WLS_scores)


#compare with ML
simu_ml=fits_random[["data"]][["data1"]][[2]]

fit_num <- lavaan::cfa(model, data = simu_ml, estimator = "ML" ) 
ML_scores <- lavaan::lavScores(fit_num)
colSums(ML_scores)







################################################################################
############################## Sandbox ######################################### 
################################################################################

