
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

fit_ord <- cfa(model, data = dt, ordered = TRUE, estimator = "WLS" ) #parameterization = "theta"
WLS_scores <- estfun.WLS(fit_ord)

object = fit_ord

################################################################################
######################### Ordinal Data ##################################### 
################################################################################

source("univ_simu.R")

model = 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5'
fits_random <- datagen(model = model, 
                       schwellen = 1,
                       rmsea_cutoff = .05,
                       ID=500,
                       times=1,
                       items=5) 

simu=fits_random[["data"]][["data1"]][[1]]

fit_ord <- cfa(model, data = simu, ordered = TRUE, estimator = "WLS" ) 
WLS_scores <- estfun.WLS(fit_ord)

################################################################################
############################## Sandbox ######################################### 
################################################################################


mvtnorm::pmvnorm(lower= 0.07526986, upper = 1.205072 , sigma = 0.3807874)
pbivnorm::pbivnorm(x = 1.405072 , y = 0.07526986, rho = 0.3807874, recycle = TRUE)
pbv::pbvnorm(x = 1.405072 , y = 0.07526986, rho = 0.3807874)
pbv::pbvnorm(x = -1.5 , y = Inf, rho = 0.7)




q1 = VGAM::probitlink(1.405072, inverse=T)
q2 = VGAM::probitlink(0.07526986, inverse=T)
VGAM::pbinorm(q2, q1, cov12 = 0.3807874)

