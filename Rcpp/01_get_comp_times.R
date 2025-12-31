setwd("C:\\Users\\classe\\Desktop\\Diss\\Paper3\\estfun_WLS")
source("application/multi_simu.R") #simulate data (multivariate model)


################################################################################
########################## GEE estimation non-binary ###########################
################################################################################

model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'
model = model_lav
fits_random <- datagen(model=model_lav, schwellen = 4, ID=2000, times=1, items=3,latvar=3) 
Data = fits_random[["data"]][["data1"]]

#estimate parmeters
fit.wls <- lavaan::cfa(model_lav, data = Data, ordered = TRUE, estimator = "WLS", std.lv=F )

#apply gee  --> estimate time...
source("estfun_gee.R")
system.time(
  estfun.GEE(fit.wls)
) #8.7 sec
scores_r = estfun.GEE(fit.wls)

source("Rcpp/estfun_gee_rcpp.R")
system.time(
  estfun.GEE(fit.wls)
) #0.89!!! das ist ca. ein zehntel!!!
scores_rcpp = estfun.GEE(fit.wls)

#compare
any(  abs(scores_rcpp - scores_r) > 1e+9 ) #very accurate!!!

######   microbenchmarking estfun.gee
library(profvis)
profvis({
  estfun.GEE(fit.wls)
})
