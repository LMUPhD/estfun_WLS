
setwd("C:/Users/classe/Desktop/Diss/Paper3/R")
load("liss_univ_data.RData") 


## mirt:
library(mirt)
dt <- mirt::expand.table(LSAT7)
fit_mirt <- mirt(dt[,1:5], 1, itemtype = "Rasch", SE = TRUE)
mirt_scores <- estfun.AllModelClass(fit_mirt)

####################
library(lavaan)

#model ='Eta1 =~ cp08a014 + cp08a015 + cp08a016 + cp08a017 + cp08a018'
#fit <- cfa(model, data = dt, meanstructure = TRUE)   #fit with ML
#sc <- estfun.lavaan(fit)                      #calculate scores for model fit with ML



model = 'Eta1 =~ 1*Item.1 + 1*Item.2 + 1*Item.3 + 1*Item.4 + 1*Item.5'

fit_ord <- cfa(model, data = dt, ordered = TRUE, estimator = "WLS" ) #parameterization = "theta"
object = fit_ord



#Daten besorgen mit drei Antwortkategorien!


################################################################################
############################## Sandbox ######################################### 
################################################################################









