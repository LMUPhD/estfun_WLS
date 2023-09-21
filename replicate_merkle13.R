
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")
#load("liss_univ_data.RData") 
source("estfun_WLS.R")

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
WLS_scores <- estfun.WLS(object)


#Daten besorgen mit drei Antwortkategorien!


################################################################################
############################## Sandbox ######################################### 
################################################################################



p1_1 = th[1]*-1
p1_0 = th[1]
p2_1 = th[3]*-1
p2_0 = th[3]

pbivnorm::pbivnorm(x = p1_0, y = p2_0, rho = 0.2906113, recycle = TRUE) 
pbivnorm::pbivnorm(x = p1_1, y = p2_0, rho = 0.2906113, recycle = TRUE) 
pbivnorm::pbivnorm(x = p1_0, y = p2_1, rho = 0.2906113, recycle = TRUE) 
pbivnorm::pbivnorm(x = p1_1, y = p2_1, rho = 0.2906113, recycle = TRUE) 
#doesnt add up to one... 


th1 = th[1] 
th2 = th[3]
#P(y1=0, y2=0) = 0.081
sum(c(
  pbivnorm_wls(x =   th1, y = th2, rho = 0.2906113),
  pbivnorm_wls(x =  -Inf, y = th2, rho = 0.2906113)*-1,
  pbivnorm_wls(x =  th1, y = -Inf, rho = 0.2906113)*-1,
  pbivnorm_wls(x = -Inf, y = -Inf, rho = 0.2906113)
))



#P(y1=1, y2=0) = 0.3034526 ???WTF
sum(c(
  pbivnorm_wls(x = Inf,  y = th2,  rho = 0.2906113),
  pbivnorm_wls(x = th1,  y = th2,  rho = 0.2906113)*-1,
  pbivnorm_wls(x = Inf,  y = -Inf, rho = 0.2906113)*-1,
  pbivnorm_wls(x = th1,  y = -Inf, rho = 0.2906113)
))




#P(y1=0, y2=1) = 0.1334526 ???WTF
sum(c(
  pbivnorm_wls(x = th1,  y = Inf,  rho = 0.2906113),
  pbivnorm_wls(x = -Inf, y = Inf,  rho = 0.2906113)*-1,
  pbivnorm_wls(x = th1,  y = th2,  rho = 0.2906113)*-1,
  pbivnorm_wls(x = -Inf, y =  th2, rho = 0.2906113)
))


#P(y1=1, y2=1) = 0.567
sum(c(
  pbivnorm_wls(x = Inf, y = Inf, rho = 0.2906113),
  pbivnorm_wls(x = th1, y = Inf, rho = 0.2906113)*-1,
  pbivnorm_wls(x = Inf, y = th2, rho = 0.2906113)*-1,
  pbivnorm_wls(x = th1, y = th2, rho = 0.2906113)
))








sum(0.081,0.091,0.261,0.567) #muesste eigentlich passen...

###########




