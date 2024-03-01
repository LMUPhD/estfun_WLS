library(partykit)
library(mirt)
library(strucchange)
library(lavaan)

setwd("C:\\Users\\classe\\Desktop\\Diss\\Paper3\\estfun_WLS")
source("multi_simu.R")
source("estfun_WLS.R")

#mob not possible with WLS because no log-likelihood is computed!

#############################################################################
#######################   Functions   ##########################
#############################################################################
create_dataset <- function(fits_random){
  
  goodfits1 <- fits_random[["data"]][[1]] 
  goodfits2 <- fits_random[["data"]][[2]]
  goodfits3 <- fits_random[["data"]][[3]] 
  goodfits4 <- fits_random[["data"]][[4]]
  
  
  ID=nrow(goodfits1)
  
  #sample variable
  goodfits1$sample <- rep(1,ID)
  goodfits2$sample <- rep(2,ID)
  goodfits3$sample <- rep(3,ID)
  goodfits4$sample <- rep(4,ID)
  
  
  goodfits1$cat1 <- replicate(ID,sample(c(1,5),1))
  goodfits2$cat1 <- replicate(ID,sample(c(2,3,4),1))
  goodfits3$cat1 <- replicate(ID,sample(c(1,2,3,4,5),1))
  goodfits4$cat1 <- replicate(ID,sample(c(1,2,3,4,5),1))
  
  
  goodfits1$num1 <- round(runif(ID,min=1,max=100),2) 
  goodfits2$num1 <- round(runif(ID,min=1,max=100),2)
  goodfits3$num1 <- round(runif(ID,min=101,max=200),2)
  goodfits4$num1 <- round(runif(ID,min=101,max=200),2)
  
  goodfits1$ord1 <- replicate(ID,sample(c(1,2,3,4),1))
  goodfits2$ord1 <- replicate(ID,sample(c(1,2,3,4),1))
  goodfits3$ord1 <- replicate(ID,sample(c(1,2),1))
  goodfits4$ord1 <- replicate(ID,sample(c(3,4),1))
  
  
  
  ###### Whole Dataset
  simu <- rbind(goodfits1,goodfits2,goodfits3,goodfits4)
  simu <- simu[sample(nrow(simu)),] ## Shuffle rows
  
  ##### Add random partitioning variables
  simu$rand1 <-  replicate(nrow(simu),sample(c(1,2,3,4,5),1,prob=c(0.1,0.2,0.4,0.2,0.1))) #ordinal
  simu$rand2 <-  round(runif(nrow(simu),min=1,max=5),0)                                   #ordinal
  simu$rand3 <-  replicate(nrow(simu),sample(c(0,1),1,prob=c(0.7,0.3)))                   #dichotomous
  simu$rand4 <-  round(runif(nrow(simu),min=100,max=400),2)                               #numeric (univ)
  simu$rand5 <-  round(rnorm(nrow(simu),mean=30,sd=5),4)                                  #numeric (norm)
  
  
  
  ##### Prepocessing 
  simu$ord1 <- as.ordered(simu$ord1) 
  simu$rand1 <- as.ordered(simu$rand1) 
  simu$rand2 <- as.ordered(simu$rand2) 
  
  simu$cat1 <- as.factor(simu$cat1) 
  simu$rand3 <- as.factor(simu$rand3) 
  
  simu$num1 <- as.numeric(simu$num1) 
  simu$rand5 <- as.numeric(simu$rand5) 
  simu$rand4 <- as.numeric(simu$rand4)
  
  return(simu)
}

estfun.MHRM <- function(modelobj){
  model.mhrm.values <- mirt::mirt(data=modelobj@Data[["data"]],model=model_mirt, itemtype='graded', pars = "values", method='MHRM') 
  ## We extract the estimated parameters:
  params <- coef(modelobj, simplify=TRUE)
  # 1. The rowwise entries of $items
  v1 <- as.vector(t(params$items))
  # 2. The entries of $means
  v2 <- as.vector(params$means)
  # 3. The lower diagonal of $cov
  v3 <- params$cov[lower.tri(params$cov, diag = T)]
  
  # We concatenate these vectors and insert them into model.recoded.values:
  model.mhrm.values$value <- c(v1, v2, v3)
  
  # We now evaluate the model with EM at the estimates.
  model.mhrm.EM <- mirt::mirt(data=modelobj@Data[["data"]],model=model_mirt, itemtype='graded', method='EM', pars = model.mhrm.values, TOL = NaN)
  
  # Again, we get the score contributions:
  return(mirt::estfun.AllModelClass(model.mhrm.EM))
  
}

mirt_fit <- function(model) {
  
  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE) {
    
    fit_mirt <- mirt(data = y[,1:9], model = model, itemtype="graded",method="MHRM")
    list(
      estfun = estfun.MHRM(fit_mirt), #Schritt 1: Variable selection
      objfun = fit_mirt@Fit[["logLik"]], #Schritt 2: Split point selection
      coefficients = stats4::coef(fit_mirt),#?ndert nichts an Tree-Growind --> f?r kompakte Darstellung im plot --> nur Varianzen der latenten Variablen
      object = if(object) fit_mirt else NULL #lavaan Objekt wird ausgegeben --> man kann sich alle fitmeasures ausgeben!!
    )
  }
}

lav_fit <- function(model) {
  
  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE) {
    
    fit_ord <- lavaan::cfa(model = model, data=y, std.lv=F, ordered=T, estimator = "WLS", control=list(iter.max=100)) 
    fit_num <- lavaan::cfa(model = model, data=y, std.lv=F, estimator = "ML", control=list(iter.max=100)) 
    
    print("iteration")
    
    list(
      estfun = estfun.WLS(fit_ord), #Schritt 1: Variable selection
      objfun = -as.numeric(fit_num@loglik$loglik), #Schritt 2: Split point selection
      coefficients = stats4::coef(fit_ord),
      object = if(object) fit_ord else NULL #lavaan Objekt wird ausgegeben
    )
  }
}


#############################################################################
#######################   Dataset   ##########################
#############################################################################

model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'
fits_random <- datagen(model = model_lav, mode='all', schwellen = 1, ID=500, times=4, items=3, latvar = 3)
simu = create_dataset(fits_random)






#############################################################################
##############################     MOB     ##################################
#############################################################################

#MIRT
start_time <- Sys.time()
mob_mirt <- mob(formula=" simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5 + simuvar6 + simuvar7 + simuvar8 + simuvar9 ~ cat1 + num1 + ord1 + rand1 + rand2 + rand3 + rand4 + rand5",
                data = simu,
                fit = mirt_fit(model_mirt),
                control = mob_control(ytype = "data.frame",bonferroni = T, maxdepth = Inf,minsize=100)) 
end_time <- Sys.time()
runtime_mob <- end_time - start_time #takes forever!
#save.image("mob_mirt_result.RData")





#lavaan
start_time <- Sys.time()
mob_lav <- mob(formula=" simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5 + simuvar6 + simuvar7 + simuvar8 + simuvar9 ~ cat1 + num1 + ord1 + rand1 + rand2 + rand3 + rand4 + rand5",
                data = simu,
                fit = lav_fit(model_lav),
                control = mob_control(ytype = "data.frame",bonferroni = T, maxdepth = Inf,minsize=100)) 
end_time <- Sys.time()
runtime_mob <- end_time - start_time #takes forever!
#save.image("mob_mirt_result.RData")



nodefun <-function(i) c(
  paste("n =", i$n)
)
plot(mob_lav,terminal_panel = node_terminal,tp_args = list(FUN = nodefun)) #YES!























#############################################################################
#############################     Mirt     ##################################
#############################################################################

model_mirt <- mirt::mirt.model('
  Eta1 = 1-3
  Eta2 = 4-6
  Eta3 = 7-9
  COV=Eta1*Eta2*Eta3') 


fit_mirt <- mirt(data = simu[,1:9], model = model_mirt, itemtype="graded",method="MHRM")
scores <- estfun.MHRM(fit_mirt)

stab_mirt <- gefp(fit_mirt, fit=NULL, scores=estfun.MHRM, order.by = simu[,"num1"])
sctest(stab_mirt,functional=  supLM(from = 0.15, to = NULL)) #suplm


#############################################################################
#############################     lavaan    ##################################
#############################################################################

#ord
fit_ord <- lavaan::cfa(model = model_lav, data=simu, std.lv=F, ordered=T, estimator = "WLS") 
scores <- estfun.WLS(fit_ord)

stab_wls <- gefp(fit_ord, fit=NULL, scores=estfun.WLS, order.by = simu[,"num1"])
sctest(stab_wls,functional=  supLM(from = 0.15, to = NULL)) #suplm
stab_wls <- gefp(fit_ord, fit=NULL, scores=estfun.WLS, order.by = simu[,"cat1"])
sctest(stab_wls,functional = strucchange::catL2BB(stab_wls)) #catL2BB
stab_wls <- gefp(fit_ord, fit=NULL, scores=estfun.WLS, order.by = simu[,"ord1"])
sctest(stab_wls, functional = strucchange::ordL2BB(stab_wls))  #LMo


#num --> funktioniert nicht
fit_num <- lavaan::cfa(model = model, data=simu, std.lv=F,meanstructure=T, estimator = "ML", control=list(iter.max=10000))

stab_num <- gefp(fit_num, fit=NULL, scores=lavScores, order.by = simu[,"num1"])
sctest(stab_num,functional=  supLM(from = 0.15, to = NULL)) #suplm
stab_num <- gefp(fit_num, fit=NULL, scores=lavScores, order.by = simu[,"cat1"])
sctest(stab_num,functional = strucchange::catL2BB(stab_wls)) #catL2BB
stab_num <- gefp(fit_num, fit=NULL, scores=lavScores, order.by = simu[,"ord1"])
sctest(stab_num, functional = strucchange::ordL2BB(stab_wls))  #LMo

