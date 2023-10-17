
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")

estfun.MHRM <- function(modelobj){
  model.mhrm.values <- mirt::mirt(data=simu[,1:modelobj@Data[["nitems"]]],model=model_mirt, itemtype='graded', pars = "values", method='MHRM') 
  
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
  model.mhrm.EM <- mirt::mirt(data=simu[,1:modelobj@Data[["nitems"]]],model=model_mirt, itemtype='graded', method='EM', pars = model.mhrm.values, TOL = NaN)
  
  # Again, we get the score contributions:
  return(mirt::estfun.AllModelClass(model.mhrm.EM))
  
}


source("estfun_WLS2.R")
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
start <- Sys.time()
fit_ord <- lavaan::cfa(model_lav, data = simu, ordered = TRUE, estimator = "WLS",std.lv=T ) #estimate variances and covariances of latent variables
t_lav = as.numeric(Sys.time() - start)

fitMeasures(fit_ord)[['rmsea']]
fit_ord@Model@nx.free


#mirt
model_mirt <- mirt::mirt.model('
Eta1 = 1-3
Eta2 = 4-6
Eta3 = 7-9
COV=Eta1*Eta2*Eta3') 
fit_mirt <- mirt::mirt(data=simu, model=model_mirt, itemtype="graded",method='MHRM' ) #pars="values" #'EM',TOL = NaN 
t_mirt = fit_mirt@time[["TOTAL:"]]

mirt::extract.mirt(fit_mirt, what = "RMSEA")
fit_mirt@Model[["nest"]]



################################################################################
########################## try strucchange #####################################
################################################################################
simu$rand1 <- replicate(nrow(simu),sample(c(1,2,3,4,5),1,prob=c(0.1,0.2,0.4,0.2,0.1))) #ordninal
simu$rand2 <- round(runif(nrow(simu),min=0,max=200),2)                               #numeric (univ)


flucs1 <- strucchange::gefp(fit_mirt, fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand1"]) #ordinal
flucs2 <- strucchange::gefp(fit_mirt, fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand2"]) #numeric

strucchange::sctest(flucs2, functional = strucchange::maxBB)$p.value #dm    #--> passt... andere stats funktionieren mit >25 parametern nicht
strucchange::sctest(flucs1, functional = strucchange::ordwmax(flucs1))$p.value[1]  #wdmo

















################################################################################
############################## try semtree #####################################
################################################################################
fit_num <- lavaan::cfa(model = model_lav, data=simu, std.lv=T, meanstructure=T,estimator = "ML",do.fit=FALSE) #, do.fit=FALSE
semtree_simu <- semtree::semtree(model=fit_ord, data=simu, predictors=c("rand1","rand2"), control=semtree::semtree.control(method="score",bonferroni=TRUE)) 

#plot.semtree aus github kopieren...
nodeFunSemtree<-function(x, labs, digits, varlen)
{
  paste("n=",x$frame$n)
}
plot.semtree(semtree_simu)


####Fazit: SEMTree geht noch nicht!

