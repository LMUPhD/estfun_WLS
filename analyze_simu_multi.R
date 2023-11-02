
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")

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


source("estfun_WLS2.R")
source("multi_simu.R")


################################################################################
######################### Ordinal Data ######################################### 
################################################################################






getDataSimulation <- function(model,latvars,items,n,schwellen){
  
  #no fluctuation
  fits_random <- datagen(model = model, schwellen = schwellen, ID=n, times=1, items=items, latvar = latvars)
  simu=fits_random[["data"]][["data1"]]
  simu$rand1 <- replicate(nrow(simu),sample(c(1,2,3,4,5),1,prob=c(0.1,0.2,0.4,0.2,0.1))) #ordinal
  simu$rand2 <- round(runif(nrow(simu),min=100,max=400),2)                               #metric
  simu$rand3 <- replicate(nrow(simu),sample(c(1,2,3,4),1,prob=c(0.15,0.35,0.35,0.15)))   #nominal
  
  
  simu$rand1 <- factor(simu$rand1, ordered =T, levels = c(1,2,3,4,5))
  simu$rand2 <- as.numeric(simu$rand2)
  simu$rand3 <- factor(simu$rand3, ordered=F)
  simu_nonfluc <- simu[sample(nrow(simu)),] ## Shuffle rows
  
  
  
  #fluctuation
  fits_random <- datagen(model = model, schwellen = schwellen, ID=n/2, times=2, items=items, latvar = latvars)
  simu1 <- fits_random[["data"]][["data1"]]
  simu1$rand1 <- replicate(nrow(simu1),sample(c(1,2,3),1,prob=c(0.25,0.5,0.25))) #ordinal
  simu1$rand2 <- round(runif(nrow(simu1),min=100,max=200),2)                     #metric
  simu1$rand3 <- replicate(nrow(simu1),sample(c(1,4),1,prob=c(0.5,0.5)))         #ordinal
  
  simu2 <- fits_random[["data"]][["data2"]]
  simu2$rand1 <- replicate(nrow(simu2),sample(c(4,5,6),1,prob=c(0.25,0.5,0.25))) #ordinal
  simu2$rand2 <- round(runif(nrow(simu2),min=200.1,max=400),2)                   #metric
  simu2$rand3 <- replicate(nrow(simu2),sample(c(2,3),1,prob=c(0.5,0.5)))         #nominal
  
  simu <- rbind(simu1,simu2)
  simu$rand1 <- factor(simu$rand1, ordered =T, levels = c(1,2,3,4,5,6))
  simu$rand2 <- as.numeric(simu$rand2)
  simu$rand3 <- factor(simu$rand3, ordered=F)
  simu_fluc <- simu[sample(nrow(simu)),] ## Shuffle rows
  
  #also export lavaan model
  model_lav = fits_random[["model"]][[1]]
  
  
  return( list(simu_nonfluc, simu_fluc, model, items*latvars)  )
  #[[1]] nonfluc
  #[[2]] fluc
}


################################################################################
############################### Models ######################################### 
################################################################################


getModelsSimulation <- function(x_simus,model_mirt){
  
  models <- list()
  
  for(i in 1:2){
    ############ mirt
    fit_ord <- lavaan::cfa(x_simus[[3]], data = x_simus[[i]], ordered = TRUE, estimator = "WLS",std.lv=TRUE )

    ############ mirt
    fit_mirt <- mirt::mirt(x_simus[[i]][,1:x_simus[[4]]], model=model_mirt, itemtype="graded",method="MHRM" )  

    
    models[[i]] <- list(fit_mirt, fit_ord) 
    
  }
  return(models) 
  #[[1]][[1]]: nonfluc mirt 
  #[[1]][[2]]: nonfluc WLS 
  #[[2]][[1]]: fluc mirt 
  #[[2]][[2]]: fluc WLS

}




################################################################################
########################## try strucchange #####################################
################################################################################



##################### Fluctuation test

#x=1 --> nonfluc 
#x=2 --> fluc

getFlucsSimulation <- function(models,simus){
  
  flucs <- lapply(1:2,function(y) {
    flucs1 <- list()
    flucs2 <- list()
    flucs3 <- list()
    simu = simus[[y]]
    
    #mirt
    flucs1[[1]] <- strucchange::gefp(models[[y]][[1]], fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand1"]) #ordinal/categorical
    flucs2[[1]] <- strucchange::gefp(models[[y]][[1]], fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand2"]) #numeric
    flucs3[[1]] <- strucchange::gefp(models[[y]][[1]], fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand3"]) #nominal
    
    
    #lavaan
    flucs1[[2]] <- strucchange::gefp(models[[y]][[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand1"]) #ordinal/categorical
    flucs2[[2]] <- strucchange::gefp(models[[y]][[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand2"]) #numeric
    flucs3[[2]] <- strucchange::gefp(models[[y]][[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand3"]) #nominal
    
    
    return(list(flucs1,flucs2,flucs3))
    
  })
  
  return(flucs) 
  #[[1]][[1]][[1]] nonfluc,ordinal,mirt
  #[[1]][[1]][[2]] nonfluc,ordinal,WLS
  #[[1]][[2]][[1]] nonfluc,metric,mirt
  #[[1]][[2]][[2]] nonfluc,metric,WLS
  #[[1]][[3]][[1]] nonfluc,nominal,mirt
  #[[1]][[3]][[2]] nonfluc,nominal,WLS 
  
  #[[2]][[1]][[1]] fluc,ordinal,mirt
  #[[2]][[1]][[2]] fluc,ordinal,WLS
  #[[2]][[2]][[1]] fluc,metric,mirt
  #[[2]][[2]][[2]] fluc,metric,WLS
  #[[2]][[3]][[1]] fluc,nominal,mirt
  #[[2]][[3]][[2]] fluc,nominal,WLS 
}





## Tests
getResultsSimulation <- function(g,.flucs){
  x=unlist(g[2])
  z=unlist(g[1])
  
  p_maxBB    <- strucchange::sctest(.flucs[[x]][[2]][[z]], functional = strucchange::maxBB)$p.value #dm
  p_ordwmax  <- strucchange::sctest(.flucs[[x]][[1]][[z]], functional = strucchange::ordwmax(.flucs[[x]][[1]][[z]]))$p.value[1]  #wdmo 
  p_catl2bb  <- strucchange::sctest(.flucs[[x]][[3]][[z]], functional = strucchange::catL2BB(.flucs[[x]][[3]][[z]]))$p.value  #LMuo
  
  c(p_maxBB,p_ordwmax,p_catl2bb)
}

















################################################################################

model_mirt <- mirt::mirt.model('
Eta1 = 1-3
Eta2 = 4-6
Eta3 = 7-9
COV=Eta1*Eta2*Eta3') 
model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'
iters = 10
for(i in 1:iters){
  ii=3;lvv=3;nn=1000;sc=4
  simus = getDataSimulation(model=model_lav,items=ii,latvars=lvv,n=nn,schwellen=sc)
  models <- getModelsSimulation(simus,model_mirt)
  flucs <- getFlucsSimulation(models,simus)
  grid = expand.grid(1:2,1:2)
  res = apply(grid,1L,function(g) {getResultsSimulation(g,.flucs=flucs)} )
  ### get computation times
  ms = apply(grid,1L,function(g) {models[[g[2]]][[g[1]]]} )
  tms = c(ms[[1]]@time[["TOTAL:"]],ms[[2]]@timing$total,ms[[3]]@time[["TOTAL:"]],ms[[4]]@timing$total)
  res = rbind(res,tms)
  ###
  
  colnames(res) <- c("nonfluc/mirt","nonfluc/WLS","fluc/mirt","fluc/WLS")
  rownames(res) <- c("maxBB","ordmax","catL2BB","times")
  print(res)
}





