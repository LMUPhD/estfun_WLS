
library(doParallel)
library(foreach)
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


source("estfun_WLS.R")
source("multi_simu.R")


################################################################################
######################### Ordinal Data ######################################### 
################################################################################

getDataSimulation <- function(model,mode,latvars,items,n,schwellen){

  
  
  #fluctuation
  fits_random <- datagen(model = model, mode = mode, schwellen = schwellen, ID=n/2, times=2, items=items, latvar = latvars)
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
  
  
  return( list(simu_fluc, model, items*latvars)  )
  #[[1]] fluc
}


################################################################################
############################### Models ######################################### 
################################################################################

getModelsSimulation <- function(x_simus,model_mirt){
  
  models <- list()
  
  
  ############ mirt
  fit_ord <- lavaan::cfa(x_simus[[2]], data = x_simus[[1]], ordered = TRUE, estimator = "WLS",std.lv=TRUE )
  
  ############ mirt
  fit_mirt <- mirt::mirt(x_simus[[1]][,1:x_simus[[3]]], model=model_mirt, itemtype="graded",method="MHRM" )  
  
  
  models <- list(fit_mirt, fit_ord) 
  
  
  return(models) 
  #[[1]]: fluc mirt 
  #[[2]]: fluc WLS
  }


################################################################################
########################## strucchange #####################################
################################################################################


getFlucsSimulation <- function(models,simus){
  
  flucs1 <- list()
  flucs2 <- list()
  flucs3 <- list()
  simu = simus[[1]]
  
  #mirt
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[1]], fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand1"]) #ordinal/categorical
  time <- as.numeric( Sys.time() - start_time )
  flucs1[[1]] <- list(fluc,time)
  
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[1]], fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand2"]) #numeric
  time <- as.numeric( Sys.time() - start_time )
  flucs2[[1]] <- list(fluc,time)
  
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[1]], fit=NULL, scores = estfun.MHRM, order.by = simu[,"rand3"]) #nominal
  time <- as.numeric( Sys.time() - start_time )
  flucs3[[1]] <- list(fluc,time)
  
  #lavaan
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand1"]) #ordinal/categorical
  time <- as.numeric( Sys.time() - start_time )
  flucs1[[2]] <- list(fluc,time)
  
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand2"]) #numeric
  time <- as.numeric( Sys.time() - start_time )
  flucs2[[2]] <- list(fluc,time)
  
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand3"]) #nominal
  time <- as.numeric( Sys.time() - start_time )
  flucs3[[2]] <- list(fluc,time)
  
  
  return(list(flucs1,flucs2,flucs3))
  
  
  
  #[[1]][[1]][[1]] ordinal,mirt
  #[[1]][[2]][[1]] ordinal,WLS
  #[[2]][[1]][[1]] metric,mirt
  #[[2]][[2]][[1]] metric,WLS
  #[[3]][[1]][[1]] nominal,mirt
  #[[3]][[2]][[1]] nominal,WLS 
}
## Tests
getResultsSimulation <- function(g,.flucs){
  
  p_maxBB    <- strucchange::sctest(.flucs[[2]][[g]][[1]], functional = strucchange::maxBB)$p.value #dm
  p_meanL2BB <- strucchange::sctest(.flucs[[2]][[g]][[1]], functional = strucchange::meanL2BB)$p.value #cvm 
  p_sumlm    <- strucchange::sctest(.flucs[[2]][[g]][[1]], functional=  strucchange::supLM(from = 0.15, to = NULL))$p.value #suplm
  p_ordwmax  <- strucchange::sctest(.flucs[[1]][[g]][[1]], functional = strucchange::ordwmax(.flucs[[1]][[g]][[1]]))$p.value[1]  #wdmo
  p_ordL2BB  <- strucchange::sctest(.flucs[[1]][[g]][[1]], functional = strucchange::ordL2BB(.flucs[[1]][[g]][[1]]))$p.value[1]  #LMo
  p_catl2bb  <- strucchange::sctest(.flucs[[3]][[g]][[1]], functional = strucchange::catL2BB(.flucs[[3]][[g]][[1]]))$p.value  #LMuo
  
  c(p_maxBB,p_meanL2BB,p_sumlm,p_ordL2BB,p_ordwmax,p_catl2bb)
}

## Times fluc
getResultsTime <- function(g,.flucs){

  ord_time    <- .flucs[[1]][[g]][[2]]
  num_time    <- .flucs[[2]][[g]][[2]]
  cat_time    <- .flucs[[3]][[g]][[2]]
  c(num_time,ord_time,cat_time)
}



################################################################################
################################################################################
################################################################################
for(nn in c(500,1000,2000)){
  for(sc in c(1,2,4,6)) {
    
    ### hyperparameters:
    ii=3;lvv=3
    ###
    
    
    model_mirt <- mirt::mirt.model('
  Eta1 = 1-3
  Eta2 = 4-6
  Eta3 = 7-9
  COV=Eta1*Eta2*Eta3') 
    model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'
    
    iters = 1000
    cl <- parallel::makeCluster(64) 
    doParallel::registerDoParallel(cl)
    
    
    result <- foreach(x=1:iters,.packages=c("lavaan","psych","resample"),  .errorhandling="pass") %dopar% {
      
      simus = getDataSimulation(model=model_lav,mode='betas',items=ii,latvars=lvv,n=nn,schwellen=sc)
      models <- getModelsSimulation(simus,model_mirt)
      flucs <- getFlucsSimulation(models,simus)
      res = sapply(1:2,function(g) {getResultsSimulation(g,.flucs=flucs)} )
      ###get computation times for flucs...
      res_t = sapply(1:2,function(g) {getResultsTime(g,.flucs=flucs)} )
      ###...and for models
      ms = lapply(1:2,function(g) {models[[g]]} )
      tms = c(ms[[1]]@time[["TOTAL:"]],ms[[2]]@timing$total)
      res_t = rbind(res_t,tms)
      
      
      
      colnames(res) <- colnames(res_t) <- c("medfluc/mirt","medfluc/WLS")
      rownames(res) <- c("maxBB","meanL2BB","supLM","ordL2BB","ordwmax","catL2BB")
      rownames(res_t) <- c("fluc_num","fluc_ord","fluc_cat","model")
      return(list(res,res_t))
      
    }
    
    parallel::stopCluster(cl)
    
    
    
    ######
    res = lapply(result,"[[",1)
    res_t = lapply(result,"[[",2)
    err = which(sapply(res,is.character))   #exclude errors
    if(length(err)!=0) {
      res = res[-err] 
      res_t = res_t[-err]
    }
    
    
    iters = length(res)
    res_arr <- array(as.numeric(unlist(res)), dim=c(6, 2, iters))
    resu <- data.frame()
    for(i in 1:6){
      r1 = sum(  res_arr[i,1,] > .05 )/iters
      r2 = sum(  res_arr[i,2,] > .05 )/iters
      resu <- rbind(resu,c(r1,r2))
    }
    
    
    colnames(resu) <- c("medfluc/mirt","medfluc/WLS")
    rownames(resu) <- c("maxBB","meanL2BB","supLM","ordL2BB","ordwmax","catL2BB")
    
    
    
    ######
    res_t_arr <- array(as.numeric(unlist(res_t)), dim=c(4, 2, iters))
    resu_t <- data.frame()
    for(i in 1:4){
      for(j in 1:2){
        r = mean(  res_t_arr[i,j,] )
        resu_t[i,j] <- r
      }
    }
    
    rownames(resu_t) <- c("fluc_num","fluc_ord","fluc_cat","model")
    colnames(resu_t) <- c("medfluc/mirt","medfluc/WLS")
    
    ######
    
    save(result, resu, resu_t,ii,lvv,nn,sc, file=paste0("240131_",sc,"_",nn,"_result_multi_medfluc_simpleWLS.RData"))
    
    
  }
}








