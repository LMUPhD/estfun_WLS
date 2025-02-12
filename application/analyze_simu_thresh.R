library("psych")
library("resample")
library("mvtnorm")
#####

library(doParallel)
library(foreach)
setwd("//dss/dsshome1/0B/ra35tik2/paper3")
source("estfun_WLS.R")
source("univ_simu.R")

################################################################################
######################### Ordinal Data ######################################### 
################################################################################

getDataSimulation <- function(items,n,schwellen,mode,dif_items){

  #fluctuation
  fits_random <- datagen(mode = mode, dif_items = dif_items, schwellen = schwellen, ID=n/2, times=2, items=items)
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
  
  
  return( list(simu_fluc, model_lav, items)  )
  #[[1]] fluc
}

################################################################################
############################### Models ######################################### 
################################################################################

############ Lavaan

getModelsSimulation <- function(x_simus){
  
  models <- list()
  
  
  ############ mirt
  fit_ord <- lavaan::cfa(x_simus[[2]], data = x_simus[[1]], ordered = TRUE, estimator = "WLS",std.lv=TRUE )
  
  ############ mirt
  fit_mirt <- mirt::mirt(x_simus[[1]][,1:x_simus[[3]]], 1, itemtype="graded",method="EM" )  
  
  
  models <- list(fit_mirt, fit_ord) 
  
  
  return(models) 
  #[[1]]: fluc mirt 
  #[[2]]: fluc WLS
}


################################################################################
########################## try strucchange #####################################
################################################################################


getFlucsSimulation <- function(models,simus){
  
  flucs1 <- list()
  flucs2 <- list()
  flucs3 <- list()
  simu = simus[[1]]
  
  #mirt
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[1]], fit=NULL, scores = mirt::estfun.AllModelClass, order.by = simu[,"rand1"]) #ordinal/categorical
  time <- as.numeric( Sys.time() - start_time )
  flucs1[[1]] <- list(fluc,time)
  
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[1]], fit=NULL, scores = mirt::estfun.AllModelClass, order.by = simu[,"rand2"]) #numeric
  time <- as.numeric( Sys.time() - start_time )
  flucs2[[1]] <- list(fluc,time)
  
  start_time <- Sys.time()
  fluc <- strucchange::gefp(models[[1]], fit=NULL, scores = mirt::estfun.AllModelClass, order.by = simu[,"rand3"]) #nominal
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
for(nn in c(1000,2000)){  #c(500,1000,2000)
  for(sc in c(1,2,4,6)){  #c(1,2,4,6)
    ##hyperparamters:
    ii=5
    ##
    
    iters = 1000
    
    cl <- parallel::makeCluster(20)
    doParallel::registerDoParallel(cl)
    result <- foreach(x=1:iters,.packages=c("lavaan","psych","resample"),  .errorhandling="pass") %dopar% {
      
      simus = getDataSimulation(items=ii,n=nn,schwellen=sc,mode="thresholds",dif_items=1)
      models <- getModelsSimulation(simus)
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
      r1 = sum( res_arr[i,1,] > .05 )/iters
      r2 = sum( res_arr[i,2,] > .05 )/iters
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
    save(result, resu, resu_t,nn,sc, file=paste0("250201_",nn,"_",sc,"_univ_thresh20.RData"))
    
    
  }
}







