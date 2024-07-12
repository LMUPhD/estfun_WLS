library("psych")
library("resample")
library("mvtnorm")
#####

library(doParallel)
library(foreach)
setwd("C:\\Users\\classe\\Desktop\\Diss\\Paper3\\estfun_WLS")
source("estfun_WLS.R")
source("application\\univ_simu.R")

################################################################################
######################### Ordinal Data ######################################### 
################################################################################





getDataSimulation <- function(items,n,schwellen){
  
  #no fluctuation
  fits_random <- datagen(schwellen = schwellen, ID=n, times=1, items=items)
  simu=fits_random[["data"]][["data1"]]
  simu$rand1 <- replicate(nrow(simu),sample(c(1,2,3,4,5),1,prob=c(0.1,0.2,0.4,0.2,0.1))) #ordinal
  simu$rand2 <- round(runif(nrow(simu),min=100,max=400),2)                               #metric
  simu$rand3 <- replicate(nrow(simu),sample(c(1,2,3,4),1,prob=c(0.15,0.35,0.35,0.15)))   #nominal
  
  
  simu$rand1 <- factor(simu$rand1, ordered =T, levels = c(1,2,3,4,5))
  simu$rand2 <- as.numeric(simu$rand2)
  simu$rand3 <- factor(simu$rand3, ordered=F)
  simu_nonfluc <- simu[sample(nrow(simu)),] ## Shuffle rows
  
  
  
  #fluctuation
  fits_random <- datagen(schwellen = schwellen, ID=n/2, times=2, items=items)
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
  
  
  return( list(simu_nonfluc, simu_fluc, model_lav, items)  )
  #[[1]] nonfluc
  #[[2]] fluc
}

################################################################################
############################### Models ######################################### 
################################################################################

############ Lavaan

getModelsSimulation <- function(x_simus){
  
  models <- list()
  
  for(i in 1:2){
    ############ lavaan
    fit_ord <- lavaan::cfa(x_simus[[3]], data = x_simus[[i]], ordered = TRUE, estimator = "WLS",std.lv=TRUE )
    
    ############ mirt
    fit_mirt <- mirt::mirt(x_simus[[i]][,1:x_simus[[4]]], 1, itemtype="graded",method="EM" )  
    
    
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
    start_time <- Sys.time()
    fluc <- strucchange::gefp(models[[y]][[1]], fit=NULL, scores = mirt::estfun.AllModelClass, order.by = simu[,"rand1"]) #ordinal/categorical
    time <- as.numeric( Sys.time() - start_time )
    flucs1[[1]] <- list(fluc,time)
    
    start_time <- Sys.time()
    fluc <- strucchange::gefp(models[[y]][[1]], fit=NULL, scores = mirt::estfun.AllModelClass, order.by = simu[,"rand2"]) #numeric
    time <- as.numeric( Sys.time() - start_time )
    flucs2[[1]] <- list(fluc,time)
    
    start_time <- Sys.time()
    fluc <- strucchange::gefp(models[[y]][[1]], fit=NULL, scores = mirt::estfun.AllModelClass, order.by = simu[,"rand3"]) #nominal
    time <- as.numeric( Sys.time() - start_time )
    flucs3[[1]] <- list(fluc,time)
    
    #lavaan
    start_time <- Sys.time()
    fluc <- strucchange::gefp(models[[y]][[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand1"]) #ordinal/categorical
    time <- as.numeric( Sys.time() - start_time )
    flucs1[[2]] <- list(fluc,time)
    
    start_time <- Sys.time()
    fluc <- strucchange::gefp(models[[y]][[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand2"]) #numeric
    time <- as.numeric( Sys.time() - start_time )
    flucs2[[2]] <- list(fluc,time)
    
    start_time <- Sys.time()
    fluc <- strucchange::gefp(models[[y]][[2]], fit=NULL, scores = estfun.WLS, order.by = simu[,"rand3"]) #nominal
    time <- as.numeric( Sys.time() - start_time )
    flucs3[[2]] <- list(fluc,time)
    
    
    return(list(flucs1,flucs2,flucs3))
    
  })
  
  return(flucs) 
  #[[1]][[1]][[1]][[1]] nonfluc,ordinal,mirt
  #[[1]][[1]][[2]][[1]] nonfluc,ordinal,WLS
  #[[1]][[2]][[1]][[1]] nonfluc,metric,mirt
  #[[1]][[2]][[2]][[1]] nonfluc,metric,WLS
  #[[1]][[3]][[1]][[1]] nonfluc,nominal,mirt
  #[[1]][[3]][[2]][[1]] nonfluc,nominal,WLS 
  
  #[[2]][[1]][[1]][[1]] fluc,ordinal,mirt
  #[[2]][[1]][[2]][[1]] fluc,ordinal,WLS
  #[[2]][[2]][[1]][[1]] fluc,metric,mirt
  #[[2]][[2]][[2]][[1]] fluc,metric,WLS
  #[[2]][[3]][[1]][[1]] fluc,nominal,mirt
  #[[2]][[3]][[2]][[1]] fluc,nominal,WLS 
}



## Tests
getResultsSimulation <- function(g,.flucs){
  x=g[2]
  z=g[1]
  
  p_maxBB    <- strucchange::sctest(.flucs[[x]][[2]][[z]][[1]], functional = strucchange::maxBB)$p.value #dm
  p_meanL2BB <- strucchange::sctest(.flucs[[x]][[2]][[z]][[1]], functional = strucchange::meanL2BB)$p.value #cvm 
  p_sumlm    <- strucchange::sctest(.flucs[[x]][[2]][[z]][[1]], functional=  strucchange::supLM(from = 0.15, to = NULL))$p.value #suplm
  p_ordL2BB  <- strucchange::sctest(.flucs[[x]][[1]][[z]][[1]], functional = strucchange::ordL2BB(.flucs[[x]][[1]][[z]][[1]]))$p.value[1]  #LM #istead of fluc1 ... factor(simu[,"rand1"])
  p_ordwmax  <- strucchange::sctest(.flucs[[x]][[1]][[z]][[1]], functional = strucchange::ordwmax(.flucs[[x]][[1]][[z]][[1]]))$p.value[1]  #wdmo 
  p_catl2bb  <- strucchange::sctest(.flucs[[x]][[3]][[z]][[1]], functional = strucchange::catL2BB(.flucs[[x]][[3]][[z]][[1]]))$p.value  #LMuo
  
  
  c(p_maxBB,p_meanL2BB,p_sumlm,p_ordL2BB,p_ordwmax,p_catl2bb)
}

## Times fluc
getResultsTime <- function(g,.flucs){
  x=unlist(g[2])
  z=unlist(g[1])
  
  ord_time    <- .flucs[[x]][[1]][[z]][[2]]
  num_time    <- .flucs[[x]][[2]][[z]][[2]]
  cat_time    <- .flucs[[x]][[3]][[z]][[2]]
  
  c(num_time,ord_time,cat_time)
}







################################################################################
for(h in c(3000)){
  for(s in c(1,2,4,6)){
    ##hyperparamters:
    ii=5;nn=h;sc=s
    ##
    
    iters = 1000
    
    cl <- parallel::makeCluster(20)
    doParallel::registerDoParallel(cl)
    result <- foreach(x=1:iters,.packages=c("lavaan","psych","resample"),  .errorhandling="pass") %dopar% {
      
      simus = getDataSimulation(items=ii,n=nn,schwellen=sc)
      models <- getModelsSimulation(simus)
      flucs <- getFlucsSimulation(models,simus)
      grid = expand.grid(1:2,1:2)
      res = apply(grid,1L,function(g) {getResultsSimulation(g,.flucs=flucs)} )
      ###get computation times for flucs...
      res_t = apply(grid,1L,function(g) {getResultsTime(g,.flucs=flucs)} )
      ###...and for models
      ms = apply(grid,1L,function(g) {models[[g[2]]][[g[1]]]} )
      tms = c(ms[[1]]@time[["TOTAL:"]],ms[[2]]@timing$total,ms[[3]]@time[["TOTAL:"]],ms[[4]]@timing$total)
      res_t = rbind(res_t,tms)
      
      
      
      colnames(res) <- colnames(res_t) <- c("nonfluc/mirt","nonfluc/WLS","fluc/mirt","fluc/WLS")
      rownames(res_t) <- c("fluc_num","fluc_ord","fluc_cat","model")
      rownames(res) <- c("maxBB","meanL2BB","sumLM","ordL2BB","ordwmax","catL2BB")
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
    res_arr <- array(as.numeric(unlist(res)), dim=c(6, 4, iters))
    resu <- data.frame()
    for(i in 1:6){
      r1 = sum(  res_arr[i,1,] < .05 )/iters
      r2 = sum(  res_arr[i,2,] < .05 )/iters
      r3 = sum( res_arr[i,3,] > .05 )/iters
      r4 = sum( res_arr[i,4,] > .05 )/iters
      resu <- rbind(resu,c(r1,r2,r3,r4))
      
    }
    
    colnames(resu) <- c("nonfluc/mirt","nonfluc/WLS","fluc/mirt","fluc/WLS")
    rownames(resu) <- c("maxBB","meanL2BB","sumLM","ordL2BB","ordwmax","catL2BB")
    
    ######
    res_t_arr <- array(as.numeric(unlist(res_t)), dim=c(4, 4, iters))
    resu_t <- data.frame()
    
    for(i in 1:4){
      for(j in 1:4){
        r = mean(  res_t_arr[i,j,] )
        resu_t[i,j] <- r
      }
    }
    rownames(resu_t) <- c("fluc_num","fluc_ord","fluc_cat","model")
    colnames(resu_t) <- c("nonfluc/mirt","nonfluc/WLS","fluc/mirt","fluc/WLS")
    ######
    save(result, resu, resu_t,nn,sc, file=paste0("231206_",h,"_",s,"_result_univ_complexWLS.RData"))
    
    
  }
}







