
setwd("C:/Users/classe/Desktop/Diss/Paper3/estfun_WLS")
source("estfun_WLS.R")
source("univ_simu.R")

################################################################################
######################### Ordinal Data ######################################### 
################################################################################





getDataSimulation <- function(items,n,schwellen){
  
  #no fluctuation
  fits_random <- datagen(schwellen = schwellen, ID=n, times=1, items=items)
  simu=fits_random[["data"]][["data1"]][[1]]
  simu$rand1 <- replicate(nrow(simu),sample(c(1,2,3,4,5),1,prob=c(0.1,0.2,0.4,0.2,0.1))) #ordninal
  simu$rand2 <- round(runif(nrow(simu),min=100,max=400),2)                               #numeric (univ)
  
  simu$rand1 <- as.factor(simu$rand1)
  simu$rand2 <- as.numeric(simu$rand2)
  
  simu_nonfluc <- simu[sample(nrow(simu)),] ## Shuffle rows
  
  
  
  #fluctuation
  fits_random <- datagen(schwellen = schwellen, ID=n/2, times=2, items=items) 
  simu1 <- fits_random[["data"]][["data1"]][[1]]
  simu1$rand1 <- replicate(nrow(simu1),sample(c(1,2,3),1,prob=c(0.25,0.5,0.25))) #ordninal
  simu1$rand2 <- round(runif(nrow(simu1),min=100,max=200),2)                               #numeric (univ)
  
  simu2 <- fits_random[["data"]][["data2"]][[1]]
  simu2$rand1 <- replicate(nrow(simu2),sample(c(4,5,6),1,prob=c(0.25,0.5,0.25))) #ordninal
  simu2$rand2 <- round(runif(nrow(simu2),min=200.1,max=400),2)                               #numeric (univ)
  
  simu <- rbind(simu1,simu2)
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
    
    fit_ord <- lavaan::cfa(x_simus[[3]], data = x_simus[[i]], ordered = TRUE, estimator = "WLS",std.lv=TRUE )
    #fit_ord@Fit@npar
    #WLS_scores <- estfun.WLS(fit_ord) #scaling bringt nichts fuer strucchange
    #colSums(WLS_scores)
    
    
    
    ############ mirt
    #model_mirt <- mirt.model('Eta1 = 1-4')
    fit_mirt <- mirt::mirt(x_simus[[i]][,1:x_simus[[4]]], 1, itemtype="graded",method="EM" )  
    #fit_mirt@Model$nest #passt zu oben
    #coef(fit_mirt)
    #mirt_scores <- estfun.AllModelClass(fit_mirt)
    #colSums(mirt_scores)
    
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
  
  flucs <- lapply(1:2,function(x) {
    flucs1 <- list()
    flucs2 <- list()
    
    #mirt
    flucs1[[1]] <- strucchange::gefp(models[[x]][[1]], fit=NULL, scores = estfun.AllModelClass, order.by = simus[[x]][,"rand1"]) #ordinal/categorical
    flucs2[[1]] <- strucchange::gefp(models[[x]][[1]], fit=NULL, scores = estfun.AllModelClass, order.by = simus[[x]][,"rand2"]) #numeric
    
    #lavaan
    flucs1[[2]] <- strucchange::gefp(models[[x]][[2]], fit=NULL, scores = estfun.WLS, order.by = simus[[x]][,"rand1"]) #ordinal/categorical
    flucs2[[2]] <- strucchange::gefp(models[[x]][[2]], fit=NULL, scores = estfun.WLS, order.by = simus[[x]][,"rand2"]) #numeric
    
    return(list(flucs1,flucs2))
    
  })
  
  return(flucs) 
  #[[1]][[1]][[1]] nonfluc,ordinal,mirt
  #[[1]][[1]][[2]] nonfluc,ordinal,WLS
  #[[1]][[2]][[1]] nonfluc,numeric,mirt
  #[[1]][[2]][[2]] nonfluc,numeric,WLS
  #[[2]][[1]][[1]] fluc,ordinal,mirt
  #[[2]][[1]][[2]] fluc,ordinal,WLS
  #[[2]][[2]][[1]] fluc,numeric,mirt
  #[[2]][[2]][[2]] fluc,numeric,WLS
}





## Tests
getResultsSimulation <- function(g,.flucs){
  x=g[2]
  z=g[1]
  
  p_maxBB    <- strucchange::sctest(.flucs[[x]][[2]][[z]], functional = strucchange::maxBB)$p.value #dm
  p_meanL2BB <- strucchange::sctest(.flucs[[x]][[2]][[z]], functional = strucchange::meanL2BB)$p.value #cvm 
  p_sumlm    <- strucchange::sctest(.flucs[[x]][[2]][[z]], functional=  strucchange::supLM(from = 0.15, to = NULL))$p.value #suplm
  p_ordL2BB  <- strucchange::sctest(.flucs[[x]][[1]][[z]], functional = strucchange::ordL2BB(.flucs[[x]][[1]][[z]]))$p.value[1]  #LM #istead of fluc1 ... factor(simu[,"rand1"])
  p_ordwmax  <- strucchange::sctest(.flucs[[x]][[1]][[z]], functional = strucchange::ordwmax(.flucs[[x]][[1]][[z]]))$p.value[1]  #wdmo 
  
  c(p_maxBB,p_meanL2BB,p_sumlm,p_ordL2BB,p_ordwmax)
}




results <- lapply(1:100, function(x,ii=5,nn=1000,sc=4){
  
  simus = getDataSimulation(items=ii,n=nn,schwellen=sc)
  models <- getModelsSimulation(simus)
  flucs <- getFlucsSimulation(models,simus)
  grid = expand.grid(1:2,1:2)
  res = apply(grid,1L,function(g) {getResultsSimulation(g,.flucs=flucs)} )
  colnames(res) <- c("nonfluc/mirt","nonfluc/WLS","fluc/mirt","fluc/WLS")
  rownames(res) <- c("maxBB","meanL2BB","sumLM","ordL2BB","ordwmax")
  return(res)
  
})



res_arr <- array(as.numeric(unlist(results)), dim=c(5, 4, 100))
resu <- data.frame()
for(i in 1:5){
  r1 = sum(  res_arr[i,1,] < .05 )/100
  r2 = sum(  res_arr[i,2,] < .05 )/100
  r3 = sum( res_arr[i,3,] > .05 )/100
  r4 = sum( res_arr[i,4,] > .05 )/100
  resu <- rbind(resu,c(r1,r2,r3,r4))
  
}

colnames(resu) <- c("nonfluc/mirt","nonfluc/WLS","fluc/mirt","fluc/WLS")
rownames(resu) <- c("maxBB","meanL2BB","sumLM","ordL2BB","ordwmax")
################################################################################
############################## Sandbox ######################################### 
################################################################################




fits_random <- datagen(schwellen = 4, ID=1000, times=1, items=5)
simu=fits_random[["data"]][["data1"]][[1]]


### dat with fluctuation
fits_random <- datagen(schwellen = 4, ID=1000, times=2, items=5) 
simu1 <- fits_random[["data"]][["data1"]][[1]]
simu1$rand1 <- replicate(nrow(simu1),sample(c(1,2,3),1,prob=c(0.25,0.5,0.25))) #ordninal
simu1$rand2 <- round(runif(nrow(simu1),min=100,max=200),2)                               #numeric (univ)

simu2 <- fits_random[["data"]][["data2"]][[1]]
simu2$rand1 <- replicate(nrow(simu2),sample(c(4,5,6),1,prob=c(0.25,0.5,0.25))) #ordninal
simu2$rand2 <- round(runif(nrow(simu2),min=200.1,max=400),2)                               #numeric (univ)

simu <- rbind(simu1,simu2)
simu <- simu[sample(nrow(simu)),] ## Shuffle rows


### model
model = fits_random[["model"]][[1]]
fit_ord <- lavaan::cfa(model, data = simu, ordered = TRUE, estimator = "WLS",std.lv=TRUE )
colSums(estfun.WLS(fit_ord))
object=fit_ord


### fluctuation test
flucs1 <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand1"]) #ordinal/categorical
flucs2 <- strucchange::gefp(fit_ord, fit=NULL, scores = estfun.WLS, order.by = simu[,"rand2"]) #numeric

strucchange::sctest(flucs2, functional = strucchange::maxBB)$p.value #dm
strucchange::sctest(flucs2, functional = strucchange::meanL2BB)$p.value #cvm 
strucchange::sctest(flucs2, functional=  strucchange::supLM(from = 0.15, to = NULL))$p.value #suplm
strucchange::sctest(flucs1, functional = strucchange::ordL2BB(flucs1))$p.value[1]  #LM #istead of fluc1 ... factor(simu[,"rand1"])
strucchange::sctest(flucs1, functional = strucchange::ordwmax(flucs1))$p.value[1]  #wdmo 

