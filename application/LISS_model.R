library(haven)
library(tidyverse)
devtools::install_github("yrosseel/lavaan")
library(lavaan)


setwd("/dss/dsshome1/0B/ra35tik2/paper2/real_data/rawdata")

backvars <- read_sav("avars_200808_EN_2.0p.sav")
dt08a = read_sav("cp08a_1p_EN.sav")
dt09b = read_sav("cp09b_1.0p_EN.sav")
#dt10c = read_sav("cp10c_1.0p_EN.sav")
dt11d = read_sav("cp11d_1.0p_EN.sav")
#dt12e = read_sav("cp12e_1.0p_EN.sav")
dt13f = read_sav("cp13f_EN_1.0p.sav")
#dt14g = read_sav("cp14g_EN_1.0p.sav")



################################################################################
###################  Datensatz mehrere Wellen erstellen ########################
################################################################################

ids <-  Reduce(intersect,list(backvars$nomem_encr,dt08a$nomem_encr,dt09b$nomem_encr,dt11d$nomem_encr,dt13f$nomem_encr))   
print(length(ids)) #2893

#SL --> happiness
ye = c("08a","09b","11d","13f")
it = paste0("0",14:18) #Satisfaction with life
#it_se = paste0("0",70:79) #Self esteem (Rosenberg)
#it = c(it_ls,it_se)

data  <-  vars <- var <- va <- list()
for (i in 1:length(ye)){
  for (j in 1:length(it)){
    var[j] <-  paste0("cp",ye[i],it[j])  
  }
  vars[[i]] <- var 
}


#multi data
data[[1]] <- backvars 
data[[2]] <- dt08a[, c("nomem_encr",unlist(vars[[1]]))] 
data[[3]] <- dt09b[, c("nomem_encr",unlist(vars[[2]]))]  
data[[4]] <- dt11d[, c("nomem_encr",unlist(vars[[3]]))] 
data[[5]] <- dt13f[, c("nomem_encr",unlist(vars[[4]]))] 


dt = data %>% reduce(inner_join, by = "nomem_encr")
dt = dt %>% drop_na(unlist(vars)) %>% as.data.frame() %>% mutate_at(unlist(vars), as.numeric)

#backvars formatieren
for(i in c("aantalhh","aantalki","lftdcat","oplmet","nettocat","oplzon","sted")){ #ordered factor
  dt[,i] <- as.character(dt[,i])
  dt[is.na(dt[,i]),i] <- "-99"
  dt[,i] <- factor(dt[,i], ordered = TRUE)
}

for(i in c("geslacht","partner","burgstat","woning","woonvorm","belbezig")){ #factor 
  dt[,i] <- as.character(dt[,i])
  dt[is.na(dt[,i]),i] <- "-99"
  dt[,i] <- as.factor(dt[,i])
}

for(i in c("brutoink","nettoink","leeftijd")){ #numeric
  dt[,i] <- as.character(dt[,i])
  dt[is.na(dt[,i]),i] <- "-99"
  dt[,i] <- as.numeric(dt[,i])
}


##################### Listen 
models <- mfit <- mfit_time <- flucs <- flucs_time <- mirtmodels <- pvals <- list()

################################################################################
#######################  Sehr Kleines model LS  ################################
################################################################################



model_multi <- "
etaLS1=~ cp08a014 + cp08a015 + cp08a016 + cp08a017+ cp08a018
"
fit_num <- lavaan::cfa(model = model_multi, data=dt, estimator="MLR", std.lv=F) 
fit_ord <- lavaan::cfa(model = model_multi, data=dt, ordered = TRUE, estimator="WLS", std.lv=F) 
#summary(fit_ord, fit.measures=T) #35 params #RMSEA=0.122
#summary(fit_num, fit.measures=T) #10 params #RMSEA=0.077 (scaled)


dt_mirt=dt[,unlist(vars[[1]])]
fit_mirt <- mirt::mirt(dt_mirt, 1, itemtype="graded",method="EM" )  
#coef(fit_mirt, simplify=T)


###write
mfit[[1]] <- fitMeasures(fit_num)["rmsea.robust"]
mfit[[2]] <- fitMeasures(fit_ord)["rmsea"]
mfit[[3]] <- mirt::extract.mirt(fit_mirt, "RMSEA")

mfit_time[[1]] <- fit_num@timing["total"][[1]]
mfit_time[[2]] <- fit_ord@timing["total"][[1]]
mfit_time[[3]] <- fit_mirt@time[["TOTAL:"]]

models[[1]] <- fit_num
models[[2]] <- fit_ord
models[[3]] <- fit_mirt
names(mfit) <- names(mfit_time) <- names(models) <- c("MLR","WLS","MML")

mfit_m1 = mfit
mfit_time_m1 = mfit_time
models_m1 = models
################################################################################
##########################  Kleines model LS  ##################################
################################################################################


model_multi <- "
etaLS1=~ 1*cp08a014 + 1*cp08a015 + 1*cp08a016 + 1*cp08a017+ 1*cp08a018
etaLS2=~ 1*cp09b014 + 1*cp09b015 + 1*cp09b016 + 1*cp09b017+ 1*cp09b018
"
fit_num <- lavaan::cfa(model = model_multi, data=dt, estimator="MLR", std.lv=F) 
fit_ord <- lavaan::cfa(model = model_multi, data=dt, ordered = TRUE, estimator="WLS", std.lv=F) 
#summary(fit_ord, fit.measures=T) #63 params #RMSEA=0.127
#summary(fit_num, fit.measures=T) #13 params #RMSEA=0.100 (scaled)


dt_mirt=dt[,unlist(vars[1:2])]
model_mirt_m2 <- mirt::mirt.model('
  Eta1 = 1-5
  Eta2 = 6-10
  
  CONSTRAIN = (1, a1, 1), (2, a1, 1), (3, a1, 1), (4, a1, 1),(5, a1, 1),
  (6, a2, 1), (7, a2, 1), (8, a2, 1), (9, a2, 1), (10, a2, 1),

  COV=Eta1*Eta2') 
fit_mirt <- mirt::mirt(dt_mirt, model=model_mirt_m2, itemtype="graded",method="MHRM" )  
#coef(fit_mirt, simplify=T)


###write
mfit[[1]] <- fitMeasures(fit_num)["rmsea.robust"]
mfit[[2]] <- fitMeasures(fit_ord)["rmsea"]
mfit[[3]] <- mirt::extract.mirt(fit_mirt, "RMSEA")

mfit_time[[1]] <- fit_num@timing["total"][[1]]
mfit_time[[2]] <- fit_ord@timing["total"][[1]]
mfit_time[[3]] <- fit_mirt@time[["TOTAL:"]]

models[[1]] <- fit_num
models[[2]] <- fit_ord
models[[3]] <- fit_mirt
names(mfit) <- names(mfit_time) <- names(models) <- c("MLR","WLS","MML")

mfit_m2 = mfit
mfit_time_m2 = mfit_time
models_m2 = models



################################################################################
#############################  PIEG model LS  ##################################
################################################################################


model_multi <- "
etaLS1=~ 1*cp08a014 + 1*cp08a015 + 1*cp08a016 + 1*cp08a017+ 1*cp08a018
etaLS2=~ 1*cp09b014 + 1*cp09b015 + 1*cp09b016 + 1*cp09b017+ 1*cp09b018
etaLS3=~ 1*cp11d014 + 1*cp11d015 + 1*cp11d016 + 1*cp11d017+ 1*cp11d018
etaLS4=~ 1*cp13f014 + 1*cp13f015 + 1*cp13f016 + 1*cp13f017+ 1*cp13f018
deltaLS015=~1*cp08a015+1*cp09b015+1*cp11d015+1*cp13f015
deltaLS016=~1*cp08a016+1*cp09b016+1*cp11d016+1*cp13f016
deltaLS017=~1*cp08a017+1*cp09b017+1*cp11d017+1*cp13f017
deltaLS018=~1*cp08a018+1*cp09b018+1*cp11d018+1*cp13f018
"
fit_num <- lavaan::cfa(model = model_multi, data=dt, estimator="MLR", std.lv=F) 
fit_ord <- lavaan::cfa(model = model_multi, data=dt, ordered = TRUE, estimator="WLS", std.lv=F) 
#summary(fit_ord, fit.measures=T) #156 params #RMSEA=0.051
#summary(fit_num, fit.measures=T) #65 params #RMSEA=0.036 (scaled)



dt_mirt=dt[,unlist(vars)]
model_mirt_m3 <- mirt::mirt.model('
  Eta1 = 1-5
  Eta2 = 6-10
  Eta3 = 11-15
  Eta4 = 16-20
  Beta2 = 2,7,12,17
  Beta3 = 3,8,13,18
  Beta4 = 4,9,14,19
  Beta5 = 5,10,15,20
  
  CONSTRAIN = (1, a1, 1), (2, a1, 1), (3, a1, 1), (4, a1, 1),(5, a1, 1),
  (6, a2, 1), (7, a2, 1), (8, a2, 1), (9, a2, 1), (10, a2, 1),
  (11, a3, 1), (12, a3, 1), (13, a3, 1), (14, a3, 1),(15, a3, 1),
  (16, a4, 1), (17, a4, 1), (18, a4, 1), (19, a4, 1),(20, a4, 1),
  (2, a5, 1), (7, a5, 1), (12, a5, 1),(17, a5, 1),
  (3, a6, 1), (8, a6, 1), (13, a6, 1),(18, a6, 1),
  (4, a7, 1), (9, a7, 1), (14, a7, 1),(19, a7, 1),
  (5, a8, 1), (10, a8, 1), (15, a8, 1),(20, a8, 1)
  
  COV=Eta1*Eta2*Eta3*Eta4*Beta2*Beta3*Beta4*Beta5') 
fit_mirt <- mirt::mirt(dt_mirt, model=model_mirt_m3, itemtype="graded",method="MHRM" )  
#coef(fit_mirt, simplify=T)



###write
mfit[[1]] <- fitMeasures(fit_num)["rmsea.robust"]
mfit[[2]] <- fitMeasures(fit_ord)["rmsea"]
mfit[[3]] <- mirt::extract.mirt(fit_mirt, "RMSEA")

mfit_time[[1]] <- fit_num@timing["total"][[1]]
mfit_time[[2]] <- fit_ord@timing["total"][[1]]
mfit_time[[3]] <- fit_mirt@time[["TOTAL:"]]

models[[1]] <- fit_num
models[[2]] <- fit_ord
models[[3]] <- fit_mirt
names(mfit) <- names(mfit_time) <- names(models) <- c("MLR","WLS","MML")

mfit_m3 = mfit
mfit_time_m3 = mfit_time
models_m3 = models


### All models list
mirtmodels = list(1,model_mirt_m2,model_mirt_m3)
mfit = list(mfit_m1,mfit_m2,mfit_m3)
mfit_time = list(mfit_time_m1,mfit_time_m2,mfit_time_m3)
models = list(models_m1,models_m2,models_m3)
names(mfit) <- names(mfit_time) <- names(models) <- names(mirtmodels) <- c("unidim","multidim","PIEG")


################################################################################
#################################  Flucs  ######################################
################################################################################

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
for(i in 1:3){
  
  #lavaan num
  start_time <- Sys.time()
  fluc_cat <- strucchange::gefp(models[[i]][[1]],fit=NULL, order.by = dt[,"geslacht"]) #categorical
  pval_cat <- strucchange::sctest(fluc_cat, functional = strucchange::catL2BB(fluc_cat))$p.value #LMuo
  time_cat <- difftime(Sys.time(), start_time, units="secs")
  
  start_time <- Sys.time()
  fluc_met <- strucchange::gefp(models[[i]][[1]],fit=NULL, order.by = dt[,"leeftijd"]) #metric
  pval_met <- strucchange::sctest(fluc_met, functional = strucchange::maxBB)$p.value #DM
  time_met <-  difftime(Sys.time(), start_time, units="secs") 
  
  start_time <- Sys.time()
  fluc_ord <- strucchange::gefp(models[[i]][[1]],fit=NULL, order.by = dt[,"sted"]) #ordered
  pval_ord <- strucchange::sctest(fluc_ord, functional = strucchange::ordwmax(fluc_ord))$p.value #WDMo
  time_ord <-  difftime(Sys.time(), start_time, units="secs")
  
  
  times_mlr = mean(as.numeric(time_cat),as.numeric(time_met),as.numeric(time_ord))
  flucs_mlr = list(fluc_cat,fluc_met,fluc_ord)
  pvals_mlr = list(pval_cat,pval_met,pval_ord)
  names(flucs_mlr) = names(pvals_mlr) = c("categorical","metric","ordinal")
  
  
  
  #lavaan ord
  start_time <- Sys.time()
  fluc_cat <- strucchange::gefp(models[[i]][[2]],fit=NULL, order.by = dt[,"geslacht"]) #categorical
  pval_cat <- strucchange::sctest(fluc_cat, functional = strucchange::catL2BB(fluc_cat))$p.value #LMuo
  time_cat <- difftime(Sys.time(), start_time, units="secs")
  
  start_time <- Sys.time()
  fluc_met <- strucchange::gefp(models[[i]][[2]],fit=NULL, order.by = dt[,"leeftijd"]) #metric
  pval_met <- strucchange::sctest(fluc_met, functional = strucchange::maxBB)$p.value #DM
  time_met <-  difftime(Sys.time(), start_time, units="secs") 
  
  start_time <- Sys.time()
  fluc_ord <- strucchange::gefp(models[[i]][[2]],fit=NULL, order.by = dt[,"sted"]) #ordered
  pval_ord <- strucchange::sctest(fluc_ord, functional = strucchange::ordwmax(fluc_ord))$p.value #WDMo
  time_ord <-  difftime(Sys.time(), start_time, units="secs")
  
  
  times_wls = mean(as.numeric(time_cat),as.numeric(time_met),as.numeric(time_ord))
  flucs_wls = list(fluc_cat,fluc_met,fluc_ord)
  pvals_wls = list(pval_cat,pval_met,pval_ord)
  names(flucs_wls) = names(pvals_wls) = c("categorical","metric","ordinal")
  
  #mirt
  model_mirt = mirtmodels[[i]]
  start_time <- Sys.time()
  fluc_cat <- strucchange::gefp(models[[i]][[3]],fit=NULL, scores=estfun.MHRM, order.by = dt[,"geslacht"]) #categorical
  pval_cat <- strucchange::sctest(fluc_cat, functional = strucchange::catL2BB(fluc_cat))$p.value #LMuo
  time_cat <- difftime(Sys.time(), start_time, units="secs")
  
  start_time <- Sys.time()
  fluc_met <- strucchange::gefp(models[[i]][[3]],fit=NULL, scores=estfun.MHRM, order.by = dt[,"leeftijd"]) #metric
  pval_met <- strucchange::sctest(fluc_met, functional = strucchange::maxBB)$p.value #DM
  time_met <-  difftime(Sys.time(), start_time, units="secs") 
  
  start_time <- Sys.time()
  fluc_ord <- strucchange::gefp(models[[i]][[3]],fit=NULL, scores=estfun.MHRM, order.by = dt[,"sted"]) #ordered
  pval_ord <- strucchange::sctest(fluc_ord, functional = strucchange::ordwmax(fluc_ord))$p.value #WDMo
  time_ord <-  difftime(Sys.time(), start_time, units="secs")
  
  
  times_mml = mean(as.numeric(time_cat),as.numeric(time_met),as.numeric(time_ord))
  flucs_mml = list(fluc_cat,fluc_met,fluc_ord)
  pvals_mml = list(pval_cat,pval_met,pval_ord)
  names(flucs_mml) = names(pvals_mml) = c("categorical","metric","ordinal")
  
  #write
  flucs_time[[i]] = list(times_mlr,times_wls,times_mml)
  flucs[[i]] = list(flucs_mlr,flucs_wls,flucs_mml)
  pvals[[i]] = list(pvals_mlr,pvals_wls,pvals_mml)
  names(flucs_time[[i]])  = names(flucs[[i]])= names(pvals[[i]]) = c("MLR","WLS","MML")
  
}
names(flucs_time) = names(flucs) = names(pvals) = names(models)
results = list(models,mfit,mfit_time,flucs_time,flucs,pvals)
save(dt,results, file="/dss/dsshome1/0B/ra35tik2/paper3/res/240605_results_liss.RData")























################################################################################
################################  Sandbox  #####################################
################################################################################

################################################################################
#############################  PIEG model SE  ##################################
################################################################################
model_multi <- "

etaSE1 =~ 1*cp08a070 + 1*cp08a071 + 1*cp08a072 + 1*cp08a073 + 1*cp08a074 + 1*cp08a075 + 1*cp08a076 + 1*cp08a077 + 1*cp08a078 + 1*cp08a079 #+ 1*cp08a080 + 1*cp08a081 + 1*cp08a082  
etaSE2 =~ 1*cp09b070 + 1*cp09b071 + 1*cp09b072 + 1*cp09b073 + 1*cp09b074 + 1*cp09b075 + 1*cp09b076 + 1*cp09b077 + 1*cp09b078 + 1*cp09b079 #+ 1*cp09b080 + 1*cp09b081 + 1*cp09b082
etaSE3 =~ 1*cp11d070 + 1*cp11d071 + 1*cp11d072 + 1*cp11d073 + 1*cp11d074 + 1*cp11d075 + 1*cp11d076 + 1*cp11d077 + 1*cp11d078 + 1*cp11d079 #+ 1*cp11d080 + 1*cp11d081 + 1*cp11d082
etaSE4 =~ 1*cp13f070 + 1*cp13f071 + 1*cp13f072 + 1*cp13f073 + 1*cp13f074 + 1*cp13f075 + 1*cp13f076 + 1*cp13f077 + 1*cp13f078 + 1*cp13f079 #+ 1*cp13f080 + 1*cp13f081 + 1*cp13f082


deltaSE071=~1*cp08a071+1*cp09b071+1*cp11d071+1*cp13f071
deltaSE072=~1*cp08a072+1*cp09b072+1*cp11d072+1*cp13f072
deltaSE073=~1*cp08a073+1*cp09b073+1*cp11d073+1*cp13f073
deltaSE074=~1*cp08a074+1*cp09b074+1*cp11d074+1*cp13f074
deltaSE075=~1*cp08a075+1*cp09b075+1*cp11d075+1*cp13f075
deltaSE076=~1*cp08a076+1*cp09b076+1*cp11d076+1*cp13f076
deltaSE077=~1*cp08a077+1*cp09b077+1*cp11d077+1*cp13f077
deltaSE078=~1*cp08a078+1*cp09b078+1*cp11d078+1*cp13f078
deltaSE079=~1*cp08a079+1*cp09b079+1*cp11d079+1*cp13f079

#deltaSE080=~1*cp08a080+1*cp09b080+1*cp11d080+1*cp13f080
#deltaSE081=~1*cp08a081+1*cp09b081+1*cp11d081+1*cp13f081
#deltaSE082=~1*cp08a082+1*cp09b082+1*cp11d082+1*cp13f082
"
fit_ord <- lavaan::cfa(model = model_multi, data=dt, ordered = TRUE, estimator="WLS", std.lv=F)
summary(fit_ord, fit.measures=T) #331 params #RMSEA=0.052
fit_ord@timing["total"]


dt_mirt = dt[,unlist(vars)]
model_mirt <- mirt::mirt.model('
  Eta1 = 1-10
  Eta2 = 11-20
  Eta3 = 21-30
  Eta4 = 31-40
  Beta2 = 2,12,22,32
  Beta3 = 3,13,23,33
  Beta4 = 4,14,24,34
  Beta5 = 5,15,25,35
  Beta6 = 6,16,26,36
  Beta7 = 7,17,27,37
  Beta8 = 8,18,28,38
  Beta9 = 9,19,29,39
  Beta10 = 10,20,30,40
  
  CONSTRAIN = (1, a1, 1), (2, a1, 1), (3, a1, 1), (4, a1, 1),(5, a1, 1),(6, a1, 1),(7, a1, 1),(8, a1, 1),(9, a1, 1),(10, a1, 1),
  (11, a2, 1), (12, a2, 1), (13, a2, 1), (14, a2, 1),(15, a2, 1),(16, a2, 1),(17, a2, 1),(18, a2, 1),(19, a2, 1),(20, a2, 1),
  (21, a3, 1), (22, a3, 1), (23, a3, 1), (24, a3, 1),(25, a3, 1),(26, a3, 1),(27, a3, 1),(28, a3, 1),(29, a3, 1),(30, a3, 1),
  (31, a4, 1), (32, a4, 1), (33, a4, 1), (34, a4, 1),(35, a4, 1),(36, a4, 1),(37, a4, 1),(38, a4, 1),(39, a4, 1),(40, a4, 1),
  (2, a5, 1), (12, a5, 1), (22, a5, 1), (32, a5, 1), 
  (3, a6, 1), (13, a6, 1), (23, a6, 1), (33, a6, 1), 
  (4, a7, 1), (14, a7, 1), (24, a7, 1), (34, a7, 1),
  (5, a8, 1), (15, a8, 1), (25, a8, 1), (35, a8, 1),
  (6, a9, 1), (16, a9, 1), (26, a9, 1), (36, a9, 1),
  (7, a10, 1), (17, a10, 1), (27, a10, 1), (37, a10, 1),
  (8, a11, 1), (18, a11, 1), (28, a11, 1), (38, a11, 1),
  (9, a12, 1), (19, a12, 1), (29, a12, 1), (39, a12, 1),
  (10, a13, 1), (20, a13, 1), (30, a13, 1), (40, a13, 1)

  COV=Eta1*Eta2*Eta3*Eta4*Beta2*Beta3*Beta4*Beta5*Beta6*Beta7*Beta8*Beta9*Beta10') 


fit_mirt <- mirt::mirt(dt_mirt, model=model_mirt, itemtype="graded",method="MHRM" )  
coef(fit_mirt, simplify=T)
#params_mirt = mirt::mirt(dt_mirt, model=model_mirt, itemtype="graded",pars="values" )  
#nrow(params_mirt[params_mirt$est==T,])
















################################################################################
##########################  Large Lavaan model  ################################
################################################################################

model_multi <- "
etaLS1=~ cp08a014 + cp08a015 + cp08a016 + cp08a017+ cp08a018
etaLS2=~ cp09b014 + cp09b015 + cp09b016 + cp09b017+ cp09b018
etaLS3=~ cp11d014 + cp11d015 + cp11d016 + cp11d017+ cp11d018

etaSE1 =~ cp08a070 + cp08a071 + cp08a072 + cp08a073 + cp08a074 + cp08a075 + cp08a076 + cp08a077 + cp08a078 + cp08a079 #+ cp08a080 + cp08a081 + cp08a082  
etaSE2 =~ cp09b070 + cp09b071 + cp09b072 + cp09b073 + cp09b074 + cp09b075 + cp09b076 + cp09b077 + cp09b078 + cp09b079 #+ cp09b080 + cp09b081 + cp09b082
etaSE3 =~ cp11d070 + cp11d071 + cp11d072 + cp11d073 + cp11d074 + cp11d075 + cp11d076 + cp11d077 + cp11d078 + cp11d079 #+ cp11d080 + cp11d081 + cp11d082


thetaLS =~ etaLS1 + etaLS2 + etaLS3
thetaSE =~ etaSE1 + etaSE2 + etaSE3

etaLS1 ~~ 0*etaSE1
etaLS1 ~~ 0*etaSE2
etaLS1 ~~ 0*etaSE3
etaLS2 ~~ 0*etaSE1
etaLS2 ~~ 0*etaSE2
etaLS2 ~~ 0*etaSE3
etaLS3 ~~ 0*etaSE1
etaLS3 ~~ 0*etaSE2
etaLS3 ~~ 0*etaSE3


#Regression
thetaLS  ~ thetaSE
"

fit_ord <- lavaan::cfa(model = model_multi, data=dt, ordered = TRUE, estimator="WLS", std.lv=T) #control=list(iter.max=100000)
summary(fit_ord,fit.measures=T) #RMSEA= 0.056!
fit_ord@timing[["total"]]

#fit_mml <- lavaan::cfa(model = model_multi, data=dt, ordered = TRUE, estimator="MML", std.lv=T) #takes forever...





####
specific <- c(rep(1,5),rep(2,5), rep(3,5))
model <- "G = 1-15
          CONSTRAIN = (1, a1, a2), (2, a1, a2), (3, a1, a2), (4, a1, a2),(5, a1, a2),
            (6, a1, a3), (7, a1, a3), (8, a1, a3), (9, a1, a3), (10, a1, a3),
            (11, a1, a4), (12, a1, a4), (13, a1, a4), (14, a1, a4),(15, a1, a4)
          COV = S1*S1, S2*S2, S3*S3"
fit_bfactor <- mirt::bfactor(dt_mirt, specific, model)
params_bfactor = mirt::bfactor(dt_mirt, specific, model, pars="values")
coef(fit_bfactor, simplify=T)



model_multi <- "
etaLS1=~ 1*cp08a014 + 1*cp08a015 + 1*cp08a016 + 1*cp08a017+ 1*cp08a018
etaLS2=~ 1*cp09b014 + 1*cp09b015 + 1*cp09b016 + 1*cp09b017+ 1*cp09b018
etaLS3=~ 1*cp11d014 + 1*cp11d015 + 1*cp11d016 + 1*cp11d017+ 1*cp11d018

etaSE1 =~ cp08a070 + cp08a071 + cp08a072 + cp08a073 + cp08a074 + cp08a075 + cp08a076 + cp08a077 + cp08a078 + cp08a079 #+ cp08a080 + cp08a081 + cp08a082  
etaSE2 =~ cp09b070 + cp09b071 + cp09b072 + cp09b073 + cp09b074 + cp09b075 + cp09b076 + cp09b077 + cp09b078 + cp09b079 #+ cp09b080 + cp09b081 + cp09b082
etaSE3 =~ cp11d070 + cp11d071 + cp11d072 + cp11d073 + cp11d074 + cp11d075 + cp11d076 + cp11d077 + cp11d078 + cp11d079 #+ cp11d080 + cp11d081 + cp11d082

##reference item 16
deltaLS014=~1*cp08a014+1*cp09b014+1*cp11d014
deltaLS015=~1*cp08a015+1*cp09b015+1*cp11d015
deltaLS017=~1*cp08a017+1*cp09b017+1*cp11d017
deltaLS018=~1*cp08a018+1*cp09b018+1*cp11d018

thetaLS =~ etaLS1 + etaLS2 + etaLS3
thetaSE =~ etaSE1 + etaSE2 + etaSE3

etaLS1 ~~ 0*etaSE1
etaLS1 ~~ 0*etaSE2
etaLS1 ~~ 0*etaSE3
etaLS2 ~~ 0*etaSE1
etaLS2 ~~ 0*etaSE2
etaLS2 ~~ 0*etaSE3
etaLS3 ~~ 0*etaSE1
etaLS3 ~~ 0*etaSE2
etaLS3 ~~ 0*etaSE3
deltaLS014 ~~ 0*etaSE1
deltaLS015 ~~ 0*etaSE1
deltaLS017 ~~ 0*etaSE1
deltaLS018 ~~ 0*etaSE1
deltaLS014 ~~ 0*etaSE2
deltaLS015 ~~ 0*etaSE2
deltaLS017 ~~ 0*etaSE2
deltaLS018 ~~ 0*etaSE2
deltaLS014 ~~ 0*etaSE3
deltaLS015 ~~ 0*etaSE3
deltaLS017 ~~ 0*etaSE3
deltaLS018 ~~ 0*etaSE3

#Regression
thetaLS  ~ thetaSE
"
fit_ord <- lavaan::cfa(model = model_multi, data=dt, ordered = TRUE, estimator="WLS", std.lv=T)

