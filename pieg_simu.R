setwd("/dss/dsshome1/lxc08/di67gug/paper1/simulation/additive_model/sctree")

#PIEG model without mean structure and with standard LV-variance
#Simu 1 all parameters different

datagen_pieg <- function(model,times=4,ID=250,schwellen=4,items=3,timepoints=3,rmsea_cutoff=.04){ 
  
  
  randparams <- function(N,num.min,rnd=2,num.max=0.99) {
    random.numbers <- round(runif(min=num.min,max=num.max, N),rnd)
    random.signs   <- sample(c(1,-1), N, replace=TRUE)      
    return(random.numbers * random.signs)
  }
  
  
  
  library(plyr)
  library(psych)
  library(lavaan)
  library(faux)
  library(metaSEM)
  library(resample)
  
  for(i in c("saved_pvalues","saved_rmseas","saved_kappas","saved_data","saved_covs","saved_vars")){assign(i,list())}
  
  #####################################################################################
  ################# Algorithmus
  #####################################################################################
  h=1;g=1
  while(h<(times+1)){ 
    l <- schwellen #Anzahl an Schwellen
    m <- items #Anzahl an Items pro Level-1 latente Variable
    n <- timepoints #Anzahl an Level-1 latente Variablen
    nlatvar = n+(m-1)
    
    ## Simu 1: Alle Parameter unterschiedlich, keine random data (obviously: keine konfundierung)
    kappa_shift <- round(rnorm(n = n*m, mean=0, sd = 0.6),2) #"Erwartungswert"-Shift pro Item
    perz_kappa <- c(0.10,0.35,0.65,0.90)#Perzentile, um die kappa-Parameter zu bestimmen
    
    #Variablen erstellen

    
    cm <- matrix(randparams(nlatvar^2,num.min=0.1,num.max=0.4), ncol=nlatvar)
    cmat <- t(cm) %*% cm

    data_mtx <- rnorm_multi(n = ID, 
                       mu = 0, #no meanstructure
                       r = cmat, 
                       varnames = sapply(1:nlatvar,function(y){paste0("latvar",y)}),
                       empirical = F,
                       as.matrix = T)

    ## Variablennamen 
    Eta <- data_mtx[,1:n] 
    Beta <- cbind(rep(0,ID),  data_mtx[,(n+1):nlatvar] )
    

    
    ## Kappa Matrix erstellen
    kappa <- matrix(ncol = n*m,nrow = l);v=1 #ncol = L-1 Variablen * Items ; nrow = schwellen
    for(t in 1:n){ #time points
      for(i in 1:m){ # items
        quant <- as.vector(quantile(Eta[,t]  + Beta[,i], probs= perz_kappa)) - kappa_shift[v] 
        kappa[1:l,v] <- quant 
        v <- v+1
      }}
    kappa <- round(kappa,2)
  
    
    ## Daten generieren
    probitY <- array(NA,c(ID,l,m,n)) #(ID, Schwellen, Items pro L-1 Var., L-1 Var.)
    PY <-  array(NA,c(ID,l,m,n)) 
    Pis <-  array(NA,c(ID,l+1,m,n))
    Y <- array(NA,c(ID,m,n));v=1;j=1 #(ID, Schwellen, Items pro L-1 Var., L-1 Var.)
    
    for (t in 1:n) { #Time points
      for (i in 1:m) { #Items 
        for (k in 1:l) { #Schwellen
          probitY[,k,i,t] <-  Eta[,t] + Beta[,i] - kappa[k,v] 
          PY[,k,i,t] <- probitlink(probitY[,k,i,t],inverse=T) 
          if (k == 1) Pis[,k,i,t] <- 1 - PY[,k,i,t] #P(Y=0)
          if (k > 1) Pis[,k,i,t] <- PY[,k-1,i,t] - PY[,k,i,t] 
          if (k == l) Pis[,k+1,i,t] <- PY[,k,i,t]
        }
        for(j in 1:ID){Y[j,i,t] <- sample(1:(l+1), size=1, replace=TRUE, prob=c(Pis[j,1,i,t],Pis[j,2,i,t],Pis[j,3,i,t],Pis[j,4,i,t],Pis[j,5,i,t]))} 
        v <- v+1
      }
    }
    

    #### Tabelle erzeugen & Deskriptive Statistik
    table <- matrix(ncol = n*m,nrow = ID)
    v=1
    for (t in 1:n){
      for (i in 1:m){
        table[,v] <- Y[,i,t]
        v <- v+1
      }
    }
    
    
    ############################           Ueberpruefen            ##########################################
    
    output <- sapply(1:(n*m),function(y){paste0("simuvar",y)})
    table <- as.data.frame(table)
    colnames(table) <- c(output)
    
   
    
    print(paste("Iteration ",g)) #wieviele Iterationen braucht es?
    g=g+1
    ##Ergebnisse Zusammenfügen!
    #fit_ord <- lavaan::cfa(model = model, data=table, ordered = output, estimator="WLS", do.fit=T, std.lv=F, control=list(iter.max=500))
    #summary(fit_ord,fit.measures=T)
    #View(lavInspect(fit_ord, "cov.lv"))
    
    fit_ord <- tryCatch({lavaan::cfa(model = model, data=table, ordered = output, estimator="WLS", do.fit=T, std.lv=F, control=list(iter.max=500))},warning = function(w){return(NA)},error = function(e){return(NA)})
    
    if( suppressWarnings(!is.na(fit_ord))  ){
      if(fitMeasures(fit_ord)["rmsea"] < rmsea_cutoff){
        
        saved_pvalues[[h]] <- fitMeasures(fit_ord)["pvalue"]; names(saved_pvalues)[[h]] <- paste0("pvalue",h)
        saved_rmseas[[h]] <- fitMeasures(fit_ord)["rmsea"]; names(saved_rmseas)[[h]] <- paste0("rmsea",h)
        saved_vars[[h]] <- colVars(data_mtx); names(saved_vars)[[h]] <- paste0("var",h)
        saved_kappas[[h]] <- kappa; names(saved_kappas)[[h]] <- paste0("kappa",h)
        saved_covs[[h]] <- cov(data_mtx); names(saved_covs)[[h]] <- paste0("cov",h)
        saved_data[[h]] <- table; names(saved_data)[[h]] <- paste0("data",h)
        
        h=h+1 
        
      } else {next}
    } else {next}
  } 
  
  fits <- list(saved_pvalues,saved_rmseas,saved_vars,saved_kappas,saved_covs,saved_data)
  names(fits) <- c("pvalues","rmseas","vars","kappas","covs","data")
  return(fits)
}


fits <- datagen1(
  model = '
  Eta1 =~ 1*simuvar1 + 1*simuvar2 + 1*simuvar3 
  Eta2 =~ 1*simuvar4 + 1*simuvar5 + 1*simuvar6 
  Eta3 =~ 1*simuvar7 + 1*simuvar8 + 1*simuvar9 
  Beta2 =~ 1*simuvar2 + 1*simuvar5 + 1*simuvar8
  Beta3 =~ 1*simuvar3 + 1*simuvar6 + 1*simuvar9',
  rmsea_cutoff = .03,
  ID=500,
  times=4
) 

goodfits1 <- fits[["data"]][[1]] 
goodfits2 <- fits[["data"]][[2]]
goodfits3 <- fits[["data"]][[3]] 
goodfits4 <- fits[["data"]][[4]]

#############################################################################
#######################   Partitioning Variables   ##########################
#############################################################################
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

goodfits1$cat2 <- replicate(ID,sample(c(1,2,3,4),1))
goodfits2$cat2 <- replicate(ID,sample(c(1,2,3,4),1))
goodfits3$cat2 <- replicate(ID,sample(c(1,2),1))
goodfits4$cat2 <- replicate(ID,sample(c(3,4),1))





###### Whole Dataset
simu1 <- rbind(goodfits1,goodfits2,goodfits3,goodfits4)
simu1 <- simu1[sample(nrow(simu1)),] ## Shuffle rows

##### Add random partitioning variables
simu1$rand1 <-  replicate(nrow(simu1),sample(c(1,2,3,4,5),1,prob=c(0.1,0.2,0.4,0.2,0.1))) #ordinal
simu1$rand2 <-  round(runif(nrow(simu1),min=1,max=5),0)                                   #ordinal
simu1$rand3 <-  replicate(nrow(simu1),sample(c(0,1),1,prob=c(0.7,0.3)))                   #dichotomous
simu1$rand4 <-  round(runif(nrow(simu1),min=100,max=400),2)                               #numeric (univ)
simu1$rand5 <-  round(rnorm(nrow(simu1),mean=30,sd=5),4)                                  #numeric (norm)

  
  
##### Prepocessing 
simu1$cat1 <- as.factor(simu1$cat1) 
simu1$cat2 <- as.factor(simu1$cat2) 
simu1$rand1 <- as.factor(simu1$rand1) 
simu1$rand2 <- as.factor(simu1$rand2) 
simu1$rand3 <- as.factor(simu1$rand3) 

simu1$num1 <- as.numeric(simu1$num1) 
simu1$rand4 <- as.numeric(simu1$rand5) 
  
##### Es hat geklappt!!!
save(simu1,fits,file="simu1_data_pieg.RData")

















