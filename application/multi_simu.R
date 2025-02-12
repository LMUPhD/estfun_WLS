

#PIEG model without mean structure and with standard LV-variance
datagen <- function(model,mode="all",dif_items=9,times=1,ID=250,korr="random",schwellen=4,items=3,latvar=3,rmsea_cutoff=.05){ 
  l <- schwellen #Anzahl an Schwellen
  m <- items #Anzahl an Items pro latente Variable
  n <- latvar #Anzahl an latente Variablen
  nlatvar = n 
  
  randparams <- function(N,num.min,rnd=2,num.max=0.99) {
    random.numbers <- round(runif(min=num.min,max=num.max, N),rnd)
    random.signs   <- sample(c(1,-1), N, replace=TRUE)      
    return(random.numbers * random.signs)
  }
  
  if(mode=="none"){dif_items = 0}
  
  schwellen_probs <- function(schwellen){
    if(schwellen == 6) return(c(0.10,0.25,0.40,0.60,0.75,0.90)  )
    if(schwellen == 5) return(c(0.10,0.30,0.50,0.70,0.90) )
    if(schwellen == 4) return(c(0.15,0.40,0.60,0.85)  )
    if(schwellen == 3) return(c(0.25,0.5,0.75)  )
    if(schwellen == 2) return(c(0.30,0.70)  )
    if(schwellen == 1) return(c(0.50))
    else stop("Error: Function can only produce data with 2 to 7 categories.")
  }
  
  library(plyr)
  library(psych)
  library(lavaan)
  library(faux)
  library(resample)
  
  for(i in c("saved_pvalues","saved_rmseas","saved_kappas","saved_data","saved_covs","saved_betas","saved_vars","saved_latvars")){assign(i,list())}
  
  #unterschiedliche modes:
  if(mode == "thresholds" | mode == "betas" | mode == "none"){
    if(korr=="high"){mulmin=0.8;mulmax=0.9;errvar=.1} else if(korr=="low"){mulmin=0.1;mulmax=0.3;errvar=1} else if(korr=="random"){mulmin=0.1;mulmax=0.9;errvar=runif(1,min = .3, max = .9)}
    mulpis <- randparams(num.min=mulmin,num.max=mulmax, nlatvar) 
    
    if(mode == "thresholds" | mode == "none"){
      beta <- sapply(1:n, function(x)  c(1,round(runif(min=0.3,max=1.6, m-1),2))) #Discrimination parameter
    }
    
    if(mode == "betas" | mode == "none"){
      kappa_shift <- round(rnorm(n = n*m, mean=0, sd = 1),2) #"Erwartungswert"-Shift pro Item
      perz_kappa <- matrix(NA,ncol = l, nrow = n*m) #Perzentile, um die kappa-Parameter zu bestimmen
      for(i in 1:(n*m)){perz_kappa[i,] <- round(schwellen_probs(l) + randparams(num.min=0.01,num.max=0.05, l) ,2)}
    }
  }

  
  #####################################################################################
  ################# Algorithmus
  #####################################################################################
  h=1;g=1
  while(h<(times+1) & g<11){ #max. 10 Widerholungen...

    
    #Alle Parameter unterschiedlich
    if(mode == "all"){
      if(korr=="high"){mulmin=0.8;mulmax=0.9;errvar=.1} else if(korr=="low"){mulmin=0.1;mulmax=0.3;errvar=1} else if(korr=="random"){mulmin=0.1;mulmax=0.9;errvar=runif(1,min = .3, max = .9)}
      mulpis <- randparams(num.min=mulmin,num.max=mulmax, nlatvar) 
    }
    if(mode=="betas" | mode == "all"){
      beta <- sapply(1:n, function(x)  c(1,round(runif(min=0.3,max=1.6, m-1),2))) #Discrimination parameter
    }
    if(mode == "thresholds" | mode == "all"){
      kappa_shift <- round(rnorm(n = n*m, mean=0, sd = 1),2) #"Erwartungswert"-Shift pro Item
      perz_kappa <- matrix(NA,ncol = l, nrow = n*m) #Perzentile, um die kappa-Parameter zu bestimmen
      for(i in 1:(n*m)){perz_kappa[i,] <- round(schwellen_probs(l) + randparams(num.min=0.01,num.max=0.05, l) ,2)}
    }
    
    #Variablen erstellen
    

    Psi <- rnorm(ID,0,1)
    
    data_mtx <- mulpis %*% t(Psi)  
    data_mtx <- t(data_mtx) + matrix( rnorm(ID*nlatvar,mean=0,sd=sqrt(errvar)), ID, nlatvar)  #kleine Fehlervarianz --> hohe Korrelation 
    
    
    ## Variablennamen 

    Eta <- data_mtx[,1:n] 

    
    
    ## Kappa Matrix erstellen
    kappa <- matrix(ncol = n*m,nrow = l);v=1 #ncol = L-1 Variablen * Items ; nrow = schwellen
    for(t in 1:n){ #time points
      for(i in 1:m){ # items
        quant <- as.vector(quantile(beta[i,t]*Eta[,t]  , probs= perz_kappa[v,])) - kappa_shift[v] 
        kappa[1:l,v] <- quant 
        v <- v+1
      }}
    kappa <- round(kappa,2)
    
    
    
    
    ## Number of dif items
    if(h>1 & dif_items > 0 & dif_items < 5){
      if(mode == "thresholds" | mode == "all"){
        kappa_new = kappa
        kappa = saved_kappas[[h-1]]
        kappa[,1:dif_items] = kappa_new[,1:dif_items]
      }
      if(mode=="betas" | mode == "all"){ #at least dif_items=2 needed 
        beta_new = beta
        beta = saved_betas[[h-1]]
        beta[1:dif_items] = beta[1:dif_items]
      }
    }
    
    
    
    
    ## Daten generieren
    probitY <- array(NA,c(ID,l,m,n)) #(ID, Schwellen, Items, Lat.Var.)
    PY <-  array(NA,c(ID,l,m,n)) 
    Pis <-  array(NA,c(ID,l+1,m,n))
    Y <- array(NA,c(ID,m,n));v=1 
    
    for (t in 1:n) { #Latvars
      for (i in 1:m) { #Items 
        for (k in 1:l) { #Schwellen
          probitY[,k,i,t] <-   beta[i,t]*Eta[,t] - kappa[k,v] 
          PY[,k,i,t] <- VGAM::probitlink(probitY[,k,i,t],inverse=T) 
          if (k == 1) Pis[,k,i,t] <- 1 - PY[,k,i,t] #P(Y=0)
          if (k > 1) Pis[,k,i,t] <- PY[,k-1,i,t] - PY[,k,i,t] 
          if (k == l) Pis[,k+1,i,t] <- PY[,k,i,t]
        }
        for(j in 1:ID){Y[j,i,t] <- sample(1:(l+1), size=1, replace=TRUE, prob=sapply(1:(l+1), function(x) Pis[j,x,i,t])  )} 
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
    
    
    ## Model schätzen
    fit_ord <- tryCatch({lavaan::cfa(model = model, data=table, ordered = output, estimator="WLS", do.fit=T, std.lv=T, control=list(iter.max=500))},warning = function(w){return(NA)},error = function(e){return(NA)})
    
    ##Ergebnisse Zusammenfügen!
    if( suppressWarnings(!is.na(fit_ord))  ){
      
      
      if(fitMeasures(fit_ord)["rmsea"] < rmsea_cutoff  ){ 
        
        
        
        saved_pvalues[[h]] <- fitMeasures(fit_ord)["pvalue"]; names(saved_pvalues)[[h]] <- paste0("pvalue",h)
        saved_rmseas[[h]] <- fitMeasures(fit_ord)["rmsea"]; names(saved_rmseas)[[h]] <- paste0("rmsea",h)
        saved_vars[[h]] <- colVars(data_mtx); names(saved_vars)[[h]] <- paste0("var",h)
        saved_kappas[[h]] <- kappa; names(saved_kappas)[[h]] <- paste0("kappa",h)
        saved_covs[[h]] <- cov(data_mtx); names(saved_covs)[[h]] <- paste0("cov",h)
        saved_betas[[h]] <- beta; names(saved_betas)[[h]] <- paste0("beta",h)
        saved_latvars[[h]] <- data_mtx
        saved_data[[h]] <- table; names(saved_data)[[h]] <- paste0("data",h)
        
        h=h+1 
        
      } else {next}
    } else {next}
  } 
  
  fits <- list(saved_pvalues,saved_rmseas,saved_vars,saved_kappas,saved_covs,saved_betas,saved_latvars,saved_data)
  names(fits) <- c("pvalues","rmseas","vars","kappas","covs","betas","latvars","data")
  return(fits)
}


