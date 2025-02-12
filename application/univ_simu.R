
#no clear structure --> forest needeed
datagen <- function(times=1,mode='all',dif_items=5,ID=250,schwellen=6,items=5,rmsea_cutoff=.05){ 
  l <- schwellen #Anzahl an Schwellen
  m <- items #Anzahl an Items pro latente Variable
  nlatvar = 1
  
  model = paste(sapply(1:items, function(x) paste0("simuvar",x," +")), collapse = "")
  model = paste0("Eta1 =~ ",substr(model,1,nchar(model)-2))
  
  
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
  library(resample)
  
  for(i in c("saved_pvalues","saved_rmseas","saved_kappas","saved_betas","saved_data","saved_vars","saved_latvars","saved_means","saved_model")){assign(i,list())}
  
  #unterschiedliche modes:
  if(mode == "thresholds" | mode == "betas" | mode == "none"){
    
    var <- round(runif(min=0.3,max=1.2, 1),2)
    
    if(mode == "thresholds" | mode == "none"){
      beta <- c(1,round(runif(min=0.3,max=1.6, m-1),2)) #Discrimination parameter
    }
    
    if(mode == "betas" | mode == "none"){
      kappa_shift <- round(rnorm(n = m, mean=0, sd = 1),2) #"Erwartungswert"-Shift pro Item
      perz_kappa <- matrix(NA,ncol = l, nrow = m) #Perzentile, um die kappa-Parameter zu bestimmen
      for(i in 1:(m)){perz_kappa[i,] <- round(schwellen_probs(l) + randparams(num.min=0.01,num.max=0.05, l) ,2)}
    }
  }
  
  
  #####################################################################################
  ################# Algorithmus
  #####################################################################################
  h=1;g=1
  while(h<(times+1) & g<11){ 
    
    #Alle Parameter unterschiedlich
    if(mode == "all"){
      var <- round(runif(min=0.3,max=1.2, 1),2)
    }
    if(mode=="betas" | mode == "all"){
      beta <- c(1,round(runif(min=0.3,max=1.6, m-1),2)) #Discrimination parameter
    }
    if(mode == "thresholds" | mode == "all"){
      kappa_shift <- round(rnorm(n = m, mean=0, sd = 0.5),2) #"Erwartungswert"-Shift pro Item
      perz_kappa <- matrix(NA,ncol = l, nrow = m) #Perzentile, um die kappa-Parameter zu bestimmen
      for(i in 1:(m)){perz_kappa[i,] <- round(schwellen_probs(l) + randparams(num.min=0.01,num.max=0.05, l) ,2)}
    }
    

    

    #Variable erstellen
    Psi <- rnorm(ID,0,sqrt(var)) + rnorm(ID*nlatvar,mean=0,sd=sqrt(runif(1,min = .3, max = .9))) #psi erstmal mittelwert 0 + Fehlervarianz
    

    ###### Categorical Variables
    
    ## Kappa Matrix erstellen
    kappa <- matrix(ncol = m,nrow = l)  #ncol = L-1 Variablen * Items ; nrow = schwellen
    for(i in 1:m){ # items
      quant <- as.vector(quantile(beta[i]*Psi , probs= perz_kappa[i,])) - kappa_shift[i] 
      kappa[,i] <- quant
    }
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
    probitY <- array(NA,c(ID,l,m)) #(ID, Schwellen, Items pro L-1 Var., L-1 Var.)
    PY <-  array(NA,c(ID,l,m)) 
    Pis <-  array(NA,c(ID,l+1,m))
    Y_ord <- array(NA,c(ID,m));v=1 #(ID, Schwellen, Items pro L-1 Var., L-1 Var.)
    
    
    for (i in 1:m) { #Items 
      for (k in 1:l) { #Schwellen
        probitY[,k,i] <-  beta[v]*Psi  - kappa[k,v] 
        PY[,k,i] <- VGAM::probitlink(probitY[,k,i],inverse=T) 
        if (k == 1) Pis[,k,i] <- 1 - PY[,k,i] #P(Y=0)
        if (k > 1) Pis[,k,i] <- PY[,k-1,i] - PY[,k,i] 
        if (k == l) Pis[,k+1,i] <- PY[,k,i]
      }
      for(j in 1:ID){Y_ord[j,i] <- sample(1:(l+1), size=1, replace=TRUE, prob=sapply(1:(l+1), function(x) Pis[j,x,i])  )} 
      v <- v+1
    }
    
    
    
    #### Tabelle erzeugen & Deskriptive Statistik
    table_ord <- matrix(ncol = m,nrow = ID)

    for (i in 1:m){
      table_ord[,i] <- Y_ord[,i]
    }
    
    
    ############################           Ueberpruefen            ##########################################
    
    output <- sapply(1:(m),function(y){paste0("simuvar",y)})
    table_ord <- as.data.frame(table_ord)
    colnames(table_ord) <-  c(output) #colnames(table_num)<-
    
    
    
    print(paste("Iteration ",g)) #wieviele Iterationen braucht es?
    g=g+1
    
    
    ## Model schätzen
    fit_ord <- tryCatch({lavaan::cfa(model = model, data=table_ord, ordered = output, estimator="WLS", do.fit=T, std.lv=F, control=list(iter.max=1000))},warning = function(w){NA},error = function(e){NA})

    
    ##Ergebnisse Zusammenfügen!
    if( suppressWarnings(!is.na(fit_ord))  ){
      
      if((fitMeasures(fit_ord)["rmsea"] < rmsea_cutoff)){ 
        
        saved_pvalues[[h]] <- fitMeasures(fit_ord)["pvalue"]; names(saved_pvalues)[[h]] <- paste0("pvalue",h)
        saved_rmseas[[h]] <- fitMeasures(fit_ord)["rmsea"]; names(saved_rmseas)[[h]] <- paste0("rmsea",h)
        saved_vars[[h]] <- var; names(saved_vars)[[h]] <- paste0("var",h)
        saved_kappas[[h]] <- kappa; names(saved_kappas)[[h]] <- paste0("kappa",h)
        saved_betas[[h]] <- beta; names(saved_betas)[[h]] <- paste0("beta",h)
        saved_latvars[[h]] <- Psi; names(saved_latvars)[[h]] <- paste0("true_var",h)
        saved_data[[h]] <- table_ord; names(saved_data)[[h]] <- paste0("data",h)
        saved_model[[h]] <- model; names(saved_model)[[h]] <- paste0("model",h)
        
        h=h+1 
        
      } else {next}
    } else {next}
  } 
  
  fits <- list(saved_pvalues,saved_rmseas,saved_vars,saved_kappas,saved_betas,saved_latvars,saved_data,saved_model) 
  names(fits) <- c("pvalues","rmseas","vars","kappas","betas","latvars","data","model")
  return(fits)
}


################################################################################
################################# Sandbox ######################################
################################################################################

#fits_random <- datagen(mode = 'betas', schwellen = 5, ID=1000, times=5, items=8)
