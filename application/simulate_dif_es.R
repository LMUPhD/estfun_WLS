################################################################################
################ Noncompensatory DIF estimate for Item 2 #######################
########################## Unidimensional Model ################################
################################################################################
setwd("/dss/dsshome1/0B/ra35tik2/paper3")
library(stats)

randparams <- function(N,num.min,rnd=2,num.max=0.99) {
  random.numbers <- round(runif(min=num.min,max=num.max, N),rnd)
  random.signs   <- sample(c(1,-1), N, replace=TRUE)      
  return(random.numbers * random.signs)
}

schwellen_probs <- function(schwellen){
  if(schwellen == 6) return(c(0.10,0.25,0.40,0.60,0.75,0.90)  )
  if(schwellen == 5) return(c(0.10,0.30,0.50,0.70,0.90) )
  if(schwellen == 4) return(c(0.15,0.40,0.60,0.85)  )
  if(schwellen == 3) return(c(0.25,0.5,0.75)  )
  if(schwellen == 2) return(c(0.30,0.70)  )
  if(schwellen == 1) return(c(0.50))
  else stop("Error: Function can only produce data with 2 to 7 categories.")
}


#always five items!
#only DIF estimate for item 2 reported!


ID = 1000;iter=1
means = list()
for(mode in c("all","thresholds","betas")){
  for(l in c(1,2,4,6)){
    
    NCDIFs = list()
    for(h in 1:1000){
      
      scoreks = list()
      psis = list()
      
      if(mode == "thresholds" | mode == "betas" | mode == "none"){
        
        var <- round(runif(min=0.3,max=1.2, 1),2)
        
        if(mode == "thresholds" | mode == "none"){
          beta <- c(1,round(runif(min=0.3,max=1.6, 5-1),2)) #Discrimination parameter
        }
        
        if(mode == "betas" | mode == "none"){
          kappa_shift <- round(rnorm(n = 5, mean=0, sd = 1),2) #"Erwartungswert"-Shift pro Item
          perz_kappa <- matrix(NA,ncol = l, nrow = 5) #Perzentile, um die kappa-Parameter zu bestimmen
          for(i in 1:5){perz_kappa[i,] <- round(schwellen_probs(l) + randparams(num.min=0.01,num.max=0.05, l) ,2)}
        }
      }
      
      for(j in 1:2){
        probitY <- array(NA,c(ID,l)) #(ID, Schwellen, Items pro L-1 Var., L-1 Var.)
        PY <-  array(NA,c(ID,l)) 
        Pis <-  array(NA,c(ID,l+1))
        
        ##### Model Params
        if(mode == "all"){
          var <- round(runif(min=0.3,max=1.2, 1),2)
        }
        if(mode=="betas" | mode == "all"){
          beta <- c(1,round(runif(min=0.3,max=1.6, 5-1),2)) #Discrimination parameter
        }
        if(mode == "thresholds" | mode == "all"){
          kappa_shift <- round(rnorm(n = 5, mean=0, sd = 0.5),2) #"Erwartungswert"-Shift pro Item
          perz_kappa <- matrix(NA,ncol = l, nrow = 5) #Perzentile, um die kappa-Parameter zu bestimmen
          for(i in 1:5){perz_kappa[i,] <- round(schwellen_probs(l) + randparams(num.min=0.01,num.max=0.05, l) ,2)}
        }
        
        
        ##Psi berechnen
        Psi <- rnorm(ID,0,sqrt(var)) + rnorm(ID,mean=0,sd=sqrt(runif(1,min = .3, max = .9))) #psi erstmal mittelwert 0 + Fehlervarianz
        ## Kappa Matrix erstellen
        kappa <- matrix(ncol = 5,nrow = l)  #ncol = L-1 Variablen * Items ; nrow = schwellen
        for(i in 1:5){ # items
          quant <- as.vector(quantile(beta[i]*Psi , probs= perz_kappa[i,])) - kappa_shift[i] 
          kappa[,i] <- quant
        }
        kappa <- round(kappa,2)
        #####
        
        
        for (k in 1:l) { #Schwellen
          probitY[,k] <-  beta[2]*Psi  - kappa[k,2] #Item2!!!
          PY[,k] <- VGAM::probitlink(probitY[,k],inverse=T) 
          if (k == 1) Pis[,k] <- 1 - PY[,k] #P(Y=0)
          if (k > 1) Pis[,k] <- PY[,k-1] - PY[,k] 
          if (k == l) Pis[,k+1] <- PY[,k]
        }
        
        scorek = array(0,ID)
        for (k in 1:(l+1)){
          scorek = scorek + (k-1) * Pis[,k]
        }
        
        psis[[j]] = Psi
        scoreks[[j]] = scorek
      }
      
      ### Score function interpolation!
      app2 = approx(psis[[2]],scoreks[[2]], xout =psis[[1]])
      
      #chalmers_2023 p.322: compute SIDS
      NCDIF = sum(na.exclude(as.vector(scoreks[[1]])-app2$y)^2)/ID
      NCDIFs = append(NCDIFs,NCDIF)
    }
    
    means[[iter]] = mean(unlist(NCDIFs)) #yeah!
    names(means)[[iter]] = paste0("Mode: ",mode,"; k=",l)
    iter = iter+1
    
  }
}

means_univ = means
save(means_univ, file="res/250122_results_NCDIF_univ.RData")







################################################################################
################################ Sandbox #######################################
################################################################################
plot(psis[[1]],scoreks[[1]], col="red")
points(psis[[2]],scoreks[[2]], col="green")

points(app2$x,app2$y, col="green")


#vardiff = fits_random$vars[[1]] - fits_random$vars[[2]]
#betadiff = fits_random$betas[[1]][2] - fits_random$betas[[2]][2]
#kappadiff = fits_random$kappas[[1]][1] - fits_random$kappas[[2]][1]

#vardiffs = append(vardiffs,vardiff)
#betadiffs = append(betadiffs,betadiff)
#kappadiffs = append(kappadiffs,kappadiff)
  

#chalmers_2023 Equation 1 #what about density???
#appdiff = approx(app1$x, app1$y-app2$y, n=ID)
#appdifffun = approxfun(app1$x, app1$y-app2$y)
#ES = integrate(appdifffun, lower = min(appdiff$x), upper = max(appdiff$x))$value


