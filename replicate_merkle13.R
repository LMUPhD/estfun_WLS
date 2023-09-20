
setwd("C:/Users/classe/Desktop/Diss/Paper3/R")
load("liss_univ_data.RData") 


## mirt:
library(mirt)
dt <- mirt::expand.table(LSAT7)
fit_mirt <- mirt(dt[,1:5], 1, itemtype = "Rasch", SE = TRUE)
mirt_scores <- estfun.AllModelClass(fit_mirt)

####################
library(lavaan)

#model ='Eta1 =~ cp08a014 + cp08a015 + cp08a016 + cp08a017 + cp08a018'
#fit <- cfa(model, data = dt, meanstructure = TRUE)   #fit with ML
#sc <- estfun.lavaan(fit)                      #calculate scores for model fit with ML



model = 'Eta1 =~ 1*Item.1 + 1*Item.2 + 1*Item.3 + 1*Item.4 + 1*Item.5'

fit_ord <- cfa(model, data = dt, ordered = TRUE, estimator = "WLS" ) #parameterization = "theta"
object = fit_ord



#Daten besorgen mit drei Antwortkategorien!


################################################################################
############################## Sandbox ######################################### 
################################################################################


item_probs1 = unlist(mus[[1]][2]) #variable 1
item_probs2 = unlist(mus[[2]][2]) #variable 2
p1 = VGAM::probitlink(item_probs1[2])
p2 = VGAM::probitlink(item_probs2[1])

pbivnorm::pbivnorm(x = p1, y = p2, rho = 0.2265593, recycle = TRUE) 

th1 = -0.9462914 
th2 = -0.4070109


#P(y1=0, y2=0) = 0.081
pbivnorm::pbivnorm(x =   th1, y = th2, rho = 0.2265593, recycle = TRUE) #alle anderen müssen Null sein...
pbivnorm::pbivnorm(x =  -Inf, y = th2, rho = 0.2265593, recycle = TRUE)
pbivnorm::pbivnorm(x =  th1, y = -Inf, rho = 0.2265593, recycle = TRUE)
pbivnorm::pbivnorm(x = -Inf, y = -Inf, rho = 0.2265593, recycle = TRUE) #-Inf, -Inf muss null sein!
#müsste passen...



#P(y1=1, y2=1) = 0.567
pbivnorm::pbivnorm(x = Inf, y = Inf, rho = 0.2265593, recycle = TRUE) #Inf, Inf muss eins sein!
pbivnorm::pbivnorm(x = th1, y = Inf, rho = 0.2265593, recycle = TRUE)
pbivnorm::pbivnorm(x = Inf, y = th2, rho = 0.2265593, recycle = TRUE)
pbivnorm::pbivnorm(x = th1, y = th2, rho = 0.2265593, recycle = TRUE)


#P(y1=0, y2=1) = 0.1334526 ???WTF
pbivnorm::pbivnorm(x = th1,  y = Inf,   rho = 0.2265593, recycle = TRUE) 
pbivnorm::pbivnorm(x = -Inf, y = Inf,   rho = 0.2265593, recycle = TRUE)
pbivnorm::pbivnorm(x = th1,  y = th2,   rho = 0.2265593, recycle = TRUE)
pbivnorm::pbivnorm(x = -Inf, y =  th2,  rho = 0.2265593, recycle = TRUE)


#P(y1=1, y2=0) = 0.3034526 ???WTF
pbivnorm::pbivnorm(x = Inf,  y = th2,  rho = 0.2265593, recycle = TRUE) 
pbivnorm::pbivnorm(x = th1,  y = th2,  rho = 0.2265593, recycle = TRUE)
pbivnorm::pbivnorm(x = Inf,  y = -Inf, rho = 0.2265593, recycle = TRUE)
pbivnorm::pbivnorm(x = th1,  y = -Inf, rho = 0.2265593, recycle = TRUE)


pbivnorm_wls <- function(x,y,rho){
  if(x==Inf & y==Inf){return(1)}
  if(x==-Inf & y==-Inf){return(0)}
  else {return(  pbivnorm::pbivnorm(x = x, y = y, rho = rho, recycle = TRUE)   )} 
}


sum(c(
  pbivnorm_wls(x = Inf,  y = th2,  rho = 0.2265593),
  pbivnorm_wls(x = th1,  y = th2,  rho = 0.2265593)*-1,
  pbivnorm_wls(x = Inf,  y = -Inf, rho = 0.2265593)*-1,
  pbivnorm_wls(x = th1,  y = -Inf, rho = 0.2265593)
))

sum(0.423,0.091,0.261,0.567) #doesnt add up to one...

###########


get_joint_exp <- function(c){  #catvals, th, lv
  
  #-> Ebene: Item zu Item
  cat_combs = expand.grid(1:lv[c[1]],1:lv[c[2]])
  
  vals_var1 = unlist(catvals[c[1]])
  vals_var2 = unlist(catvals[c[2]])
  
  th_var1 = c(-Inf,th[c[1]],Inf) #*-1
  th_var2 = c(-Inf,th[c[2]],Inf) #*-1
  
  #--> Ebene Kategorie-zu-Kategorie
  
  mu_joint = sum( apply(cat_combs, 1L, function(x){
    s = unlist(x+1)
    x=unlist(x)
    p_katkat = sum(c(pbivnorm::pbivnorm(x = th_var1[s[1]], y =th_var2[s[2]], rho = -c[3], recycle = TRUE), 
                     pbivnorm::pbivnorm(x = th_var1[s[1]-1], y =th_var2[s[2]], rho = c[3], recycle = TRUE)*-1,
                     pbivnorm::pbivnorm(x = th_var1[s[1]], y =th_var2[s[2]-1], rho = c[3], recycle = TRUE)*-1,
                     pbivnorm::pbivnorm(x = th_var1[s[1]-1], y =th_var2[s[2]-1], rho = c[3], recycle = TRUE)),na.rm=T)
    vals_var1[x[1]]*vals_var2[x[2]]*p_katkat
  }) )
  
  return(mu_joint)
}




