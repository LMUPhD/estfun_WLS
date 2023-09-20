
setwd("C:/Users/classe/Desktop/Diss/Paper3/R")
load("liss_univ_data.RData") 


## mirt:
library(mirt)
dt <- mirt::expand.table(LSAT7)
#fit_mirt <- mirt(dt[,1:5], 1, itemtype = "Rasch", SE = TRUE)
#mirt_scores <- estfun.AllModelClass(fit_mirt)

####################
library(lavaan)

model ='Eta1 =~ cp08a014 + cp08a015 + cp08a016 + cp08a017 + cp08a018'
model = 'Eta1 =~ 1*Item.1 + 1*Item.2 + 1*Item.3 + 1*Item.4 + 1*Item.5'
#fit <- cfa(model, data = dt, meanstructure = TRUE)   #fit with ML
#sc <- estfun.lavaan(fit)                      #calculate scores for model fit with ML


fit_ord <- cfa(model, data = dt, ordered = TRUE, estimator = "WLS" ) #parameterization = "theta"
object = fit_ord







################################################################################
############################## Sandbox ######################################### 
################################################################################





getCols <- function(lv,nvar){
  maxcols = mincols = c()
  maxcol = mincol = 0
  for(i in 1:nvar){
    mincol = maxcol+1
    maxcol = maxcol + lv[i]-1
    mincols = c(mincols,mincol)
    maxcols = c(maxcols,maxcol)
  }
  return(cbind(mincols,maxcols))
}


get_probit_categories <- function(var){  #, th, lv, selcols
  
  cols = selcols[var,1]:selcols[var,2]
  p.item = th[cols] 
  thps.item = sapply(1:lv[var], function(x){
    if(x==1) {
      probit_cat = p.item[1]
    }
    else if(x==lv[var]) {
      probit_cat = tail(p.item, n=1)*-1
    } else{ 
      probit_cat = VGAM::probitlink( VGAM::probitlink(p.item[x], inverse = T) - 
                                       VGAM::probitlink(p.item[x-1], inverse = T)   )
    }
    return(probit_cat)
  }) 
  return(thps.item)
}



get_mu <- function(probit_categories,catvals,lv){  
  
  sapply(1:length(probit_categories), function(x){
    sum(sapply(1:lv[x], function(y){
      catvals[[x]][y] * VGAM::probitlink(probit_categories[[x]][y] ,inverse=T)
    }))
  })
  
}


get_joint_exp <- function(c){  #catvals, probit_categories
  
  #-> Ebene: Item zu Item
  cat_combs = expand.grid(1:lv[c[1]],1:lv[c[2]])
  
  vals_var1 = unlist(catvals[c[1]])
  vals_var2 = unlist(catvals[c[2]])
  
  probit_var1 = unlist(probit_categories[c[1]])
  probit_var2 = unlist(probit_categories[c[2]])
  
  #--> Ebene Kategorie-zu-Kategorie
  
  mu_joint = sum( apply(cat_combs, 1L, function(x){
    x = unlist(x)
    p_katkat = pbivnorm::pbivnorm(x = probit_var1[x[1]], y = probit_var2[x[2]], rho = c[3], recycle = TRUE)
    vals_var1[x[1]]*vals_var2[x[2]]*p_katkat
  }) )
  
  return(mu_joint)
}


####

selcols = getCols(lv,nvar)
catvals = lapply(1:nvar, function(x) as.numeric(unlist( strsplit(lavdata@ov[["lnam"]][x],"|", fixed = T) ))  )
probit_categories = lapply(1:nvar, get_probit_categories)
mus = get_mu(probit_categories,catvals,lv)
y_minus_mu = t( apply(X, 1L, function(x) x - mus ) ) 


combs = rbind(  combn(1:nvar,2), lav_matrix_vech(polychors,diagonal=FALSE) ) 
joint_exps = apply(combs, 2L, get_joint_exp) #E(y1y2)
sigma =  joint_exps - t(  lav_matrix_vech(tcrossprod(mus) ,diagonal=FALSE) )

s_vech = t(apply(y_minus_mu, 1L, function(i){    lav_matrix_vech(tcrossprod(i) ,diagonal=FALSE) })) #s=c( (y1-mu1)(y2-mu2)....

e2 = t( apply(s_vech, 1L, function(x) x - sigma ) ) 
####









#sigma = apply(combs,2L,function(x) pbivnorm::pbivnorm(x = th.p[,x[1]], y = th.p[,x[2]], rho = x[3], recycle = TRUE)  ) #E(y1y2)






