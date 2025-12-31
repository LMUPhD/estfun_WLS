setwd("C:\\Users\\classe\\Desktop\\Diss\\Paper3\\estfun_WLS")



################################################################################
################################################################################

# C implementation #use R 4.4.3!!!!
#toms462
#Rcpp::sourceCpp("Rcpp/toms462_arma.cpp") #compile
#pbivnor(0.5, 0.5, 0.3) #doesnt work...
#R wrapper
#pbivnor_cpp <- function(x, y, r) {
#  pbivnor(x*-1, y*-1, r)
#}

#pbv rcpp
Rcpp::sourceCpp("Rcpp/gee_rcpp_support.cpp")
pbv_rcpp_pbvnorm(0.5, 0.5, 0.3) #that also works...
pbivnorm__rcpp_wls(c(0.5,0.5), c(Inf,0.5), c(0.3,0.3)) #that also works...



# make use of the fact that pbivnor can handle vectors!
# for this let pnorm be added to the c++ file
# then let get_mus and get_joint_exp be written directly in C++





################################################################################
################################################################################



source("application/multi_simu.R") #simulate data (multivariate model)
model_lav = '
  Eta1 =~ simuvar1 + simuvar2 + simuvar3 
  Eta2 =~ simuvar4 + simuvar5 + simuvar6 
  Eta3 =~ simuvar7 + simuvar8 + simuvar9'
model = model_lav
fits_random <- datagen(model=model_lav, schwellen = 4, ID=2000, times=1, items=3,latvar=3) 
Data = fits_random[["data"]][["data1"]]

#estimate parmeters
fit.wls <- lavaan::cfa(model_lav, data = Data, ordered = TRUE, estimator = "WLS", std.lv=F )
object = fit.wls


#from estfun_gee
lavsamplestats <- object@SampleStats
lavdata        <- object@Data
nvar <- ncol(lavsamplestats@cov[[1]])
X <- lavdata@X[[1]]
polychors = object@Fit@Sigma.hat[[1]] 
th = object@Fit@TH[[1]]
th.pr = pnorm( th*-1)                                       #add the diagonal of the model implied matrix!                   
lv = lavdata@ov[["nlev"]]
combs = rbind(  combn(1:nvar,2), lavaan::lav_matrix_vech(polychors,diagonal=FALSE) ) 
catvals = lapply(1:nvar, function(x)  as.integer(names(table(X[,x]))) )



#should be the same result
Rcpp::sourceCpp("Rcpp/gee_rcpp_support.cpp")
system.time(
  apply_get_joint_exp(combs, th, lv, nvar, catvals)
)
print(apply_get_joint_exp(combs, th, lv, nvar, catvals))


source("support.R")
system.time(
  apply(combs, 2L, function(x) get_joint_exp(x, th, lv, nvar, catvals)  ) #E(y1y2)
)
print(apply(combs, 2L, function(x) get_joint_exp(x, th, lv, nvar, catvals)  ) )

















################################################################################
############################### Sandbox ########################################
################################################################################

################################################ redefine function...

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




get_joint_exp <- function(c, th, lv, nvar, catvals){
  
  selcols = getCols(lv,nvar)
  
  #-> Ebene: Item zu Item
  cat_combs = expand.grid(1:lv[c[1]],1:lv[c[2]])
  
  vals_var1 = unlist(catvals[c[1]])
  vals_var2 = unlist(catvals[c[2]])
  
  wth1=selcols[c[1],1]:selcols[c[1],2]
  wth2=selcols[c[2],1]:selcols[c[2],2]
  th_var1 = c(-Inf,th[wth1],Inf)
  th_var2 = c(-Inf,th[wth2],Inf)
  
  #--> Ebene Kategorie-zu-Kategorie
  
  mu_joint = sum( apply(cat_combs, 1L, function(x){
    s = unlist(x+1)
    x = unlist(x)
    p_katkat = sum(pbivnorm__rcpp_wls(x = th_var1[s[1]], y =th_var2[s[2]], rho = c[3]), 
                   pbivnorm__rcpp_wls(x = th_var1[s[1]-1], y =th_var2[s[2]], rho = c[3])*-1,
                   pbivnorm__rcpp_wls(x = th_var1[s[1]], y =th_var2[s[2]-1], rho = c[3])*-1,
                   pbivnorm__rcpp_wls(x = th_var1[s[1]-1], y =th_var2[s[2]-1], rho = c[3]))
    #print(p_katkat)
    vals_var1[x[1]]*vals_var2[x[2]]*p_katkat
  }) )
  
  return(mu_joint)
}

