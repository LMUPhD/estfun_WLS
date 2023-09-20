
setwd("C:/Users/classe/Desktop/Diss/Paper3/R")
source("support.R")


# shortcuts
lavdata        <- object@Data
lavmodel       <- object@Model
lavsamplestats <- object@SampleStats
lavoptions     <- object@Options


ntab <- unlist(lavdata@norig)
ntot <- sum(ntab)
npar <- lav_object_inspect_npar(object)
nvar <- ncol(lavsamplestats@cov[[1]])


moments <- fitted(object)
N1 <- 1
X <- lavdata@X[[1]]

Score.mat <- matrix(NA, ntot, npar) #empty matrix



################################################################################
################################# WLS ##########################################
################################################################################

#polychoric corr
polychors = object@SampleStats@cov[[1]]
#s = lavaan::fitted(object)$cov  #doesnt work...

th = lavsamplestats@th[[1]]
#th.p = t(  th %*% t(rep(-1,ntot))  ) #ohne mu_star (weil keine x-Variablen)... 
th.pr = VGAM::probitlink( th*-1,inverse=T) 

#dummies
lv = lavdata@ov[["nlev"]]
Xd = do.call(cbind, lapply(1:nvar, function(i) doDummySingleVar(X,lv,i)  ))


###e1
e1 = t( apply(Xd, 1L, function(x) x-th.pr ) ) 


###e2       
selcols = getCols(lv,nvar)
#catvals = lapply(1:nvar, function(x) as.numeric(unlist( strsplit(lavdata@ov[["lnam"]][x],"|", fixed = T) ))  )
catvals = lapply(1:nvar, function(x)  as.numeric(names(table(X[,x]))) )
mus = unlist(lapply(1:nvar, get_mus))
y_minus_mu = t( apply(X, 1L, function(x) x - mus ) ) 
#im object stecken die falschen Daten drin (1/2 anstatt 0/1)

combs = rbind(  combn(1:nvar,2), lav_matrix_vech(polychors,diagonal=FALSE) ) 
joint_exps = apply(combs, 2L, get_joint_exp) #E(y1y2)
sigma =  joint_exps - t(  lav_matrix_vech(tcrossprod(mus) ,diagonal=FALSE) )  #E(y1y2)-mu1mu2

s_vech = t(apply(y_minus_mu, 1L, function(i){    lav_matrix_vech(tcrossprod(i) ,diagonal=FALSE) })) #s=c( (y1-mu1)(y2-mu2)....

e2 = t( apply(s_vech, 1L, function(x) x - sigma ) ) 



###e
e = cbind(e1,e2)

#weigthing matrix
W = lavsamplestats@WLS.V[[1]] 

#Delta
Delta <- computeDelta(lavmodel = lavmodel)[[1]] #should also work for WLS...

### combine matrices
Score.mat = t( t(Delta) %*% W %*% t(e)  )  #works without mu*...



################################################################################


#Score.mat <- (-1/ntot) * Score.mat #scaling

# provide column names
colnames(Score.mat) <- names(lav_object_inspect_coef(object,
                                                     type = "free", add.labels = TRUE))
