
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


moments <- lavaan::fitted(object)
N1 <- 1
X <- lavdata@X[[1]]

Score.mat <- matrix(NA, ntot, npar) #empty matrix



################################################################################
################################# WLS ##########################################
################################################################################

#polychoric corr
s = object@SampleStats@cov[[1]]

th = lavsamplestats@th[[1]]
th.p = t(  th %*% t(rep(-1,ntot))  ) #ohne mu_star... geht das?
th.pr = VGAM::probitlink( th.p,inverse=T) 

#dummies
lv = lavdata@ov[["nlev"]]
Xd = do.call(cbind, lapply(1:nvar, function(i) doDummySingleVar(X,lv,i)  ))

###e1
e1 = Xd-th.pr

###e2        #adjust for ordinal!
combs = rbind(  combn(1:nvar,2), lav_matrix_vech(s,diagonal=FALSE) ) 
sigma = apply(combs,2L,function(x) pbivnorm::pbivnorm(x = th.p[,x[1]], y = th.p[,x[2]], rho = x[3], recycle = TRUE)  ) #E(y1y2)
sigma = sigma - t(apply(th.pr,1L, function(i){    lav_matrix_vech(tcrossprod(i) ,diagonal=FALSE) })) #sigma = E(y1y2)-mu1mu2
s_vech = t(apply(e1,1L, function(i){    lav_matrix_vech(tcrossprod(i) ,diagonal=FALSE) })) #s=c( (y1-mu1)(y2-mu2)....
e2 = s_vech - sigma

###e
e = cbind(e1,e2)

#weigthing matrix
W = lavsamplestats@WLS.V[[1]] 

#Delta
Delta <- computeDelta(lavmodel = lavmodel)[[1]] #should also work for WLS...

#Compute score matrix
Score.mat = t( t(Delta) %*% W %*% t(e)  )  


################################################################################




#Score.mat <- (-1/ntot) * Score.mat #scaling


# provide column names
colnames(Score.mat) <- names(lav_object_inspect_coef(object,
                                                     type = "free", add.labels = TRUE))


