
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
s = object@SampleStats@cov[[1]]

#thresholed probabilites
yhat = lavPredict(object, type = "yhat")
th = lavsamplestats@th[[1]]
th.pr = sapply(th, function(x) VGAM::probitlink(x*-1,inverse=T)   )

#dummies
lv = lavdata@ov[["nlev"]]
Xd = do.call(cbind, lapply(1:nvar, function(i) doDummySingleVar(X,lv,i)  ))

#item probabilities...
th.idx = lavsamplestats@th.idx[[1]]
thx = cbind(th.idx,th.pr)
s.pr = do.call(cbind, lapply(1:nvar, function(i) getProbItem(X,lv,thx,i)  ))

#sigma_jk (as in Muthen 1997)
th_idx = cbind(th.idx,th)
combs = combn(1:nvar,2)


###!
e1 = t(apply(Xd,1L,  function(x) x - th.pr   )) #may work pretty well without multiplying W...
#e2 <- t(  apply(s.pr, 1L,function(x) lav_matrix_vech((tcrossprod(x) - Sigma.hat),diagonal=FALSE))   )
#e2 = do.call(rbind,  lapply(1:nrow(X),function(i){    lav_matrix_vech(tcrossprod(s.pr[i,]),diagonal=FALSE) - getSigmajk(X,lv,th_idx,s.pr,combs,i)    })  ) #doesnt really work...
e2 = apply(s.pr,1L, function(i){    lav_matrix_vech(tcrossprod(i) - s  ,diagonal=FALSE) }) #doesnt really work...

e = cbind(e1,e2)
e=e*(1/ntot) #????

#weigthing matrix
W = lavsamplestats@WLS.V[[1]] #do we really need this?


Delta <- computeDelta(lavmodel = lavmodel)[[1]] #should also work for WLS...
scores.H1 = e %*% W   #means are close to zero but not sums....

wi <- lavdata@case.idx[[1]]
Score.mat[wi,] <- -scores.H1 %*% Delta #Zwischenfazit (21.08.): Threshold derivatives eigentlich gar nicht so schlecht... Probleme mit allen anderen Parametern...


################################################################################




Score.mat[wi,] <- (-1/ntot) * Score.mat[wi,] #scaling


# provide column names
colnames(Score.mat) <- names(lav_object_inspect_coef(object,
                                                     type = "free", add.labels = TRUE))
