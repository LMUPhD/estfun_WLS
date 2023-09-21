
################################ my code



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


get_mus <- function(var, th, lv, nvar,  catvals){  
  
  selcols = getCols(lv,nvar)
  cols = selcols[var,1]:selcols[var,2]
  p.item = th[cols] 
  catprobs.item = sapply(1:lv[var], function(x){
    if(x==1) {
      prob_cat = VGAM::probitlink(p.item[1], inverse=T)
    }
    else if(x==lv[var]) {
      prob_cat = VGAM::probitlink(tail(p.item, n=1)*-1, inverse=T)
    } else{ 
      prob_cat = VGAM::probitlink(p.item[x], inverse = T) - 
        VGAM::probitlink(p.item[x-1], inverse = T)   
    }
    return(prob_cat)
  }) 
  
  
  
  mu = sum(sapply(1:lv[var], function(y){
    catvals[[var]][y] * catprobs.item[y]
  }))
  
  
  return(mu)
}


pbivnorm_wls <- function(x,y,rho){
  if(x==Inf & y==Inf){return(1)}
  else if(x==-Inf | y==-Inf){return(0)}
  else if(x==Inf & y>1){return(   pbivnorm::pbivnorm(x = Inf, y = 1, rho = rho, recycle = TRUE)    )}
  else if(x>1 & y==Inf){return(   pbivnorm::pbivnorm(x = 1, y = Inf, rho = rho, recycle = TRUE)    )}
  else {return(  pbivnorm::pbivnorm(x = x, y = y, rho = rho, recycle = TRUE)   )} 
}


get_joint_exp <- function(c, X, th, lv, nvar, catvals){
  
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
    p_katkat = sum(pbivnorm_wls(x = th_var1[s[1]], y =th_var2[s[2]], rho = c[3]), 
                   pbivnorm_wls(x = th_var1[s[1]-1], y =th_var2[s[2]], rho = c[3])*-1,
                   pbivnorm_wls(x = th_var1[s[1]], y =th_var2[s[2]-1], rho = c[3])*-1,
                   pbivnorm_wls(x = th_var1[s[1]-1], y =th_var2[s[2]-1], rho = c[3]))
    print(p_katkat)
    vals_var1[x[1]]*vals_var2[x[2]]*p_katkat
  }) )
  
  return(mu_joint)
}




doDummySingleVar <- function(X,lv,ntot,num){
  Xd = matrix(NA,nrow = ntot, ncol =  lv[num]-1 ) 
  x = X[,num]
  minx = min(x)
  categ=minx-1
  v=1
  while(categ < lv[num]-1){
    categ = categ+1
    Xd[,v] = ifelse(x > categ, 1, 0) 
    v=v+1
  }
  return(Xd)
}



################################ lavaan code
lav_object_inspect_npar <- function(object, type = "free") {
  
  if(type == "free") {
    npar <- sum(object@ParTable$free > 0L &
                  !duplicated(object@ParTable$free))
  } else {
    npar <- length(object@ParTable$lhs)
  }
  
  npar
}

inv.chol <- function(S, logdet=FALSE) {
  cS <- chol(S)
  #if( inherits(cS, "try-error") ) {
  #    print(S)
  #    warning("lavaan WARNING: symmetric matrix is not positive symmetric!")
  #}
  S.inv <- chol2inv( cS )
  if(logdet) {
    diag.cS <- diag(cS)
    attr(S.inv, "logdet") <- sum(log(diag.cS*diag.cS))
  }
  S.inv
}

computeDelta <- function(lavmodel = NULL, GLIST. = NULL,
                         m.el.idx. = NULL, x.el.idx. = NULL,
                         ceq.simple = FALSE,
                         force.conditional.x.false = FALSE) {
  
  representation   <- lavmodel@representation
  categorical      <- lavmodel@categorical
  if(.hasSlot(lavmodel, "correlation")) {
    correlation   <- lavmodel@correlation
  } else {
    correlation   <- FALSE
  }
  conditional.x    <- lavmodel@conditional.x
  group.w.free     <- lavmodel@group.w.free
  nmat             <- lavmodel@nmat
  nblocks          <- lavmodel@nblocks
  nvar             <- lavmodel@nvar
  num.idx          <- lavmodel@num.idx
  th.idx           <- lavmodel@th.idx
  nexo             <- lavmodel@nexo
  parameterization <- lavmodel@parameterization
  
  # number of thresholds per group (if any)
  nth <- sapply(th.idx, function(x) sum(x > 0L))
  
  # state or final?
  if(is.null(GLIST.))
    GLIST <- lavmodel@GLIST
  else
    GLIST <- GLIST.
  
  # type = "free" or something else?
  type <- "nonfree"
  m.el.idx <- m.el.idx.; x.el.idx <- x.el.idx.
  if(is.null(m.el.idx) && is.null(x.el.idx))
    type <- "free"
  
  # number of rows in DELTA.group
  pstar <- integer(nblocks)
  for(g in 1:nblocks) {
    pstar[g] <- as.integer(nvar[g] * (nvar[g] + 1) / 2)
    if(lavmodel@meanstructure) {
      pstar[g] <- nvar[g] + pstar[g]  # first the means, then sigma
    }
    if(categorical) {
      pstar[g] <- pstar[g] - nvar[g] # remove variances
      pstar[g] <- pstar[g] - nvar[g] # remove means
      
      pstar[g] <- pstar[g] + nth[g]  # add thresholds
      pstar[g] <- pstar[g] + length(num.idx[[g]]) # add num means
      pstar[g] <- pstar[g] + length(num.idx[[g]]) # add num vars
    } else if(correlation) {
      pstar[g] <- pstar[g] - nvar[g] # remove variances
    }
    if(conditional.x && nexo[g] > 0L) {
      pstar[g] <- pstar[g] + (nvar[g] * nexo[g]) # add slopes
    }
    if(group.w.free) {
      pstar[g] <- pstar[g] + 1L # add group weight
    }
  }
  
  
  # number of columns in DELTA + m.el.idx/x.el.idx
  if(type == "free") {
    if(.hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      NCOL <- lavmodel@nx.unco
    } else {
      NCOL <- lavmodel@nx.free
    }
    m.el.idx <- x.el.idx <- vector("list", length=length(GLIST))
    for(mm in 1:length(GLIST)) {
      m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
      if(.hasSlot(lavmodel, "ceq.simple.only") &&
         lavmodel@ceq.simple.only) {
        x.el.idx[[mm]] <- lavmodel@x.unco.idx[[mm]]
      } else {
        x.el.idx[[mm]] <- lavmodel@x.free.idx[[mm]]
      }
      # handle symmetric matrices
      if(lavmodel@isSymmetric[mm]) {
        # since we use 'x.free.idx', only symmetric elements
        # are duplicated (not the equal ones, only in x.free.free)
        dix <- duplicated(x.el.idx[[mm]])
        if(any(dix)) {
          m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
          x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
        }
      }
    }
  } else {
    ## FIXME: this does *not* take into account symmetric
    ##        matrices; hence NCOL will be too large, and empty
    ##        columns will be added
    ##        this is ugly, but it doesn't hurt
    ## alternative could be:
    ## NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
    #NCOL <- sum(unlist(lapply(m.el.idx, length)))
    NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
    # sanity check
    #nx <- sum(unlist(lapply(x.el.idx, length)))
    #stopifnot(NCOL == nx)
  }
  
  
  # compute Delta
  Delta <- vector("list", length=nblocks)
  for(g in 1:nblocks) {
    Delta.group <- matrix(0, nrow=pstar[g], ncol=NCOL)
    
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    
    # label rows of Delta.group --- FIXME!!!
    #if(categorical) {
    #    # 1. th (means interleaved?)
    #    # 2. pi
    #    # 3. var num + cor
    #} else {
    #    if(meanstructure) {
    #    }
    #}
    #if(group.w.free) {
    #}
    
    # if theta, do some preparation
    if(representation == "LISREL" && parameterization == "theta") {
      sigma.hat <- computeSigmaHat.LISREL(MLIST=GLIST[mm.in.group],
                                          delta=FALSE)
      dsigma <- diag(sigma.hat)
      # dcor/dcov for sigma
      R <- lav_deriv_cov2cor(sigma.hat, num.idx = lavmodel@num.idx[[g]])
      theta.var.idx <- lav_matrix_diagh_idx(nvar[g])
    }
    
    for(mm in mm.in.group) {
      mname <- names(lavmodel@GLIST)[mm]
      
      # skip empty ones
      if(!length(m.el.idx[[mm]])) next
      
      # get Delta columns for this model matrix
      if(representation == "LISREL") {
        
        # Sigma
        DELTA <- dxSigma <-
          derivative.sigma.LISREL(m = mname,
                                  idx = m.el.idx[[mm]],
                                  MLIST = GLIST[ mm.in.group ],
                                  delta = parameterization == "delta")
        if(categorical && parameterization == "theta") {
          DELTA <- R %*% DELTA
        }
        
        if(categorical) {
          # reorder: first variances (of numeric), then covariances
          cov.idx  <- lav_matrix_vech_idx(nvar[g])
          covd.idx <- lav_matrix_vech_idx(nvar[g], diagonal = FALSE)
          
          var.idx <- which(is.na(match(cov.idx,
                                       covd.idx)))[num.idx[[g]]]
          cor.idx <- match(covd.idx, cov.idx)
          
          DELTA <- rbind(DELTA[var.idx,,drop=FALSE],
                         DELTA[cor.idx,,drop=FALSE])
        }
        
        # correlation structure?
        if(!categorical && correlation) {
          rm.idx <- lav_matrix_diagh_idx(nvar[g])
          DELTA <- DELTA[-rm.idx, , drop = FALSE]
        }
        
        if(!categorical) {
          if(conditional.x) {
            # means/intercepts
            DELTA.mu <- derivative.mu.LISREL(m=mname,
                                             idx=m.el.idx[[mm]], MLIST=GLIST[ mm.in.group ])
            
            # slopes
            if(lavmodel@nexo[g] > 0L) {
              DELTA.pi <- derivative.pi.LISREL(m=mname,
                                               idx=m.el.idx[[mm]], MLIST=GLIST[ mm.in.group ])
              
              if(lavmodel@multilevel) {
                DELTA <- rbind(DELTA.mu, DELTA.pi, DELTA)
              } else {
                # ATTENTION: we need to change the order here
                # lav_mvreg_scores_* uses 'Beta' where the
                # the intercepts are just the first row
                # using the col-major approach, we need to
                # interweave the intercepts with the slopes!
                
                nEls <- NROW(DELTA.mu) + NROW(DELTA.pi)
                # = (nexo + 1 int) * nvar
                
                # intercepts on top
                tmp <- rbind(DELTA.mu, DELTA.pi)
                # change row index
                row.idx <- lav_matrix_vec(matrix(seq.int(nEls),
                                                 nrow = lavmodel@nexo[g] + 1L,
                                                 ncol = lavmodel@nvar[g], byrow = TRUE))
                DELTA.beta <- tmp[row.idx,,drop = FALSE]
                DELTA <- rbind(DELTA.beta, DELTA)
              }
            } else {
              DELTA <- rbind(DELTA.mu, DELTA)
            }
          } else if(!conditional.x && lavmodel@meanstructure) {
            DELTA.mu <- derivative.mu.LISREL(m=mname,
                                             idx=m.el.idx[[mm]], MLIST=GLIST[ mm.in.group ])
            DELTA <- rbind(DELTA.mu, DELTA)
          }
        }
        
        else if(categorical) {
          DELTA.th <- derivative.th.LISREL(m=mname,
                                           idx=m.el.idx[[mm]],
                                           th.idx=th.idx[[g]],
                                           MLIST=GLIST[ mm.in.group ],
                                           delta = TRUE)
          if(parameterization == "theta") {
            # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
            dDelta.dx <-
              ( dxSigma[theta.var.idx,,drop=FALSE] *
                  -0.5 / (dsigma*sqrt(dsigma)) )
            dth.dDelta <-
              derivative.th.LISREL(m = "delta",
                                   idx = 1:nvar[g],
                                   MLIST = GLIST[ mm.in.group ],
                                   th.idx = th.idx[[g]])
            # add dth.dDelta %*% dDelta.dx
            no.num.idx <- which(th.idx[[g]] > 0)
            DELTA.th[no.num.idx,] <-
              DELTA.th[no.num.idx,,drop=FALSE] +
              (dth.dDelta %*% dDelta.dx)[no.num.idx,,drop=FALSE]
          }
          if(conditional.x && lavmodel@nexo[g] > 0L) {
            DELTA.pi <-
              derivative.pi.LISREL(m=mname,
                                   idx=m.el.idx[[mm]],
                                   MLIST=GLIST[ mm.in.group ])
            if(parameterization == "theta") {
              dpi.dDelta <-
                derivative.pi.LISREL(m = "delta",
                                     idx = 1:nvar[g],
                                     MLIST = GLIST[ mm.in.group ])
              # add dpi.dDelta %*% dDelta.dx
              no.num.idx <-
                which(!seq.int(1L,nvar[g]) %in% num.idx[[g]])
              no.num.idx <- rep(seq.int(0,nexo[g]-1) * nvar[g],
                                each=length(no.num.idx)) + no.num.idx
              DELTA.pi[no.num.idx,] <-
                DELTA.pi[no.num.idx,,drop=FALSE] +
                (dpi.dDelta %*% dDelta.dx)[no.num.idx,,drop=FALSE]
            }
            DELTA <- rbind(DELTA.th, DELTA.pi, DELTA)
          } else {
            DELTA <- rbind(DELTA.th, DELTA)
          }
        }
        if(group.w.free) {
          DELTA.gw <- derivative.gw.LISREL(m=mname,
                                           idx=m.el.idx[[mm]],
                                           MLIST=GLIST[ mm.in.group ])
          DELTA <- rbind(DELTA.gw, DELTA)
        }
      } else if(representation == "RAM") {
        DELTA <- dxSigma <-
          lav_ram_dsigma(m     = mname,
                         idx   = m.el.idx[[mm]],
                         MLIST = GLIST[ mm.in.group ])
        if(lavmodel@meanstructure) {
          DELTA.mu <- lav_ram_dmu(m    = mname,
                                  idx   = m.el.idx[[mm]],
                                  MLIST = GLIST[ mm.in.group ])
          DELTA <- rbind(DELTA.mu, DELTA)
        }
      } else {
        stop("representation ", representation, " not implemented yet")
      }
      
      Delta.group[ ,x.el.idx[[mm]]] <- DELTA
    } # mm
    
    # if type == "free" take care of equality constraints
    if(type == "free" && ceq.simple &&
       .hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      Delta.group <- Delta.group %*% lavmodel@ceq.simple.K
    }
    
    Delta[[g]] <- Delta.group
    
  } # g
  
  # if multilevel, rbind levels within group
  if(.hasSlot(lavmodel, "multilevel") && lavmodel@multilevel) {
    DELTA <- vector("list", length = lavmodel@ngroups)
    for(g in 1:lavmodel@ngroups) {
      DELTA[[g]] <- rbind( Delta[[(g-1)*2 + 1]],
                           Delta[[(g-1)*2 + 2]] )
    }
    Delta <- DELTA
  }
  
  Delta
}

# dSigma/dx -- per model matrix
derivative.sigma.LISREL <- function(m     = "lambda",
                                    # all model matrix elements, or only a few?
                                    # NOTE: for symmetric matrices,
                                    # we assume that the have full size
                                    # (nvar*nvar) (but already correct for
                                    # symmetry)
                                    idx   = seq_len(length(MLIST[[m]])),
                                    MLIST = NULL,
                                    vech  = TRUE,
                                    delta = TRUE) {
  
  LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
  PSI    <- MLIST$psi
  
  # only lower.tri part of sigma (not same order as elimination matrix?)
  v.idx <- lav_matrix_vech_idx( nvar ); pstar <- nvar*(nvar+1)/2
  
  # shortcut for gamma, nu, alpha, tau,.... : empty matrix
  if(m == "nu" || m == "alpha" || m == "tau" || m == "gamma" ||
     m == "gw" || m == "cov.x" || m == "mean.x") {
    return( matrix(0.0, nrow=pstar, ncol=length(idx)) )
  }
  
  # Delta?
  delta.flag <- FALSE
  if(delta && !is.null(MLIST$delta)) {
    DELTA <- MLIST$delta
    delta.flag <- TRUE
  } else if(m == "delta") { # modindices?
    return( matrix(0.0, nrow=pstar, ncol=length(idx)) )
  }
  
  # beta?
  if(!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
  }
  
  # pre
  #if(m == "lambda" || m == "beta")
  #    IK <- diag(nvar*nvar) + lav_matrix_commutation(nvar, nvar)
  if(m == "lambda" || m == "beta") {
    L1 <- LAMBDA %*% IB.inv %*% PSI %*% t(IB.inv)
  }
  if(m == "beta" || m == "psi") {
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }
  
  # here we go:
  if(m == "lambda") {
    KOL.idx <- matrix(1:(nvar*nfac), nvar, nfac, byrow = TRUE)[idx]
    DX <- (L1 %x% diag(nvar))[,idx, drop = FALSE] +
      (diag(nvar) %x% L1)[,KOL.idx, drop = FALSE]
  } else if(m == "beta") {
    KOL.idx <- matrix(1:(nfac*nfac), nfac, nfac, byrow = TRUE)[idx]
    DX <- (L1 %x% LAMBDA..IB.inv)[,idx, drop = FALSE] +
      (LAMBDA..IB.inv %x% L1)[, KOL.idx, drop = FALSE]
    # this is not really needed (because we select idx=m.el.idx)
    # but just in case we need all elements of beta...
    DX[, which(idx %in% lav_matrix_diag_idx(nfac))] <- 0.0
  } else if(m == "psi") {
    DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv)
    # symmetry correction, but keeping all duplicated elements
    # since we depend on idx=m.el.idx
    lower.idx <- lav_matrix_vech_idx(nfac, diagonal = FALSE)
    upper.idx <- lav_matrix_vechru_idx(nfac, diagonal = FALSE)
    offdiagSum <- DX[,lower.idx] + DX[,upper.idx]
    DX[,c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
    DX <- DX[,idx, drop = FALSE]
  } else if(m == "theta") {
    #DX <- diag(nvar*nvar) # very sparse...
    DX <- matrix(0, nvar*nvar, length(idx))
    DX[cbind(idx,seq_along(idx))] <- 1
    # symmetry correction not needed, since all off-diagonal elements
    # are zero?
  } else if(m == "delta") {
    Omega <- computeSigmaHat.LISREL(MLIST, delta=FALSE)
    DD <- diag(DELTA[,1], nvar, nvar)
    DD.Omega <- (DD %*% Omega)
    A <- DD.Omega %x% diag(nvar); B <- diag(nvar) %x% DD.Omega
    DX <- A[,lav_matrix_diag_idx(nvar),drop=FALSE] +
      B[,lav_matrix_diag_idx(nvar),drop=FALSE]
    DX <- DX[,idx, drop = FALSE]
  } else {
    stop("wrong model matrix names: ", m, "\n")
  }
  
  if(delta.flag && !m == "delta") {
    DX <- DX * as.vector(DELTA %x% DELTA)
  }
  
  # vech?
  if(vech) {
    DX <- DX[v.idx,, drop=FALSE]
  }
  
  DX
}

.internal_get_IB.inv <- function(MLIST = NULL) {
  
  BETA <- MLIST$beta; nr <- nrow(MLIST$psi)
  
  if(!is.null(BETA)) {
    tmp <- -BETA
    tmp[lav_matrix_diag_idx(nr)] <- 1
    IB.inv <- solve(tmp)
  } else {
    IB.inv <- diag(nr)
  }
  
  IB.inv
}




derivative.mu.LISREL <- function(m="alpha",
                                 # all model matrix elements, or only a few?
                                 idx=seq_len(length(MLIST[[m]])),
                                 MLIST=NULL) {
  
  
  LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
  
  # shortcut for empty matrices
  if(m == "gamma" || m == "psi" || m == "theta" || m == "tau" ||
     m == "delta"|| m == "gw" || m == "cov.x" || m == "mean.x") {
    return( matrix(0.0, nrow=nvar, ncol=length(idx) ) )
  }
  
  # missing alpha
  if(is.null(MLIST$alpha))
    ALPHA <- matrix(0, nfac, 1L)
  else
    ALPHA  <- MLIST$alpha
  
  
  # beta?
  if(!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
  }
  
  if(m == "nu") {
    DX <- diag(nvar)
  } else if(m == "lambda") {
    DX <- t(IB.inv %*% ALPHA) %x% diag(nvar)
  } else if(m == "beta") {
    DX <- t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    DX[,lav_matrix_diag_idx(nfac)] <- 0.0
  } else if(m == "alpha") {
    DX <- LAMBDA %*% IB.inv
  } else {
    stop("wrong model matrix names: ", m, "\n")
  }
  
  DX <- DX[, idx, drop=FALSE]
  DX
}




derivative.th.LISREL <- function(m="tau",
                                 # all model matrix elements, or only a few?
                                 idx=seq_len(length(MLIST[[m]])),
                                 th.idx=NULL,
                                 MLIST=NULL,
                                 delta = TRUE) {
  
  
  LAMBDA <- MLIST$lambda; nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)
  TAU <- MLIST$tau; nth <- nrow(TAU)
  
  # missing alpha
  if(is.null(MLIST$alpha)) {
    ALPHA <- matrix(0, nfac, 1L)
  } else {
    ALPHA  <- MLIST$alpha
  }
  
  # missing nu
  if(is.null(MLIST$nu)) {
    NU <- matrix(0, nvar, 1L)
  } else {
    NU <- MLIST$nu
  }
  
  # Delta?
  delta.flag <- FALSE
  if(delta && !is.null(MLIST$delta)) {
    DELTA <- MLIST$delta
    delta.flag <- TRUE
  }
  
  if(is.null(th.idx)) {
    th.idx <- seq_len(nth)
    nlev <- rep(1L, nvar)
    K_nu <- diag(nvar)
  } else {
    nlev <- tabulate(th.idx, nbins=nvar); nlev[nlev == 0L] <- 1L
    K_nu <- matrix(0, sum(nlev), nvar)
    K_nu[ cbind(seq_len(sum(nlev)), rep(seq_len(nvar), times=nlev)) ] <- 1.0
  }
  
  # shortcut for empty matrices
  if(m == "gamma" || m == "psi" || m == "theta" || m == "gw" ||
     m == "cov.x" || m == "mean.x") {
    return( matrix(0.0, nrow=length(th.idx), ncol=length(idx) ) )
  }
  
  # beta?
  if(!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
  }
  
  if(m == "tau") {
    DX <- matrix(0, nrow=length(th.idx), ncol=nth)
    DX[ th.idx > 0L, ] <-  diag(nth)
    if(delta.flag)
      DX <- DX * as.vector(K_nu %*% DELTA)
  } else if(m == "nu") {
    DX <- (-1) * K_nu
    if(delta.flag)
      DX <- DX * as.vector(K_nu %*% DELTA)
  } else if(m == "lambda") {
    DX <- (-1) * t(IB.inv %*% ALPHA) %x% diag(nvar)
    DX <- K_nu %*% DX
    if(delta.flag)
      DX <- DX * as.vector(K_nu %*% DELTA)
  } else if(m == "beta") {
    DX <- (-1) * t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    DX[,lav_matrix_diag_idx(nfac)] <- 0.0
    DX <- K_nu %*% DX
    if(delta.flag)
      DX <- DX * as.vector(K_nu %*% DELTA)
  } else if(m == "alpha") {
    DX <- (-1) * LAMBDA %*% IB.inv
    DX <- K_nu %*% DX
    if(delta.flag)
      DX <- DX * as.vector(K_nu %*% DELTA)
  } else if(m == "delta") {
    DX1 <- matrix(0, nrow=length(th.idx), ncol=1)
    DX1[ th.idx > 0L, ] <-  TAU
    DX2 <- NU + LAMBDA %*% IB.inv %*% ALPHA
    DX2 <- K_nu %*% DX2
    DX <- K_nu * as.vector(DX1 - DX2)
  } else {
    stop("wrong model matrix names: ", m, "\n")
  }
  
  DX <- DX[, idx, drop=FALSE]
  DX
}



lav_object_inspect_coef <- function(object, type = "free",
                                    add.labels = FALSE, add.class = FALSE) {
  
  if(type == "user" || type == "all") {
    type <- "user"
    idx <- 1:length( object@ParTable$lhs )
  } else if(type == "free") {
    #idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
    idx <- which(object@ParTable$free > 0L)
  } else {
    stop("lavaan ERROR: argument `type' must be one of free or user")
  }
  EST <- lav_object_inspect_est(object)
  cof <- EST[idx]
  
  # labels?
  if(add.labels) {
    names(cof) <- lav_partable_labels(object@ParTable, type = type)
  }
  
  # class
  if(add.class) {
    class(cof) <- c("lavaan.vector", "numeric")
  }
  
  cof
}


lav_object_inspect_est <- function(object, unrotated = FALSE) {
  
  if(inherits(object, "lavaan")) {
    # from 0.5-19, they are in the partable
    if(!is.null(object@ParTable$est)) {
      if(unrotated) {
        OUT <- object@ParTable$est.unrotated
      } else {
        OUT <- object@ParTable$est
      }
    } else if(.hasSlot(object, "Fit")) {
      # in < 0.5-19, we should look in @Fit@est
      OUT <- object@Fit@est
    } else {
      PT <- parTable(object)
      OUT <- rep(as.numeric(NA), length(PT$lhs))
    }
  } else {
    # try generic coef()
    OUT <- coef(object, type = "user")
    if(is.matrix(OUT)) {
      # lavaanList?
      OUT <- rowMeans(OUT)
    }
  }
  
  OUT
}