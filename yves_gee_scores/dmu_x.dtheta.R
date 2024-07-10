dmu_x.dtheta <- function(x, lavmodel = NULL, x.i = NULL) {

  type <- "free"
  nvar <- lavmodel@nvar
  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  GLIST <- lavmodel@GLIST

  # number of rows in DELTA.group
  pstar <- integer(nblocks)
  for (g in 1:nblocks) {
    # only means
    pstar[g] <- nvar[g]
  }

  # number of columns in DELTA + m.el.idx/x.el.idx
  if (.hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
    NCOL <- lavmodel@nx.unco
  } else {
    NCOL <- lavmodel@nx.free
  }
  m.el.idx <- x.el.idx <- vector("list", length = length(GLIST))
  for (mm in 1:length(GLIST)) { 
    m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
    if (.hasSlot(lavmodel, "ceq.simple.only") &&
      lavmodel@ceq.simple.only) {
      x.el.idx[[mm]] <- lavmodel@x.unco.idx[[mm]]
    } else {
      x.el.idx[[mm]] <- lavmodel@x.free.idx[[mm]]
    }
    # handle symmetric matrices
    if (lavmodel@isSymmetric[mm]) {
      # since we use 'x.free.idx', only symmetric elements
      # are duplicated (not the equal ones, only in x.free.free)
      dix <- duplicated(x.el.idx[[mm]])
      if (any(dix)) {
        m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
        x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
      }
    }
  }

  # compute Delta
  Delta <- vector("list", length = nblocks)
  for (g in 1:nblocks) {
    Delta.group <- matrix(0, nrow = pstar[g], ncol = NCOL)

    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]

    for (mm in mm.in.group) {
      mname <- names(lavmodel@GLIST)[mm]

      # skip empty ones
      if (!length(m.el.idx[[mm]])) next

      DELTA.mu <- derivative.mu_xi.LISREL(
              m = mname,
              idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group],
              x.i = x.i
            )
      Delta.group[, x.el.idx[[mm]]] <- DELTA.mu
    } # mm

    if (type == "free" &&
      .hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      Delta.group <- Delta.group %*% lavmodel@ceq.simple.K
    }

    Delta[[g]] <- Delta.group
  } # g

  Delta
}

# dMu.xi/dx -- per model matrix
derivative.mu_xi.LISREL <- function(m = "alpha",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 MLIST = NULL,
                                 x.i = NULL) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  GAMMA <- MLIST$gamma
  stopifnot(!is.null(GAMMA))
  x.i <- as.matrix(x.i)
  GX <- GAMMA %*% x.i

  # shortcut for empty matrices
  if (m == "psi" || m == "theta" || m == "tau" ||
    m == "delta" || m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = nvar, ncol = length(idx)))
  }

  # missing alpha
  if (is.null(MLIST$alpha)) {
    ALPHA <- matrix(0, nfac, 1L)
  } else {
    ALPHA <- MLIST$alpha
  }

  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- lavaan:::.internal_get_IB.inv(MLIST = MLIST)
  }
  IB.inv.Gamma <- IB.inv %*% GAMMA

  if (m == "nu") {
    DX <- diag(nvar)
  } else if (m == "lambda") {
    DX <- ( t(IB.inv %*% ALPHA) %x% diag(nvar) +
            t(IB.inv.Gamma %*% x.i) %x% diag(nvar) )
  } else if (m == "beta") {
    # double check!!
    DX <- ( t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv) +
            t(IB.inv.Gamma %*% x.i) %x% (LAMBDA %*% IB.inv.Gamma) )
  } else if (m == "alpha") {
    DX <- LAMBDA %*% IB.inv
  } else if (m == "gamma") {
    DX <- t(x.i) %x% (LAMBDA %*% IB.inv)
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }
  DX <- DX[, idx, drop = FALSE]
  DX
}

