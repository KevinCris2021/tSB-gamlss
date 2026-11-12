# =============================================================================
# tsb_distribution.R — t-Simplex-Binomial, three-parameter GAMLSS family.
# Consolidated candidate version, 2026-11-11.
#
# Model: Z | Y ~ Binomial(bd,Y), Y ~ t-Simplex(mu,sigma,nu).
# GAMLSS: mu (logit), sigma=phi (log), nu original (log).
#
# Load: Rcpp::sourceCpp("tsb_core.cpp"); source("tsb_distribution.R")
# This version must pass tests/test_tsb_distribution.R in a fresh RStudio session.
# =============================================================================

.tSB_need_cpp <- function() {
  if (!exists("tSB_moments_cpp", mode = "function", inherits = TRUE) ||
      !exists("tSB_information_cpp", mode = "function", inherits = TRUE)) {
    stop("Compile tsb_core.cpp first using Rcpp::sourceCpp().", call. = FALSE)
  }
}

.tSB_recycle <- function(args, n = NULL) {
  len <- lengths(args)
  if (any(len == 0L)) return(lapply(args, function(x) numeric(0)))
  if (is.null(n)) n <- max(len)
  if (any(!(len %in% c(1L, n))))
    stop("Each argument must have length 1 or n.", call. = FALSE)
  lapply(args, rep_len, length.out = n)
}

.tSB_check <- function(mu, sigma, nu, bd) {
  if (!is.numeric(mu) || any(!is.finite(mu) | mu <= 0 | mu >= 1))
    stop("mu must lie in (0,1).", call. = FALSE)
  if (!is.numeric(sigma) || any(!is.finite(sigma) | sigma <= 0))
    stop("sigma (phi) must be positive and finite.", call. = FALSE)
  if (!is.numeric(nu) || any(!is.finite(nu) | nu <= 0))
    stop("nu must be positive and finite.", call. = FALSE)
  if (!is.numeric(bd) || any(!is.finite(bd) | bd < 1 | bd != floor(bd)))
    stop("bd must contain positive integer trial counts.", call. = FALSE)
}

# Gauss-Legendre quadrature on (0,1), using Golub-Welsch.
tSB_GL_rule <- function(K = 64L) {
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) ||
      K != floor(K) || K < 8L)
    stop("K must be an integer >= 8.", call. = FALSE)
  K <- as.integer(K)
  j <- seq_len(K - 1L)
  e <- j/sqrt(4*j^2-1)
  J <- matrix(0, K, K)
  J[cbind(j,j+1L)] <- e
  J[cbind(j+1L,j)] <- e
  des <- eigen(J, symmetric = TRUE)
  u <- (des$values+1)/2
  w <- des$vectors[1L,]^2
  ind <- order(u)
  list(nodes = as.numeric(u[ind]), weights = as.numeric(w[ind]))
}

.tSB_rule <- function(regla) {
  if (!is.list(regla) || !all(c("nodes","weights") %in% names(regla)) ||
      !is.numeric(regla$nodes) || !is.numeric(regla$weights) ||
      length(regla$nodes) < 2L || length(regla$nodes)!=length(regla$weights) ||
      any(!is.finite(regla$nodes) | regla$nodes<=0 | regla$nodes>=1) ||
      any(!is.finite(regla$weights) | regla$weights<=0) ||
      abs(sum(regla$weights)-1) > 1e-8)
    stop("Invalid Gauss-Legendre rule.", call. = FALSE)
  regla
}

# -----------------------------------------------------------------------------
# d: discrete mass. Uses the unnormalized C++ likelihood kernel.
# -----------------------------------------------------------------------------
dtSB <- function(x, mu=0.5, sigma=1, nu=1, bd=10,
                 log=FALSE, regla=tSB_GL_rule(64L)) {
  .tSB_need_cpp(); regla <- .tSB_rule(regla)
  if (!is.logical(log) || length(log)!=1L || is.na(log))
    stop("log must be TRUE or FALSE.")
  a <- .tSB_recycle(list(x=x,mu=mu,sigma=sigma,nu=nu,bd=bd))
  if (!length(a$x)) return(numeric(0))
  .tSB_check(a$mu,a$sigma,a$nu,a$bd)
  if (!is.numeric(a$x)) stop("x must be numeric.")
  out <- rep(if (log) -Inf else 0, length(a$x))
  out[is.na(a$x)] <- NA_real_
  i <- which(is.finite(a$x) & a$x>=0 & a$x<=a$bd & a$x==floor(a$x))
  if (length(i)) {
    lp <- tSB_moments_cpp(
      z=as.numeric(a$x[i]),m=as.numeric(a$bd[i]),
      mu=as.numeric(a$mu[i]),phi=as.numeric(a$sigma[i]),
      nu=as.numeric(a$nu[i]),nodes=regla$nodes,weights=regla$weights
    )$log_prob
    out[i] <- if (log) lp else exp(lp)
  }
  out
}

# Normalized weights over the nodes; the density constant cancels.
# This normalizes the CDF at the quadrature level; it does not alter the C++ PMF.
.tSB_weights <- function(mu,sigma,nu,regla) {
  u <- regla$nodes
  d <- (u-mu)^2/(u*(1-u)*mu^2*(1-mu)^2)
  t <- nu+1; k <- nu+2.5
  q <- (d/(2*sigma))/t
  lw <- log(regla$weights)-1.5*(log(u)+log1p(-u))-k*log1p(q)
  if (anyNA(lw) || any(lw==Inf) || !any(is.finite(lw)))
    stop("tSB weights cannot be represented numerically.", call. = FALSE)
  w <- exp(lw-max(lw))
  w/sum(w)
}

# -----------------------------------------------------------------------------
# p: integral of the binomial CDF, evaluated directly in both tails.
# -----------------------------------------------------------------------------
ptSB <- function(q, mu=0.5,sigma=1,nu=1,bd=10,
                 lower.tail=TRUE,log.p=FALSE,regla=tSB_GL_rule(64L)) {
  regla <- .tSB_rule(regla)
  if (!is.logical(lower.tail) || length(lower.tail)!=1L || is.na(lower.tail) ||
      !is.logical(log.p) || length(log.p)!=1L || is.na(log.p))
    stop("lower.tail and log.p must be logical scalars.")
  a <- .tSB_recycle(list(q=q,mu=mu,sigma=sigma,nu=nu,bd=bd))
  if (!length(a$q)) return(numeric(0))
  .tSB_check(a$mu,a$sigma,a$nu,a$bd)
  if (!is.numeric(a$q)) stop("q must be numeric.")
  vapply(seq_along(a$q), function(i) {
    qi <- a$q[i]
    if (is.na(qi)) return(NA_real_)
    if (qi<0) p <- if (lower.tail) 0 else 1
    else if (qi>=a$bd[i]) p <- if (lower.tail) 1 else 0
    else {
      w <- .tSB_weights(a$mu[i],a$sigma[i],a$nu[i],regla)
      p <- sum(w * stats::pbinom(floor(qi),a$bd[i],regla$nodes,
                                 lower.tail=lower.tail))
      p <- max(0,min(1,p))
    }
    if (log.p) log(p) else p
  }, numeric(1))
}

# -----------------------------------------------------------------------------
# q: generalized inverse of the discrete CDF; consistent with p.
# -----------------------------------------------------------------------------
qtSB <- function(p, mu=0.5,sigma=1,nu=1,bd=10,
                 lower.tail=TRUE,log.p=FALSE,regla=tSB_GL_rule(64L)) {
  regla <- .tSB_rule(regla)
  if (!is.logical(lower.tail) || length(lower.tail)!=1L || is.na(lower.tail) ||
      !is.logical(log.p) || length(log.p)!=1L || is.na(log.p))
    stop("lower.tail and log.p must be logical scalars.")
  a <- .tSB_recycle(list(p=p,mu=mu,sigma=sigma,nu=nu,bd=bd))
  if (!length(a$p)) return(numeric(0))
  .tSB_check(a$mu,a$sigma,a$nu,a$bd)
  if (!is.numeric(a$p)) stop("p must be numeric.")
  p0 <- if (log.p) exp(a$p) else a$p
  if (any(!is.na(p0) & (p0<0 | p0>1))) stop("p is outside [0,1].")
  vapply(seq_along(p0), function(i) {
    pi <- p0[i]
    if (is.na(pi)) return(NA_real_)
    if (lower.tail && pi<=0) return(0)
    if (lower.tail && pi>=1) return(a$bd[i])
    if (!lower.tail && pi>=1) return(0)
    if (!lower.tail && pi<=0) return(a$bd[i])
    # Evaluate the full support once in C++ instead of repeating the CDF m times.
    soporte <- 0:as.integer(a$bd[i])
    pmf <- dtSB(soporte,mu=a$mu[i],sigma=a$sigma[i],nu=a$nu[i],
                bd=a$bd[i],regla=regla)
    pmf <- pmf/sum(pmf)
    F <- cumsum(pmf)
    F[length(F)] <- 1
    if (lower.tail) {
      as.numeric(soporte[which(F>=pi)[1L]])
    } else {
      S <- pmax(0,1-F)
      as.numeric(soporte[which(S<=pi)[1L]])
    }
  }, numeric(1))
}

# -----------------------------------------------------------------------------
# r: hierarchical simulation of t-Simplex followed by binomial sampling.
# The smaller root is evaluated stably, as in the original R Markdown.
# -----------------------------------------------------------------------------
rtSB <- function(n,mu=0.5,sigma=1,nu=1,bd=10) {
  if (length(n)>1L) n <- length(n)
  if (!is.numeric(n) || length(n)!=1L || !is.finite(n) ||
      n<0 || n!=floor(n) || n>.Machine$integer.max)
    stop("n must be a nonnegative integer.")
  n <- as.integer(n)
  if (!n) return(integer(0))
  a <- .tSB_recycle(list(mu=mu,sigma=sigma,nu=nu,bd=bd),n=n)
  .tSB_check(a$mu,a$sigma,a$nu,a$bd)
  mu <- a$mu; sigma <- a$sigma; nu <- a$nu
  epsilon <- mu/(1-mu)
  G <- stats::rgamma(n, shape=nu+2,rate=1)
  W <- sigma*(nu+1)/G
  lambda <- 1/(W*(1-mu)^2)
  V <- stats::rnorm(n)^2
  v <- epsilon*V/lambda
  x_small <- epsilon/(1+v/2+sqrt(v*(1+v/4)))
  prob_small <- ((1-mu)*epsilon+mu*x_small)/(epsilon+x_small)
  small <- stats::runif(n)<prob_small
  Y <- numeric(n)
  Y[small] <- x_small[small]/(1+x_small[small])
  j <- !small
  Y[j] <- 1/(1+x_small[j]/epsilon[j]^2)
  if (any(!is.finite(Y) | Y<=0 | Y>=1))
    stop("t-Simplex generator: probability is not numerically representable.")
  stats::rbinom(n,size=a$bd,prob=Y)
}

# -----------------------------------------------------------------------------
# tSB: GAMLSS constructor. All three parameters retain their original meaning.
#     Moments and information are cached in each family instance.
#     Closures are embedded in derivative functions to withstand
#     environment reassignment performed by gamlss().
# -----------------------------------------------------------------------------
tSB <- function(mu.link="logit",sigma.link="log",nu.link="log",
                regla=tSB_GL_rule(64L)) {
  .tSB_need_cpp()
  regla <- .tSB_rule(regla)
  mstats <- gamlss.dist::checklink("mu.link","tSB",substitute(mu.link),c("logit"))
  dstats <- gamlss.dist::checklink("sigma.link","tSB",substitute(sigma.link),c("log"))
  vstats <- gamlss.dist::checklink("nu.link","tSB",substitute(nu.link),c("log"))

  moments_cache <- new.env(parent=emptyenv())
  fisher_cache <- new.env(parent=emptyenv())
  reset <- function() {
    for (e in list(moments_cache,fisher_cache)) {
      e$key <- NULL; e$value <- NULL; e$calls <- 0L; e$hits <- 0L
    }
    invisible(NULL)
  }
  reset()

  rep_n <- function(x,n,nom) {
    if (!is.numeric(x) || !(length(x) %in% c(1L,n)))
      stop("Incompatible length: ",nom)
    rep_len(as.numeric(x),n)
  }
  get_moments <- function(y,bd,mu,sigma,nu) {
    n <- length(y)
    key <- list(z=as.numeric(y),m=rep_n(bd,n,"bd"),mu=rep_n(mu,n,"mu"),
                phi=rep_n(sigma,n,"sigma"),nu=rep_n(nu,n,"nu"))
    if (!is.null(moments_cache$key) && identical(key,moments_cache$key)) {
      moments_cache$hits <- moments_cache$hits+1L
      return(moments_cache$value)
    }
    ans <- tSB_moments_cpp(z=key$z,m=key$m,mu=key$mu,phi=key$phi,nu=key$nu,
                           nodes=regla$nodes,weights=regla$weights)
    moments_cache$key <- key; moments_cache$value <- ans
    moments_cache$calls <- moments_cache$calls+1L
    ans
  }
  get_fisher <- function(bd,mu,sigma,nu) {
    n <- max(lengths(list(bd,mu,sigma,nu)))
    key <- list(m=rep_n(bd,n,"bd"),mu=rep_n(mu,n,"mu"),
                phi=rep_n(sigma,n,"sigma"),nu=rep_n(nu,n,"nu"))
    if (!is.null(fisher_cache$key) && identical(key,fisher_cache$key)) {
      fisher_cache$hits <- fisher_cache$hits+1L
      return(fisher_cache$value)
    }
    ans <- tSB_information_cpp(m=key$m,mu=key$mu,phi=key$phi,nu=key$nu,
                               nodes=regla$nodes,weights=regla$weights)$information
    fisher_cache$key <- key; fisher_cache$value <- ans
    fisher_cache$calls <- fisher_cache$calls+1L
    ans
  }
  create_score <- function(column) {
    f <- function(y,mu,sigma,nu,bd,...) NULL
    body(f) <- substitute(as.numeric(GET(y=y,bd=bd,mu=mu,sigma=sigma,nu=nu)$score[,COL]),
                          list(GET=get_moments,COL=column))
    f
  }
  create_curvature <- function(column) {
    f <- function(y,mu,sigma,nu,bd,...) NULL
    body(f) <- substitute(-as.numeric(GET(bd=bd,mu=mu,sigma=sigma,nu=nu)[,COL]),
                          list(GET=get_fisher,COL=column))
    f
  }
  fdev <- function(y,mu,sigma,nu,bd,...) NULL
  body(fdev) <- substitute(-2*GET(y=y,bd=bd,mu=mu,sigma=sigma,nu=nu)$log_prob,
                           list(GET=get_moments))

  # Discrete mean and variance under the mixture probability law.
  meanfun <- function(bd,mu,sigma,nu) bd*mu
  varfun <- function(bd,mu,sigma,nu) {
    a <- .tSB_recycle(list(bd=bd,mu=mu,sigma=sigma,nu=nu))
    .tSB_check(a$mu,a$sigma,a$nu,a$bd)
    v <- vapply(seq_along(a$mu),function(i) {
      w <- .tSB_weights(a$mu[i],a$sigma[i],a$nu[i],regla)
      sum(w*(regla$nodes-a$mu[i])^2)
    },numeric(1))
    a$bd*a$mu*(1-a$mu)+a$bd*(a$bd-1)*v
  }

  fam <- list(
    family=c("tSB","t-Simplex-Binomial","Binomial"),
    parameters=list(mu=TRUE,sigma=TRUE,nu=TRUE),nopar=3,type="Discrete",
    mu.link=as.character(substitute(mu.link)),
    sigma.link=as.character(substitute(sigma.link)),
    nu.link=as.character(substitute(nu.link)),
    mu.linkfun=mstats$linkfun,sigma.linkfun=dstats$linkfun,nu.linkfun=vstats$linkfun,
    mu.linkinv=mstats$linkinv,sigma.linkinv=dstats$linkinv,nu.linkinv=vstats$linkinv,
    mu.dr=mstats$mu.eta,sigma.dr=dstats$mu.eta,nu.dr=vstats$mu.eta,
    dldm=create_score("mu"),dldd=create_score("phi"),dldv=create_score("nu"),
    d2ldm2=create_curvature("mumu"),
    d2ldd2=create_curvature("phiphi"),
    d2ldv2=create_curvature("nunu"),
    d2ldmdd=create_curvature("muphi"),
    d2ldmdv=create_curvature("munu"),
    d2ldddv=create_curvature("phinu"),
    G.dev.incr=fdev,
    rqres=expression(rqres(pfun="ptSB",type="Discrete",ymin=0,
                          y=y,mu=mu,sigma=sigma,nu=nu,bd=bd)),
    mu.initial=expression({mu <- (y+0.5)/(bd+1)}),
    sigma.initial=expression({sigma <- rep(0.8,length(y))}),
    nu.initial=expression({nu <- rep(2,length(y))}),
    mu.valid=function(mu) all(is.finite(mu) & mu>0 & mu<1),
    sigma.valid=function(sigma) all(is.finite(sigma) & sigma>0),
    nu.valid=function(nu) all(is.finite(nu) & nu>0),
    y.valid=function(y) all(is.finite(y) & y>=0 & y==floor(y)),
    mean=meanfun,variance=varfun
  )
  fam <- structure(fam,class=c("gamlss.family","family"))
  attr(fam,"cache_control") <- list(
    reiniciar=reset,
    estadisticas=function() c(momentos_calculos=moments_cache$calls,
                              momentos_aciertos=moments_cache$hits,
                              fisher_calculos=fisher_cache$calls,
                              fisher_aciertos=fisher_cache$hits)
  )
  fam
}
