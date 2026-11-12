# =============================================================================
# tsb_tau_limits.R | Boundary extension for t-Simplex-Binomial
# INDEPENDENT layer separate from the GAMLSS tSB() constructor.
#
# Boundary parametrization:
#   eta_tau = logit(tau) = -log(nu)
#   eta_tau <= -threshold  -> SB    (tau = 0)
#   -threshold < eta_tau < threshold -> tSB (0 < tau < 1)
#   eta_tau >= threshold   -> tSB0  (tau = 1)
#
# Derivatives: scores with respect to eta_mu=logit(mu),
# eta_phi=log(phi), eta_tau. Switched regions have S_eta_tau=0.
#
# REQUIRES: source("tsb_distribution.R") and Rcpp::sourceCpp("tsb_core.cpp") beforehand.
# Does NOT replace or modify tSB(), dtSB(), ptSB(), qtSB(), or rtSB().
# The piecewise function must NOT be inserted directly into RS-GAMLSS:
# expected information for the shape predictor vanishes at the limits.
# =============================================================================

.tSB_tau_prepare <- function(args) {
  longitudes <- lengths(args)
  if (any(longitudes == 0L)) return(NULL)
  n <- max(longitudes)
  if (any(!(longitudes %in% c(1L, n)))) {
    stop("Arguments must have length 1 or n.", call. = FALSE)
  }
  lapply(args, function(x) rep_len(as.numeric(x), n))
}

.tSB_tau_check <- function(a, valor, umbral) {
  if (!is.numeric(a[[valor]]) || (valor != "x" && anyNA(a[[valor]]))) {
    stop("Argument ", valor, " must be numeric and contain no NA values.", call. = FALSE)
  }
  if (any(!is.finite(a$mu) | a$mu <= 0 | a$mu >= 1) ||
      any(!is.finite(a$phi) | a$phi <= 0) ||
      any(is.na(a$eta_tau)) ||
      any(!is.finite(a$m) | a$m < 1 | a$m != floor(a$m)) ||
      length(umbral) != 1L || !is.finite(umbral) || umbral <= 0) {
    stop("Parameters outside their domain.", call. = FALSE)
  }
}

# Moments at the SB or tSB0 limit, analytically integrated with respect
# to the mixing density and numerically using Gauss-Legendre quadrature.
.tSB_tau_endpoint_one <- function(z, m, mu, phi, lado, regla) {
  u <- regla$nodes
  log_u1 <- log(u) + log1p(-u)
  r <- u - mu
  c <- 1 / (u * (1-u) * (mu*(1-mu))^2)
  d <- c * r^2
  A <- -2/mu + 2/(1-mu)
  d_mu <- c * (-2*r + A*r^2)
  if (lado == "SB") {
    log_f <- -0.5 * (log(2*pi*phi) + 3*log_u1 + d/phi)
    S_mu <- -d_mu/(2*phi)
    S_phi <- (-0.5 + d/(2*phi))/phi
  } else {
    q <- d/(2*phi)
    log_f <- lgamma(2.5)-lgamma(2) -
      0.5*(log(2*pi*phi) + 3*log_u1) - 2.5*log1p(q)
    S_mu <- -2.5*d_mu/(2*phi*(1+q))
    S_phi <- (-0.5 + 2.5*q/(1+q))/phi
  }
  lw <- log(regla$weights) + stats::dbinom(z, m, u, log=TRUE) + log_f
  max_lw <- max(lw)
  if (!is.finite(max_lw)) stop("Non-finite limiting integrand.", call. = FALSE)
  w <- exp(lw-max_lw)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) stop("Invalid limiting integral.", call. = FALSE)
  w <- w/sw
  c(log_prob=max_lw+log(sw),
    score_eta_mu=sum(w*S_mu)*mu*(1-mu),
    score_eta_phi=sum(w*S_phi)*phi,
    score_eta_tau=0)
}

# Main evaluator on the eta_tau scale: accepts exact -Inf and +Inf.
tSB_tau_moments <- function(z, m, mu, phi, eta_tau,
                            regla=tSB_GL_rule(64L), umbral=35) {
  .tSB_need_cpp()
  regla <- .tSB_rule(regla)
  a <- .tSB_tau_prepare(list(z=z, m=m, mu=mu, phi=phi,
                             eta_tau=eta_tau))
  if (is.null(a)) {
    return(data.frame(log_prob=numeric(0), score_eta_mu=numeric(0),
                      score_eta_phi=numeric(0), score_eta_tau=numeric(0),
                      region=character(0)))
  }
  .tSB_tau_check(a, "z", umbral)
  if (any(!is.finite(a$z) | a$z < 0 | a$z > a$m | a$z != floor(a$z))) {
    stop("z must be an integer count within [0,m].", call. = FALSE)
  }
  n <- length(a$z)
  region <- ifelse(a$eta_tau <= -umbral, "SB",
                   ifelse(a$eta_tau >= umbral, "tSB0", "tSB"))
  ans <- data.frame(log_prob=numeric(n),score_eta_mu=numeric(n),
                    score_eta_phi=numeric(n),score_eta_tau=numeric(n),
                    region=region)
  for (lado in c("SB", "tSB0")) {
    ix <- which(region == lado)
    if (length(ix)) {
      val <- vapply(ix, function(i) .tSB_tau_endpoint_one(
        a$z[i], a$m[i], a$mu[i], a$phi[i], lado, regla
      ), numeric(4L))
      ans[ix,1:4] <- t(val)
    }
  }
  ix <- which(region == "tSB")
  if (length(ix)) {
    nu <- exp(-a$eta_tau[ix])
    momentos <- tSB_moments_cpp(z=a$z[ix],m=a$m[ix],mu=a$mu[ix],
                                phi=a$phi[ix],nu=nu,
                                nodes=regla$nodes,weights=regla$weights)
    ans$log_prob[ix] <- momentos$log_prob
    ans$score_eta_mu[ix] <- momentos$score[,"mu"]*a$mu[ix]*(1-a$mu[ix])
    ans$score_eta_phi[ix] <- momentos$score[,"phi"]*a$phi[ix]
    ans$score_eta_tau[ix] <- -nu*momentos$score[,"nu"]
  }
  if (any(!is.finite(as.matrix(ans[,1:4])))) {
    stop("Non-finite log-PMF or scores.",call.=FALSE)
  }
  ans
}

# Extended mass on the eta_tau scale, zero outside the support.
dtSB_eta_tau <- function(x, mu, sigma, eta_tau, bd,
                         log=FALSE, regla=tSB_GL_rule(64L), umbral=35) {
  if (!is.logical(log) || length(log)!=1L || is.na(log))
    stop("log must be TRUE or FALSE.",call.=FALSE)
  a <- .tSB_tau_prepare(list(x=x,mu=mu,phi=sigma,eta_tau=eta_tau,m=bd))
  if (is.null(a)) return(numeric(0))
  .tSB_tau_check(a,"x",umbral)
  out <- rep(if (log) -Inf else 0, length(a$x))
  out[is.na(a$x)] <- NA_real_
  ix <- which(is.finite(a$x) & a$x >= 0 & a$x <= a$m & a$x==floor(a$x))
  if (length(ix)) {
    lp <- tSB_tau_moments(z=a$x[ix],m=a$m[ix],mu=a$mu[ix],
                          phi=a$phi[ix],eta_tau=a$eta_tau[ix],
                          regla=regla,umbral=umbral)$log_prob
    out[ix] <- if (log) lp else exp(lp)
  }
  out
}

# Same model using tau; when tau is near 1 and rounded to 1,
# prefer dtSB_eta_tau() with the original unrounded predictor.
dtSB_tau <- function(x, mu, sigma, tau, bd,
                      log=FALSE, regla=tSB_GL_rule(64L), umbral=35) {
  if (!is.numeric(tau) || anyNA(tau) || any(tau < 0 | tau > 1))
    stop("tau must lie in [0,1].",call.=FALSE)
  dtSB_eta_tau(x=x,mu=mu,sigma=sigma,eta_tau=stats::qlogis(tau),
               bd=bd,log=log,regla=regla,umbral=umbral)
}

# Expected information for the predictors as the outer product of
# the score of the SAME piecewise objective function.
tSB_tau_information <- function(m, mu, phi, eta_tau,
                                regla=tSB_GL_rule(64L), umbral=35,
                                tolerancia_masa=1e-6) {
  .tSB_need_cpp()
  regla <- .tSB_rule(regla)
  a <- .tSB_tau_prepare(list(m=m,mu=mu,phi=phi,eta_tau=eta_tau))
  if (is.null(a)) stop("Empty arguments.",call.=FALSE)
  .tSB_tau_check(a,"m",umbral)
  n <- length(a$m)
  nombres <- c("eta_mu","eta_phi","eta_tau")
  I <- array(0,dim=c(3L,3L,n),dimnames=list(nombres,nombres,NULL))
  E <- matrix(0,nrow=n,ncol=3L,dimnames=list(NULL,nombres))
  masa <- numeric(n)
  for (i in seq_len(n)) {
    if (a$m[i] > .Machine$integer.max-1L)
      stop("m is too large to construct its support.",call.=FALSE)
    soporte <- 0:as.integer(a$m[i])
    moment <- tSB_tau_moments(z=soporte,m=a$m[i],mu=a$mu[i],
                              phi=a$phi[i],eta_tau=a$eta_tau[i],
                              regla=regla,umbral=umbral)
    p <- exp(moment$log_prob)
    S <- as.matrix(moment[,c("score_eta_mu","score_eta_phi","score_eta_tau")])
    colnames(S) <- nombres
    I[,,i] <- crossprod(S, S * p)
    E[i,] <- colSums(S * p)
    masa[i] <- sum(p)
  }
  if (any(!is.finite(I)) || any(!is.finite(E)) ||
      any(!is.finite(masa)) || max(abs(masa-1)) > tolerancia_masa)
    stop("Invalid information or normalization: check the GL rule.",call.=FALSE)
  list(information=I,expected_score=E,mass=masa)
}
