# =============================================================================
# MASTER TEST, reproducible from a clean R session.
# Run from the project root directory:
#   source("run_all_tests.R")                 # RStudio / R
#   Rscript run_all_tests.R                # command line
#   Rscript run_all_tests.R --fit       # include population GAMLSS fit
# Does not depend on objects saved by the research R Markdown.
# =============================================================================

necesarios <- c('Rcpp','gamlss','gamlss.dist')
ausentes <- necesarios[!vapply(necesarios,requireNamespace,
                              logical(1),quietly=TRUE)]
if (length(ausentes)) stop('Install first: ',paste(ausentes,collapse=', '))
stopifnot(all(file.exists(c('tsb_core.cpp','tsb_distribution.R','tsb_tau_limits.R',
                           'tests/test_mixed_regression.R'))))

cat('======= STANDALONE LOAD =======\n')
Rcpp::sourceCpp('tsb_core.cpp')
source('tsb_distribution.R')
source('tsb_tau_limits.R')
regla_verificacion <- tSB_GL_rule(64L)

# A. Discrete distribution and quadrature
z <- 0:30
pmf <- dtSB(z,mu=.5,sigma=.8,nu=2,bd=30,regla=regla_verificacion)
CDF_acumulada <- ptSB(z,mu=.5,sigma=.8,nu=2,bd=30,regla=regla_verificacion)
CDF_complementaria <- ptSB(z,mu=.5,sigma=.8,nu=2,bd=30,lower.tail=FALSE,
           regla=regla_verificacion)
error_masa <- abs(sum(pmf)-1)
error_colas <- max(abs(CDF_acumulada+CDF_complementaria-1))
error_CDF <- max(abs(CDF_acumulada-cumsum(pmf)))
stopifnot(error_masa<1e-6,error_colas<1e-10,error_CDF<1e-6)
cat('Interior distribution (mass, tails, CDF):',
    format(c(error_masa,error_colas,error_CDF),digits=10),'\n')

# B. Evaluation at both limits and in the interior.
casos <- data.frame(
  z=c(5,8,10,15,20,23,26),m=rep(30,7),
  mu=c(.35,.4,.45,.5,.55,.6,.65),
  phi=c(.8,.9,1,1.1,1.2,1.3,1.4),
  eta_tau=c(-Inf,-40,-35,0,35,40,Inf)
)
M <- tSB_tau_moments(
  z=casos$z,m=casos$m,mu=casos$mu,phi=casos$phi,
  eta_tau=casos$eta_tau,regla=regla_verificacion
)
stopifnot(identical(M$region,c('SB','SB','SB','tSB','tSB0','tSB0','tSB0')),
          all(is.finite(as.matrix(M[,1:4]))),
          all(M$score_eta_tau[M$region!='tSB']==0))
Fisher_limites <- tSB_tau_information(
  m=casos$m,mu=casos$mu,phi=casos$phi,
  eta_tau=casos$eta_tau,regla=regla_verificacion
)
error_masa_ext <- max(abs(Fisher_limites$mass-1))
error_score_ext <- max(abs(Fisher_limites$expected_score))
stopifnot(error_masa_ext<1e-6,error_score_ext<1e-5,
          all(Fisher_limites$information['eta_tau','eta_tau',M$region!='tSB']==0))
cat('Boundaries (mass, expected score):',
    format(c(error_masa_ext,error_score_ext),digits=10),'\n')

# C. nu/tau equivalence at an interior point and exact limits.
i <- which(M$region=='tSB')
nu_i <- exp(-casos$eta_tau[i])
cpp <- tSB_moments_cpp(
  z=casos$z[i],m=casos$m[i],mu=casos$mu[i],
  phi=casos$phi[i],nu=nu_i,
  nodes=regla_verificacion$nodes,weights=regla_verificacion$weights
)
error_interior <- max(abs(cpp$log_prob-M$log_prob[i]))
p_extremos <- dtSB_tau(
  x=c(8,15,20),mu=c(.3,.5,.7),sigma=c(.8,.9,1.2),
  tau=c(0,.5,1),bd=30,regla=regla_verificacion
)
stopifnot(error_interior<1e-10,all(is.finite(p_extremos) &
                                      p_extremos>0 & p_extremos<1))
cat('Interior equivalence:',format(error_interior,digits=12),'\n')

# D. Unmodified GAMLSS family.
familia <- tSB(regla=regla_verificacion)
stopifnot(inherits(familia,'gamlss.family'),
          identical(familia$nu.link,'log'),
          identical(familia$nopar,3))
cat('GAMLSS family: nu.link=',familia$nu.link,', nopar=',
    familia$nopar,'\n',sep='')

# E. Regression across three regions and gradient equivalence.
source('tests/test_mixed_regression.R')

# F. Population-weighted GAMLSS fit (optional, slower).
# '--fit' is the public CLI option; '--ajuste' remains a legacy alias.
# CLI flags are interpreted only when run from Rscript.
args <- commandArgs(trailingOnly=TRUE)
if (any(c('--fit', '--ajuste') %in% args)) {
  cat('\n======= POPULATION GAMLSS (OPTIONAL) =======\n')
  source('examples/fit_population_gamlss.R')
  stopifnot(exists('ajuste_5146'),!is.null(ajuste_5146),
            isTRUE(ajuste_5146$converged),
            is.finite(error_logLik),error_logLik<1e-5,
            max(abs(estimaciones_5146-parametros_5146))<0.1)
}

cat('\nALL SELECTED TESTS COMPLETED SUCCESSFULLY\n')
