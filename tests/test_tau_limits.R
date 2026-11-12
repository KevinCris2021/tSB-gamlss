# Run from the folder containing tsb_core.cpp and both R source files.
Rcpp::sourceCpp('tsb_core.cpp')
source('tsb_distribution.R')
source('tsb_tau_limits.R')
GL64 <- tSB_GL_rule(64L)

casos <- data.frame(z=c(5,8,10,15,20,23,26), m=rep(30,7),
                    mu=c(.35,.4,.45,.5,.55,.6,.65),
                    phi=c(.8,.9,1,1.1,1.2,1.3,1.4),
                    eta_tau=c(-Inf,-40,-35,0,35,40,Inf))
R <- tSB_tau_moments(casos$z,casos$m,casos$mu,casos$phi,casos$eta_tau,
                     regla=GL64)
stopifnot(identical(R$region,c('SB','SB','SB','tSB','tSB0','tSB0','tSB0')),
          all(is.finite(as.matrix(R[,1:4]))),
          all(R$score_eta_tau[R$region != 'tSB']==0))

# Equivalence with the interior kernel
medio <- which(R$region=='tSB')
nu <- exp(-casos$eta_tau[medio])
ref <- tSB_moments_cpp(z=casos$z[medio],m=casos$m[medio],
                       mu=casos$mu[medio],phi=casos$phi[medio],nu=nu,
                       nodes=GL64$nodes,weights=GL64$weights)
error_interior <- abs(R$log_prob[medio]-ref$log_prob)
stopifnot(error_interior < 1e-10)

I <- tSB_tau_information(casos$m,casos$mu,casos$phi,casos$eta_tau,
                         regla=GL64)
informacion_forma <- I$information['eta_tau','eta_tau',]
stopifnot(all(informacion_forma[R$region!='tSB']==0),
          all(is.finite(I$information)),
          max(abs(I$mass-1)) < 1e-6,
          max(abs(I$expected_score)) < 1e-5)

# Exact-tau interfaces cover both limits.
p <- dtSB_tau(x=c(8,15,20),mu=c(.3,.5,.7),sigma=c(.8,.9,1.2),
              tau=c(0,.5,1),bd=30,regla=GL64)
stopifnot(all(is.finite(p) & p>0 & p<1))

# The original GAMLSS constructor remains unchanged.
familia <- tSB(regla=GL64)
stopifnot(inherits(familia,'gamlss.family'),
          identical(familia$nu.link,'log'),
          identical(familia$nopar,3))

cat('TAU BOUNDARIES VALIDATED; GAMLSS FAMILY UNCHANGED\n')
print(data.frame(eta_tau=casos$eta_tau,region=R$region,
                 Fisher_shape=informacion_forma),row.names=FALSE)
