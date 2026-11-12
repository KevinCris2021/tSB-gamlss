# =============================================================================
# Self-contained regression-gradient test for the extended evaluator.
# Requires run_all_tests.R to compile C++ and load both R source files first.
# Does not use objects from the historical notebook or fit models.
# =============================================================================

set.seed(5156)
regla_mixta <- tSB_GL_rule(64L)
datos_mixtos <- data.frame(
  x = seq(-1, 1, length.out=9L),
  m = rep(35L,9L),
  z = c(6,9,11,13,16,18,21,24,27)
)
X_mixta <- stats::model.matrix(~ x, datos_mixtos)
theta_mixta <- c(beta0=-0.25,beta1=0.50,alpha0=log(0.9),
                 gamma0=0,gamma1=40)

# gamma parametrizes eta_tau=logit(tau), NOT log(nu).
evaluar_mixta <- function(theta) {
  eta_mu <- drop(X_mixta %*% theta[1:2])
  eta_tau <- drop(X_mixta %*% theta[4:5])
  mu <- stats::plogis(eta_mu)
  phi <- rep(exp(theta[3]),nrow(datos_mixtos))
  mom <- tSB_tau_moments(
    z=datos_mixtos$z,m=datos_mixtos$m,mu=mu,phi=phi,
    eta_tau=eta_tau,regla=regla_mixta,umbral=35
  )
  grad <- c(drop(crossprod(X_mixta,mom$score_eta_mu)),
            sum(mom$score_eta_phi),
            drop(crossprod(X_mixta,mom$score_eta_tau)))
  names(grad) <- names(theta_mixta)
  list(ell=sum(mom$log_prob),grad=grad,mom=mom,
       mu=mu,phi=phi,eta_tau=eta_tau)
}

resultado_mixto <- evaluar_mixta(theta_mixta)
conteos_regiones <- table(factor(resultado_mixto$mom$region,
                                levels=c('SB','tSB','tSB0')))
stopifnot(identical(as.integer(conteos_regiones),c(1L,7L,1L)))

# The numerical gradient is computed from the SAME log-likelihood.
# The steps are small enough NOT to cross eta_tau=+/-35.
h <- 1e-5*pmax(1,abs(theta_mixta))
grad_numerico <- vapply(seq_along(theta_mixta), function(j) {
  a <- theta_mixta; b <- theta_mixta
  a[j] <- a[j]+h[j]; b[j] <- b[j]-h[j]
  (evaluar_mixta(a)$ell-evaluar_mixta(b)$ell)/(2*h[j])
}, numeric(1))
error_gradiente_mixto <- max(abs(resultado_mixto$grad-grad_numerico))
stopifnot(is.finite(error_gradiente_mixto),error_gradiente_mixto<1e-5)

# Fisher information for eta_tau must be exactly zero at both limits.
F_mixta <- tSB_tau_information(
  m=datos_mixtos$m,mu=resultado_mixto$mu,
  phi=resultado_mixto$phi,eta_tau=resultado_mixto$eta_tau,
  regla=regla_mixta,umbral=35
)
I_tau <- as.numeric(F_mixta$information['eta_tau','eta_tau',])
stopifnot(all(I_tau[resultado_mixto$mom$region!='tSB']==0))

# Expected information for gamma=(gamma0,gamma1)
I_gamma <- crossprod(X_mixta,X_mixta*I_tau)
stopifnot(max(abs(I_gamma-t(I_gamma)))<1e-12,
          all(is.finite(I_gamma)),
          abs(I_gamma[2,2]-sum(datos_mixtos$x^2*I_tau))<1e-12)

# If every row belongs to SB, gamma coefficients do not enter the objective.
theta_SB <- theta_mixta
theta_SB[c('gamma0','gamma1')] <- c(-40,0)
resultado_SB <- evaluar_mixta(theta_SB)
stopifnot(all(resultado_SB$mom$region=='SB'),
          all(resultado_SB$grad[4:5]==0))

# Sign equivalence: eta_tau=-eta_nu at an interior point.
theta_nu <- c(beta0=-0.3,beta1=0.6,alpha0=log(0.9),
              gamma0=log(1.4),gamma1=0.9)
theta_tau <- theta_nu
theta_tau[4:5] <- -theta_nu[4:5]
R_tau <- evaluar_mixta(theta_tau)
nu <- exp(drop(X_mixta %*% theta_nu[4:5]))
R_nu <- tSB_moments_cpp(
  z=datos_mixtos$z,m=datos_mixtos$m,
  mu=R_tau$mu,phi=R_tau$phi,nu=nu,
  nodes=regla_mixta$nodes,weights=regla_mixta$weights
)
grad_nu <- c(drop(crossprod(X_mixta,R_nu$score[,'mu']*R_tau$mu*(1-R_tau$mu))),
             sum(R_nu$score[,'phi']*R_tau$phi),
             drop(crossprod(X_mixta,R_nu$score[,'nu']*nu)))
grad_nu[4:5] <- -grad_nu[4:5]
error_signo <- max(abs(R_tau$grad-grad_nu))
error_ell <- abs(R_tau$ell-sum(R_nu$log_prob))
stopifnot(error_signo<1e-9,error_ell<1e-9)

cat('\nMIXED REGRESSION: regions=',paste(names(conteos_regiones),
    as.integer(conteos_regiones),collapse='; '),'\n')
cat('Central-difference gradient error:',
    format(error_gradiente_mixto,digits=12),'\n')
cat('Log-likelihood error under nu/tau transformation:',format(error_ell,digits=12),'\n')
cat('Gradient error under nu/tau transformation:',
    format(error_signo,digits=12),'\n')
cat('Expected information for (gamma0,gamma1):\n')
print(I_gamma,digits=12)
cat('MIXED-REGRESSION TEST PASSED\n')
