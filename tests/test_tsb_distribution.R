# =============================================================================
# Validation of the standalone tsb_core.cpp + tsb_distribution.R port
# Run in a FRESH RStudio session, from this directory.
# No objects from the research R Markdown are required.
# =============================================================================

stopifnot(requireNamespace("Rcpp",quietly=TRUE),
          requireNamespace("gamlss",quietly=TRUE),
          requireNamespace("gamlss.dist",quietly=TRUE))
Rcpp::sourceCpp("tsb_core.cpp")
source("tsb_distribution.R")
GL64 <- tSB_GL_rule(64L)
GL96 <- tSB_GL_rule(96L)

# A. PMF and CDF, both tails; GL64/GL96 comparison.
m <- 30L; z <- 0:m
p64 <- dtSB(z,mu=.5,sigma=.8,nu=2,bd=m,regla=GL64)
p96 <- dtSB(z,mu=.5,sigma=.8,nu=2,bd=m,regla=GL96)
F <- ptSB(z,mu=.5,sigma=.8,nu=2,bd=m,regla=GL64)
S <- ptSB(z,mu=.5,sigma=.8,nu=2,bd=m,lower.tail=FALSE,regla=GL64)
errA <- c(masa=abs(sum(p64)-1), colas=max(abs(F+S-1)),
          CDF=max(abs(F-cumsum(p64))), GL=max(abs(p64-p96)))
cat("A. d/p functions and GL64/96:\n"); print(stats::setNames(errA, c("mass", "tails", "CDF", "GL")), digits=12)
stopifnot(errA["masa"]<1e-6,errA["colas"]<1e-10,
          errA["CDF"]<1e-6,errA["GL"]<1e-6)
q <- qtSB(c(.05,.50,.95),mu=.5,sigma=.8,nu=2,bd=m,regla=GL64)
cat("   Quantiles: ",paste(q,collapse=", "),"\n",sep="")

# B. Scores, Hessian, and Fisher information; regular and extreme nu.
casos <- data.frame(z=c(4,12,20,31,48),m=c(20,30,40,60,80),
                    mu=c(.2,.35,.5,.65,.8),phi=c(.5,.8,1,1.5,2),
                    nu=c(.1,.5,1,5,20))
res <- tSB_moments_cpp(casos$z,casos$m,casos$mu,casos$phi,casos$nu,
                       GL64$nodes,GL64$weights)
# Use GL96 for the expected information check in the challenging nu = 0.1 case.
# GL64 is retained for the observed moments above and the other tests below.
inf <- tSB_information_cpp(casos$m,casos$mu,casos$phi,casos$nu,
                           GL96$nodes,GL96$weights)
stopifnot(all(is.finite(res$log_prob)),all(is.finite(res$score)),
          all(is.finite(res$hessian)),all(is.finite(inf$information)),
          max(abs(inf$mass-1))<1e-6,
          max(abs(inf$expected_score))<1e-6)
cat("B. Moments and information: OK\n")

# C. Asymptotic behavior: test finiteness and normalization.
lim <- data.frame(m=c(60,70),mu=c(.35,.65),phi=c(.7,.7),
                  nu=c(1e-20,1e15))
lim_i <- tSB_information_cpp(lim$m,lim$mu,lim$phi,lim$nu,
                             GL64$nodes,GL64$weights)
stopifnot(all(is.finite(lim_i$information)),
          max(abs(lim_i$mass-1))<1e-6)
cat("C. Extreme parameters: OK\n")

# D. Compare hierarchical simulation against the computed marginal mean.
set.seed(5145)
N <- 20000L
s <- rtSB(N,mu=.3,sigma=.8,nu=2,bd=30)
p <- dtSB(0:30,mu=.3,sigma=.8,nu=2,bd=30,regla=GL96)
media_ref <- sum((0:30)*p)
var_ref <- sum(((0:30)-media_ref)^2*p)
err_media <- abs(mean(s)-media_ref)
cat("D. Simulated-mean error:",err_media,"\n")
stopifnot(err_media<5*sqrt(var_ref/N))

# E. GAMLSS family with constant nu (integration test).
set.seed(5146)
n <- 200L
x <- runif(n); m <- sample(30:60,n,replace=TRUE)
mu <- plogis(-.8+1.6*x)
z <- rtSB(n,mu=mu,sigma=.7,nu=7,bd=m)
datos <- data.frame(z=z,m=m,x=x)
ini <- suppressWarnings(stats::glm(cbind(z,m-z)~x,
                    family=stats::binomial(),data=datos))
fam <- tSB(regla=GL64)
fit <- gamlss::gamlss(cbind(z,m-z)~x,sigma.formula=~1,nu.formula=~1,
                      family=fam,data=datos,
                      mu.start=as.numeric(stats::fitted(ini)),
                      sigma.start=.7,nu.start=7,
                      control=gamlss::gamlss.control(n.cyc=50,c.crit=1e-6,
                                                     trace=FALSE),trace=FALSE)
mu1 <- as.numeric(stats::fitted(fit,what="mu"))
phi1 <- as.numeric(stats::fitted(fit,what="sigma"))
nu1 <- as.numeric(stats::fitted(fit,what="nu"))
ell_cpp <- sum(tSB_moments_cpp(datos$z,datos$m,mu1,phi1,nu1,
                              GL64$nodes,GL64$weights)$log_prob)
cat("E. Converged:",fit$converged,"iterations:",fit$iter,
    "LogLik:",-fit$G.deviance/2,
    "error C++:",abs(ell_cpp+fit$G.deviance/2),"\n")
stopifnot(is.finite(ell_cpp),abs(ell_cpp+fit$G.deviance/2)<1e-6)
if (!isTRUE(fit$converged))
  warning("The simulated data fit did not converge within 50 cycles; check the fit.")

# F. Quantile residuals in the interval corresponding to each count.
set.seed(5147)
r <- as.numeric(stats::residuals(fit))
a <- ptSB(datos$z-1,mu=mu1,sigma=phi1,nu=nu1,bd=datos$m,regla=GL64)
b <- ptSB(datos$z,mu=mu1,sigma=phi1,nu=nu1,bd=datos$m,regla=GL64)
u <- stats::pnorm(r)
stopifnot(length(r)==n,all(is.finite(r)),
          all(u>=a-1e-10 & u<=b+1e-10))
cat("F. RQR residuals: OK\n")
cat("EXTERNAL VALIDATION COMPLETED\n")
