# Computational methodology — t-Simplex-Binomial

## Hierarchical model

The model is defined by

\[
Z_i\mid Y_i\sim\mathrm{Binomial}(m_i,Y_i),\qquad
Y_i\sim t\mathrm{Simplex}(\mu_i,\phi_i,\nu_i).
\]

Its marginal probability mass function is evaluated as

\[
P(Z_i=z_i)=\binom{m_i}{z_i}\int_0^1 y^{z_i}(1-y)^{m_i-z_i}
 f_{tS}(y;\mu_i,\phi_i,\nu_i)\,dy.
\]

The GAMLSS `sigma` argument represents the dispersion parameter \(\phi\). The implemented predictors are

\[
\operatorname{logit}(\mu_i)=x_i^\top\beta,\qquad
\log(\phi_i)=\alpha_0,\qquad
\log(\nu_i)=x_i^\top\gamma_\nu.
\]

The separate boundary-evaluation layer uses

\[
\tau_i=(1+\nu_i)^{-1},\qquad
\eta_{\tau_i}=\operatorname{logit}(\tau_i)=-\log\nu_i.
\]

When the design matrices are identical, \(\gamma_\tau=-\gamma_\nu\). To compute \(\nu\) at an interior point, `exp(-eta_tau)` avoids the numerical saturation that can arise from reconstructing it using a rounded `plogis(eta_tau)`.

## Numerical probabilities and derivatives

- `tsb_core.cpp` evaluates log-PMFs, the three marginal scores and curvature quantities using Gauss–Legendre quadrature and analytic Fisher/Louis identities; the main estimator does not rely on finite-difference derivatives.
- The default rule uses 64 nodes; `tSB_GL_rule(96L)` permits a sensitivity comparison.
- `tSB_information_cpp()` sums expected information across the finite count support; `tSB_tau_information()` uses score outer products for the piecewise evaluator.
- The CDF uses normalized quadrature weights and binomial CDFs; its normalization may differ slightly from the unnormalized C++ PMF. This discrepancy should be checked for each scenario.
- `rtSB()` simulates using the hierarchical construction; its output should be checked against the discrete mass function for new parameter ranges.

## Shape limits and switching rule

\[
\tau=0\Rightarrow\mathrm{SB},\qquad
0<\tau<1\Rightarrow\mathrm{tSB},\qquad
\tau=1\Rightarrow\mathrm{tSB}_0.
\]

For computational stability, the separate evaluator uses the limiting distributions when \(\eta_\tau\leq-35\) or \(\eta_\tau\geq35\). These are **numerical switching thresholds**, not mathematical equalities for finite shape predictors and not uniform error bounds. Checks in research Sections 5.152–5.153 were limited to the specified grids.

Within a switched region, the shape predictor is absent from the *piecewise approximated* log-likelihood; its score and expected information are zero by construction. The joins are not necessarily exactly smooth.

## Identification and expected information

\[
\mathcal I_{\gamma_\tau\gamma_\tau}=
\sum_i\mathcal I_{\eta_{\tau_i}}x_ix_i^\top.
\]

As the predictor moves to an extreme limit, its shape information tends to zero. Consequently, a numerically well-evaluated PMF and gradient do not guarantee strong identification. In particular, reliable information for a slope requires sufficient covariate variation in the informative interior region.

## Scope and limitations

- The GAMLSS constructor `tSB()` fits the **interior** family for \(\nu>0\), with Fisher curvature used as working information in RS; the boundary evaluator is a separate numerical component, **not** an automatic three-region GAMLSS estimator.
- Algorithmic convergence, final gradient norm, and observed-information definiteness must be assessed separately. No local diagnostic proves a global maximum-likelihood optimum.
- Earlier Newton comparisons in research Section 5.150 indicated a local interior maximum for the particular simulated sample and poorly conditioned positive-definite observed information; they are not general guarantees.
- The repository does not yet implement an independent Monte Carlo Wald-coverage study or exact three-region estimation in GAMLSS. Neither capability should be claimed.
- Report the R, `Rcpp`, `gamlss`, and `gamlss.dist` versions, the compiler, operating system and quadrature rule alongside results. Stage 1 numeric outputs are recorded in `REFERENCE_RESULTS.md`.

## Suggested methods description (to be adapted for the paper)

> The marginal probabilities of the t-Simplex-Binomial distribution were evaluated by Gauss–Legendre quadrature. An Rcpp implementation computed the log-probability mass function, analytic marginal scores, and expected information, with expected information obtained by summing over the finite count support. Regression models were fitted within the GAMLSS framework using logit and log links for the mean and the positive dispersion and shape parameters, respectively. For numerical analysis of limiting distributions, we additionally used the transformation \(\tau=(1+\nu)^{-1}\) and evaluated the Simplex-Binomial and tSB₀ limits at their boundaries. A separate, numerically validated switching rule was used for sufficiently extreme shape predictors; it was not included in the GAMLSS RS estimator.
