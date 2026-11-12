# Stage 3 — Author-reported restoration and functional validation

**Date:** 2026-11-11  
**System:** Windows 10 x64, R 4.2.3, Rtools42; project copied to a short-path project directory (no particular local directory is required)  
**Scope:** A new project directory on the *same host*, rather than a second computer.

## Restoring package dependencies

After the first attempt from a deeply nested Windows location failed during Rcpp archive extraction, the author moved the project to a short path. The normalized 34-package lockfile was restored successfully by:

```r
source("renv/activate.R")
renv::restore(prompt = FALSE)
renv::status()
```

The final status was `No issues found -- the project is in a consistent state.` The restoration installed/cached Rcpp 1.0.10, gamlss 5.4-22, gamlss.dist 6.1-1, and the recorded reporting dependencies. A consistent lockfile does not itself establish that the statistical computations are correct.

## Functional comparisons with Stage 1

The author then executed:

```r
source("run_all_tests.R", encoding = "UTF-8")
source("examples/fit_population_gamlss.R", encoding = "UTF-8")
rmarkdown::render("reproducibility_report.Rmd", output_format = "html_document")
```

| Check | Author-reported result | Comparison with Stage 1 |
|---|---:|---|
| Interior PMF normalization error | `3.441691376e-15` | Matches |
| Complementary-tail consistency error | `0` | Matches |
| CDF cumulative-sum error | `3.552713679e-15` | Matches |
| Boundary mass error | `3.774758284e-15` | Matches |
| Boundary expected-score error | `5.710281315e-12` | Matches |
| Mixed-regression gradient error | `7.01216862353e-11` | Matches |
| SB / tSB / tSB0 region count | `1 / 7 / 1` | Matches |
| Population weighted support | 324 rows, total weight 4500 | Matches |
| Population maximum mass error | `3.942868e-11` | Matches |
| Population gradient max. abs. | `1.287956e-07` | Matches |
| GAMLSS fit | `converged = TRUE`, 5 iterations | Matches |
| GAMLSS/C++ log-likelihood difference | `0` | Matches |
| Final gradient max. abs. | `0.003989652` | Matches |
| Report rendering | HTML created | Success |

The parameter estimates matched all Stage 1 digits displayed:

| Parameter | Fitted estimate |
|---|---:|
| `beta0` | `-0.299999987378` |
| `beta1` | `0.599999996817` |
| `alpha0` | `-0.105352125238` |
| `gamma0` | `0.336425425278` |
| `gamma1` | `0.900010722174` |

Elapsed fit time changed from 7.29 to 7.19 seconds; runtime is not a reproducibility criterion. The final gradient maximum (`0.003989652`) should not be interpreted as an extremely strict stationarity certificate.

## Supplemental acceptance tests: reported outcomes (2026-11-11)

After the standalone distribution test initially stopped at its expected-score check under GL64, the author compared GL64 and GL96 on five parameter settings. The expected-information call **in the supplemental test only** was switched to GL96, with the original criterion `max(abs(inf$expected_score)) < 1e-6` retained.

The author then reported successfully executing:

```r
source("tests/test_tsb_distribution.R", encoding = "UTF-8")
source("tests/test_tau_limits.R", encoding = "UTF-8")
```

| Test section | Author-reported outcome |
|---|---|
| A. Discrete mass, CDF, tail and quadrature checks | Mass error `3.44169137634e-15`; CDF error `3.55271367880e-15`; GL64/GL96 comparison `2.91433543964e-16` |
| B. Moments and expected information | `OK` with GL96 in expected-information call |
| C. Extreme parameters | `OK` |
| D. Simulated mean | Absolute discrepancy `0.02445`; script's acceptance checks passed (single simulated run, not an estimator-bias assessment) |
| E. GAMLSS model | `Converged: TRUE`, 3 iterations, log-likelihood `-617.4722`, GAMLSS/C++ error `0` |
| F. Randomized quantile residuals | `OK` |
| Numerical tau switching | `TAU BOUNDARIES VALIDATED; GAMLSS FAMILY UNCHANGED` |

At `eta_tau` of `-Inf`, `-40`, and `-35`, the piecewise evaluator selected `SB`; at `0`, `tSB`; and at `35`, `40`, and `Inf`, `tSB0`. The printed shape information was zero by construction in the switched limiting branches, and `0.001086013` at the tested interior point.

These results establish that the documented acceptance checks passed **on the author's tested setup**. They are not uniform numerical error bounds, cross-operating-system validation, or a Monte Carlo robustness study. The primary fit's final absolute gradient (`0.003989652`) also remains a numerical diagnostic requiring appropriate interpretation.

## Remaining prepublication checks

- Record `source("docs/record_environment.R", encoding = "UTF-8")` from the accepted short-path environment, if an environment log is desired for the release record.
- Review the repository inventory for accidental local paths, generated artifacts, credentials, and research notebooks not intended for distribution.
- MIT license and `CITATION.cff` have been added; reconfirm software attribution, ownership and provenance before GitHub publication.
- Cross-computer restoration and broader robustness simulations may be useful additional work, but have **not** been carried out in these acceptance checks.

All observations here come from console outputs supplied by the author, not from independent execution by the documentation editor.

## Supplemental distribution test: quadrature sensitivity

During the subsequent Stage 3 checks, the author reported that
`source("tests/test_tsb_distribution.R")` stopped at section B with
`max(abs(inf$expected_score)) < 1e-6` when expected information used GL64.
The original Stage 1 standalone test had the same criterion; this is not a
translation-induced failure.

For the demanding configuration `(mu, phi, nu) = (0.2, 0.5, 0.1)`, the
author evaluated `tSB_information_cpp()` using both quadrature rules:

| Quadrature | Absolute mass error | Maximum absolute expected score |
|---|---:|---:|
| GL64 | 5.08622577389e-08 | 1.25032393670e-05 |
| GL96 | 3.49117401655e-11 | 9.62066669108e-11 |

GL96 also satisfied the expected-score criterion at the other four tested
parameter configurations. In this candidate, **only the expected-information
call within section B of the supplemental test** has been switched to GL96.
The C++ core, the model computations, and the test tolerance remain unchanged.
The author subsequently ran the patched test successfully through sections A–F and also passed the standalone tau-boundary test.
This comparison supports a numerical quadrature-sensitivity explanation,
not a general guarantee that GL96 is accurate for every parameter value.
