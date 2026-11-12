# t-Simplex-Binomial (tSB): R, Rcpp, and GAMLSS implementation

**Status:** Stage 3 functional validation accepted on the author's tested Windows/R environment (author-reported, 2026-11-11). The restored environment, main checks, population fit, HTML report, supplemental distribution test (A–F), and tau-boundary test were exercised across short-path project copies on the same computer. Cross-computer portability and comprehensive robustness evaluation have not been established. This is research code, not a formal R package or CRAN release.

This repository contains numerical tools for the t-Simplex-Binomial distribution, including its probability mass function, cumulative distribution, quantiles, random generation, analytic score calculations, expected information, a three-parameter interior GAMLSS family, and a separate numerical layer for limiting distributions.

## Scientific and computational scope

The interior model uses a positive shape parameter `nu`, with `nu.link = "log"` in the GAMLSS family. The auxiliary boundary evaluator instead uses `tau = 1/(1 + nu)` and `eta_tau = -log(nu)`. With identical model matrices, coefficients on these two predictor scales have opposite signs. The auxiliary boundary evaluator is **not** a GAMLSS regression estimator that automatically handles all three regions.

The switching thresholds `eta_tau <= -35` (SB limit) and `eta_tau >= 35` (tSB0 limit) are **numerical approximations** assessed on specified test grids, not exact equalities for finite predictors or uniform error bounds. Shape information vanishes in the limiting regions under the piecewise approximation.

## Requirements

- R with a working C++11 toolchain for `Rcpp::sourceCpp()`; Windows users may need an Rtools release compatible with their installed R version.
- Required R packages: `Rcpp`, `gamlss`, `gamlss.dist`.
- For HTML rendering: `rmarkdown`, `knitr` and a working Pandoc installation (typically bundled with RStudio).
- Author-reported test environment: R 4.2.3, Windows 10 x64, Rtools42; Rcpp 1.0.10, gamlss 5.4-22, gamlss.dist 6.1-1, rmarkdown 2.32, knitr 1.52, and renv 1.0.5. The project includes the author-provided 34-package `renv.lock` with CRAN repository-oriented records. The author successfully restored it in a separate project directory on the same Windows machine; cross-machine restoration has not been tested.

The recommended setup is the `renv` restoration shown below. The following direct installation command is a fallback for inspecting the code outside an `renv` project and **does not reproduce the locked versions**:

```r
required <- c("Rcpp", "gamlss", "gamlss.dist", "rmarkdown", "knitr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
```

## Reproducible dependency setup

The repository contains `.Rprofile`, `renv/activate.R`, `renv/settings.json`, and `renv.lock` (34 package records). **The package library under `renv/library/` is intentionally excluded.** Open a fresh R session in the project root and, when validating restoration in an isolated copy, run:

```r
source("renv/activate.R")
renv::restore(prompt = FALSE)
renv::status()
```

Restoration can require internet access and a suitable R compiler/toolchain. Do not run `renv::restore()` in the already validated Stage 1 folder. Re-check package sources, hashes and versions after restoring because the original lockfile was recorded from local package directories.

## Run the validations

Start a fresh R session and set the working directory to the repository root (the directory containing `tsb_core.cpp`). Then run the following commands **in order**:

```r
source("run_all_tests.R", encoding = "UTF-8")
source("examples/fit_population_gamlss.R", encoding = "UTF-8")
rmarkdown::render("reproducibility_report.Rmd", output_format = "html_document")
```

Outcomes confirmed by the author in the translated project on 2026-11-11:

1. The translated master test printed `ALL SELECTED TESTS COMPLETED SUCCESSFULLY`.
2. The population-weighted GAMLSS fit reports `Converged: TRUE` after 5 iterations, with matching GAMLSS and C++ log-likelihoods in the tested scenario.
3. The report is rendered to `reproducibility_report.html`.
4. The supplemental `tests/test_tsb_distribution.R` passed sections A–F after the expected-information check in section B switched to GL96 (test-only quadrature refinement; unchanged tolerance and statistical core).
5. The supplemental `tests/test_tau_limits.R` validated the numerical SB/tSB/tSB0 switching regions.

Command-line alternative (with Rscript configured):

```sh
Rscript run_all_tests.R
Rscript run_all_tests.R --fit
```

To capture versions and compiler-related information, execute:

```r
source("docs/record_environment.R", encoding = "UTF-8")
```

The generated `sessionInfo()` record is an environment report, **not** a replacement for a dependency lockfile. Do not commit a machine-specific record without reviewing it first.

## File layout

| Path | Purpose |
| --- | --- |
| `tsb_core.cpp` | Rcpp computations of log-probabilities, scores, and expected information |
| `tsb_distribution.R` | Distribution functions, Gauss–Legendre quadrature, and interior GAMLSS family |
| `tsb_tau_limits.R` | Auxiliary evaluation of SB/tSB/tSB0 regions using `tau` |
| `run_all_tests.R` | Main independent validation driver |
| `tests/test_tsb_distribution.R`, `tests/test_tau_limits.R` | Preserved legacy validation scripts |
| `tests/test_mixed_regression.R` | Mixed-region regression gradient and information checks |
| `examples/fit_population_gamlss.R` | Population-weighted GAMLSS example with variable `nu` |
| `reproducibility_report.Rmd` | English HTML validation report |
| `docs/` | Methods, reference outcomes, provenance, and stage checklists |
| `docs/STAGE2_REPRODUCIBILITY.md` | Environment and lockfile-portability documentation |
| `docs/STAGE3_VALIDATION.md` | Author-reported restoration and numerical comparisons in the renamed project |
| `LICENSE`, `AUTHORS.md`, `CITATION.cff` | MIT licensing, attribution, and software citation metadata |
| `SHA256SUMS_STAGE1.txt` | Historical Stage 1 hashes with original Spanish file paths (provenance only) |
| `SHA256SUMS_REPOSITORY.txt` | SHA-256 hashes of the current candidate (excluding the manifest itself) |

## Test evidence and limitations

Stage 1 evidence reported by the author includes: interior probability-mass error `3.441691376e-15`; mixed-regression central-difference gradient error `7.01216862353e-11`; 1 SB, 7 tSB, and 1 tSB0 regions; convergence of the population-weighted fit in 5 iterations; final maximum absolute gradient approximately `0.003989652`. The latter must not be confused with a strict numerical stationarity certificate.

This evidence is specific to the tested grid, model settings, quadrature rule, and working environment. The author restored the lockfile in a distinct project copy on the same host; this is not a cross-computer portability test or a full Monte Carlo robustness assessment. See `docs/STAGE3_VALIDATION.md`.

## Reproducibility and release status

This Stage 3 acceptance candidate translates source comments, console diagnostics, and documentation and uses descriptive English file paths. Public function and argument names remain unchanged. Mathematical formulas, numerical routines, parameter names, and computational constants were not intentionally altered. **The renamed/translated source is not byte-identical to the Stage 1 baseline. The author reported successful master tests, population-fit comparison, HTML rendering, and both supplemental standalone tests.** The only post-Stage-1 adjustment to a numerical test is to use GL96 instead of GL64 in the expected-information check for an especially demanding parameter configuration; the acceptance tolerance remains `1e-6`. The statistical computation functions were not edited for this correction.

The author generated a 34-package `renv.lock` and reported a consistent `renv::status()` result before editorial changes. The present archive includes the author-provided lockfile plus `.Rprofile` and `renv/` activation files. Because the original lockfile stored personal Windows installation paths for 29 CRAN packages, the repository candidate removes those local paths and changes their source records to CRAN repository records without changing versions. **The author subsequently restored this normalized lockfile and obtained a consistent `renv::status()` in a separate short-path project copy on the same computer.** See `docs/STAGE3_VALIDATION.md`.

## Authorship, software license, and citation

**Software author and maintainer:** Kevin Cristopher Gargate Osorio. Academic guidance is acknowledged separately in [`AUTHORS.md`](AUTHORS.md); research-paper authorship must be determined independently.

The original repository source and documentation are released under the **MIT License**; see [`LICENSE`](LICENSE). The licenses of external dependencies (including Rcpp and GAMLSS) are not changed by this repository license. Before the first public release, confirm that any incorporated third-party source is cleared for distribution under the applicable terms.

For software citation, see [`CITATION.cff`](CITATION.cff). GitHub can display citation formats from this file after publication. No DOI, ORCID, released version, repository URL, or article citation has been assigned or invented in this candidate. Once a public repository and first release exist, update the citation metadata with their verified details.

**Release note:** This is a publication candidate, not a published GitHub release. Functional checks were reported as passed in the stated Windows/R environment; cross-computer portability and a full robustness simulation study remain outside the validated scope.
