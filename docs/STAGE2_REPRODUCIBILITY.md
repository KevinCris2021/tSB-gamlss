# Stage 2 — English documentation and reproducible R environment

## Observed environment and scope

The author reported successful Stage 1 checks on **2026-11-11** using R 4.2.3 (Windows 10 x64, UCRT), Rtools42 (`make`, `g++`, `gcc` detected), Rcpp 1.0.10, gamlss 5.4-22, gamlss.dist 6.1-1, rmarkdown 2.32, knitr 1.52, and renv 1.0.5. The default validation uses a 64-node Gauss–Legendre rule (GL64).

The author successfully created a 34-package `renv.lock` and reported `renv::status()` as consistent **before** this editorial restructuring. This ZIP contains that real lockfile and its activation files, but **does not** contain the local package library, which must not be committed.

## Lockfile portability correction

The original lockfile in the author-provided ZIP contained 29 records with `Source: Local` and machine-specific `RemoteUrl` values pointing to `C:/Users/.../win-library/4.2/`. For this repository candidate, the 29 affected entries, all marked `Repository: CRAN` in the original file, were normalized to `Source: Repository` and their `RemoteType` / `RemoteUrl` fields were removed. **The 34 package names, versions, dependencies, hashes and R version are retained.** These edits are intended to support CRAN/Archive restoration and to avoid exposing private machine paths; they do not certify that historical package versions can actually be retrieved and compiled. A subsequent author-run `renv::restore()` succeeded in a separate short-path project copy on the same Windows machine; see `STAGE3_VALIDATION.md`. The earlier note described the initial state of the candidate.

## Files needed for renv

- `.Rprofile` — project startup bootstrap;
- `renv/activate.R`, `renv/settings.json`, `renv/.gitignore` — renv initialization metadata;
- `renv.lock` — exact packages recorded (34 entries).

Never publish `renv/library/`, `renv/staging/`, backup folders, local caches, temporary files, or private absolute paths.

## Stage 3 independent validation protocol

Use a **separate extracted copy** of the ZIP, without inheriting any library from the author's existing project. Start a new R session in that copy and execute:

```r
source("renv/activate.R")
renv::restore(prompt = FALSE)     # may download and compile packages
renv::status()
source("docs/record_environment.R", encoding = "UTF-8")
source("run_all_tests.R", encoding = "UTF-8")
source("examples/fit_population_gamlss.R", encoding = "UTF-8")
rmarkdown::render("reproducibility_report.Rmd", output_format = "html_document")
```

Alternatively, use `Rscript run_all_tests.R --fit` from the project root. Record `R.version`, compiler configuration, package provenance, and numerical outputs. Compare against `REFERENCE_RESULTS.md` allowing translated output headings but **not** numerical changes. The author subsequently completed the master checks, weighted population fit, and HTML render in the renamed candidate on 2026-11-11. The two supplemental standalone tests were subsequently reported as passed (including a GL96 expected-information check); cross-computer portability remains unverified. See `STAGE3_VALIDATION.md`.

## Lockfile cautions

`renv::status()` is not a substitute for `renv::restore()` in a clean installation. If a historical CRAN binary is unavailable, a source build may require compatible build tools; do not silently replace the recorded package version. If the normalized lockfile does not restore, resolve that issue transparently before publication.

## Model-specific limitations

The tau boundary layer is an auxiliary piecewise evaluator, not an exact three-region RS/GAMLSS estimator. The population-weighted GAMLSS example is not a Monte Carlo robustness assessment. A `converged = TRUE` flag must be evaluated alongside gradient and likelihood diagnostics.
