# Previously observed reference results

The following values come from the author's reported RStudio runs or prior research notebook checks. Sections 5.145 and 5.150–5.156 **have not been independently rerun in full for this translated distribution**. The main selected checks and the Section 5.146 population fit were subsequently reproduced on the renamed candidate, and supplemental distribution and tau-limit tests were passed on the author's Windows setup. They are provided as regression-test references, not universal numerical guarantees.

| Research section | Check | Previously reported result |
|---|---|---|
| 5.145 | Standalone implementation and interior GAMLSS | Constant-shape GAMLSS fit converged in 4 iterations; log-likelihood approximately −690.6458. |
| 5.146 | Population-weighted scenario with variable shape | 5 iterations; parameters recovered; GAMLSS/C++ log-likelihood difference 0; maximum final absolute gradient approximately 0.003989652. |
| 5.150 | Newton refinement for one simulated sample of 800 observations | Maximum gradient approximately 3.62×10⁻¹¹; smallest observed-information eigenvalue approximately 0.02011198; coefficient discrepancy across two starts approximately 4.00×10⁻⁹; a local interior maximum. |
| 5.152 | Transition at ητ = ±35 | Maximum PMF error approximately 8.33×10⁻¹⁷ over the checked grid. |
| 5.153 | Log-PMF and score behavior near thresholds | Maximum log-PMF error approximately 1.28×10⁻¹³; maximum mean-score error approximately 2.49×10⁻¹³; maximum shape-score error approximately 1.24×10⁻¹³. |
| 5.154 | Three-region evaluator | Log-PMF equivalence approximately 4.44×10⁻¹⁶; Fisher equivalence approximately 1.09×10⁻¹⁶; shape information exactly zero in switched limits. |
| 5.155 | Separate implementation versus research notebook | Moments and information equivalence approximately 4.44×10⁻¹⁶; GAMLSS family retained. |
| 5.156 | Regression gradient across mixed regions | Maximum error approximately 7.01×10⁻¹¹; 1 SB, 7 tSB, and 1 tSB₀ regions; zero transformation error in the tested interior comparison. |

**Important:** Reported errors depend on the tested grids, quadrature rules and fitted scenarios. They are not uniform error bounds or evidence of global identification.

## Author's 2026-11-11 Stage 1 verification

- `run_all_tests.R`: interior mass error 3.441691376e-15; mixed-regression central-difference gradient error 7.01216862353e-11; all selected tests passed.
- `examples/fit_population_gamlss.R`: converged in 5 iterations; GAMLSS/C++ log-likelihood discrepancy 0; maximum final absolute gradient 0.003989652.
- `reproducibility_report.Rmd`: HTML report rendered successfully.

## Stage 3 follow-up

The author reported exact agreement at the displayed precision for the main selected test and the population-weighted fit, successful HTML rendering, and successful standalone distribution (A–F) and tau-limit tests on 2026-11-11. Full details and scope limitations: `STAGE3_VALIDATION.md`.

## Fresh-session checklist

- [x] Compile `tsb_core.cpp` with `Rcpp::sourceCpp()` (performed by the master test; author-reported).
- [x] Run `source("run_all_tests.R")` without errors (author-reported, translated Stage 2).
- [x] Run the weighted population example (author-reported, translated Stage 2).
- [ ] Record `sessionInfo()` and effective package versions.
- [x] Rerun `tests/test_tsb_distribution.R`, `tests/test_tau_limits.R`, and the master script following the documented changes (author-reported).
- [ ] Document the Gauss–Legendre rule and numerical boundary thresholds when publishing results.
