# Stage 1 — Validation record (author-reported)

On **2026-11-11**, the author ran the original Spanish-named scripts in RStudio on Windows 10 / R 4.2.3. These checks passed in the original pre-translation project. On 2026-11-11, the author subsequently reproduced the main test, weighted population example, and HTML report in the translated, renamed candidate; see `STAGE3_VALIDATION.md`.

- `validar_todo.R` (now `run_all_tests.R`): interior mass error **3.441691376e-15**; interior CDF error **3.552713679e-15**; boundary mass error **3.774758284e-15**; central-difference mixed-regression gradient error **7.01216862353e-11**. Mixed-region composition: 1 SB, 7 tSB, 1 tSB0.
- `ejemplos/ajuste_poblacional_GAMLSS.R` (now `examples/fit_population_gamlss.R`): 324 weighted support rows, total weight 4500, convergence TRUE in 5 iterations, matching log-likelihoods of approximately -13134.8 in GAMLSS and C++ with reported error 0, final maximum absolute gradient **0.003989652**.
- `reporte_reproducible.Rmd` (now `reproducibility_report.Rmd`): HTML rendering succeeded.

The maximum final gradient is not an arbitrarily strict stationarity certificate; population weights do not constitute a Monte Carlo experiment. Verify numerical reproducibility independently in Stage 3.
