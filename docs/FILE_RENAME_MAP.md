# English file-name mapping

Names were standardized to make the repository easier to navigate for international readers. The names of exported R/Rcpp functions (`dtSB`, `ptSB`, `qtSB`, `rtSB`, `tSB`, `tSB_moments_cpp`, etc.) and model parameters were **not** changed.

| Stage 1 path | Stage 2 path | Purpose |
|---|---|---|
| `documentacion/CHECKLIST_GITHUB.md` | `docs/GITHUB_CHECKLIST.md` | Auxiliary method, test, example, or documentation |
| `documentacion/METODOLOGIA.md` | `docs/METHODOLOGY.md` | Auxiliary method, test, example, or documentation |
| `documentacion/REGISTRO_ETAPA1.md` | `docs/STAGE1_VALIDATION.md` | Auxiliary method, test, example, or documentation |
| `documentacion/RESULTADOS_REFERENCIA.md` | `docs/REFERENCE_RESULTS.md` | Auxiliary method, test, example, or documentation |
| `documentacion/STAGE2_REPRODUCIBILITY.md` | `docs/STAGE2_REPRODUCIBILITY.md` | Auxiliary method, test, example, or documentation |
| `documentacion/TRAZABILIDAD.md` | `docs/TRACEABILITY.md` | Auxiliary method, test, example, or documentation |
| `documentacion/registrar_entorno.R` | `docs/record_environment.R` | Auxiliary method, test, example, or documentation |
| `ejemplos/ajuste_poblacional_GAMLSS.R` | `examples/fit_population_gamlss.R` | Auxiliary method, test, example, or documentation |
| `pruebas/validar_regresion_mixta.R` | `tests/test_mixed_regression.R` | Auxiliary method, test, example, or documentation |
| `reporte_reproducible.Rmd` | `reproducibility_report.Rmd` | Reproducibility file |
| `tSB.R` | `tsb_distribution.R` | Primary R distribution and GAMLSS implementation |
| `tSB.cpp` | `tsb_core.cpp` | Core numeric C++ implementation |
| `tSB_limites_tau.R` | `tsb_tau_limits.R` | Auxiliary method, test, example, or documentation |
| `validar_limites_tau.R` | `tests/test_tau_limits.R` | Auxiliary method, test, example, or documentation |
| `validar_tSB.R` | `tests/test_tsb_distribution.R` | Auxiliary method, test, example, or documentation |
| `validar_todo.R` | `run_all_tests.R` | Auxiliary method, test, example, or documentation |

`LEEME.md` was omitted because the consolidated English `README.md` supersedes its information. Historical hash records still contain Stage 1 path names for auditability. The root `renv/` bootstrap paths are unchanged because renv expects those conventional names. Only file paths and user-facing messages were translated; the statistical interfaces were kept stable.

The master command-line switch is now `--fit`; `--ajuste` remains accepted as an optional legacy alias. This option only determines whether the population-weighted example is run in the same Rscript process.
