# Provenance and source traceability

The Stage 1 candidate derives from the research implementation through Section 5.156. Original named source files and historical SHA-256 hashes are preserved in the author's Stage 1 ZIP. `SHA256SUMS_STAGE1.txt` is an archival record **with the original file paths**, not a verification manifest for the renamed/translated source files.

## Relationship with the original implementation

- `tsb_core.cpp` corresponds to the original `tSB.cpp` numerical kernel.
- `tsb_distribution.R` corresponds to the original `tSB.R` interior GAMLSS family and distribution API.
- `tsb_tau_limits.R` corresponds to the original `tSB_limites_tau.R` auxiliary boundary evaluator.
- The former standalone tests and population GAMLSS example are located under `tests/` and `examples/` with English names.
- The exported functions and parameter names are retained; comments, user-facing messages, file-loading references and a legacy command-line flag were editorially adapted.

See `FILE_RENAME_MAP.md` for the complete old/new path mapping. The new `SHA256SUMS_REPOSITORY.txt` records the candidate's actual bytes.

## Reproducibility boundary

This candidate includes the author's authentic renv lockfile, normalized to avoid 29 machine-local CRAN source paths, and the renv bootstrap files. This does not assert that the packages have been restored successfully on a new computer. The author subsequently reported successful `renv::restore()`, main tests, population fit, HTML rendering, and both supplemental standalone tests on the same Windows machine; see `STAGE3_VALIDATION.md`. Cross-computer portability has not yet been verified.

The research notebook is not an executable dependency: it contains exploratory definitions that may collide with standalone functions. The publication candidate now includes an MIT `LICENSE`, an `AUTHORS.md` acknowledgment, and a root `CITATION.cff`. Confirm the attribution and authority to license included source before the first public release.
