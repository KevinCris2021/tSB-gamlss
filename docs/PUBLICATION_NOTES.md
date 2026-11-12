# Publication preparation notes

**Status:** Candidate repository bundle, not a GitHub release.

## Included release metadata

- `LICENSE`: MIT License with the 2026 copyright notice for Kevin Cristopher Gargate Osorio.
- `AUTHORS.md`: software authorship and academic acknowledgment kept separate.
- `CITATION.cff`: CFF 1.2.0 software metadata (no fabricated DOI, ORCID, release date, version, or repository URL).

## Scope of the verified work

The author reported successful `renv::restore()` in another short-path project copy on the same Windows computer, followed by the master test suite, population-weighted GAMLSS fit, HTML rendering, supplemental distribution checks A–F and tau-boundary tests. The expected-information supplemental check uses GL96 in the challenging case; its tolerance was not relaxed. No broader parameter-space accuracy guarantee, different-host portability claim or Monte Carlo robustness assessment is made.

## Actions immediately before first push

1. Confirm spelling and copyright ownership and whether any additional software contributors must be named. Academic advising is not automatically software co-authorship.
2. Review any borrowed code, examples, or external source snippets for compatible terms. R packages listed as dependencies do not become MIT-licensed by association.
3. Upload the contents of this ZIP as files within the repository root, including dotfiles `.Rprofile`, `.gitattributes`, `.gitignore` and `renv/.gitignore`; do **not** publish `renv/library` or compiled artifacts.
4. After the repository exists, update `CITATION.cff` with the actual repository URL and, after a release, its version and release date. Add a DOI only if one has been assigned.
5. Confirm that GitHub recognizes the MIT license and offers `Cite this repository`; inspect the README rendering and run the documented checks from an extracted release copy.

Neither copyright ownership verification nor acceptance on additional operating systems can be established from this bundle alone.
