# Prepublication checklist

## Stage 1 — Organization and original validation

- [x] Preserve original research implementation and its Stage 1 SHA-256 inventory.
- [x] Author reported successful master, population-weighted GAMLSS and HTML checks on 2026-11-11.

## Stage 2 — English documentation and renv

- [x] Translate research code comments, error messages, report and documentation to English (editorial candidate).
- [x] Rename repository files and folders with descriptive English names; update `source()` calls and file references.
- [x] Include real `renv.lock` with 34 package records, `.Rprofile` and `renv/` activation files.
- [x] Omit `renv/library/`, caches and backups from the repository candidate.
- [x] Normalize the 29 local-path lockfile records to repository-based CRAN records without changing versions.
- [x] Restore the portable lockfile in a separate directory on the same Windows computer.
- [ ] Validate cross-computer restoration (not yet performed).
- [ ] Confirm no remaining accidentally machine-specific file paths or generated binaries.

## Stage 3 — Acceptance tests (author-reported, passed on tested setup)

- [x] Open a fresh R session in a *new copy* of the project and run `renv::restore()` (same machine, short path).
- [x] Confirm `renv::status()` reports a consistent state after restoration.
- [ ] Capture a post-restoration package/toolchain log for the release record.
- [x] Execute `source("run_all_tests.R", encoding="UTF-8")` (author-reported).
- [x] Execute `source("examples/fit_population_gamlss.R", encoding="UTF-8")` (author-reported).
- [x] Render `reproducibility_report.Rmd` to HTML (author-reported).
- [x] Re-run the standalone distribution test (section B uses GL96 for expected information); sections A–F passed (author-reported, 2026-11-11).
- [x] Run the standalone tau-boundary test; SB, tSB, and tSB0 regions validated (author-reported, 2026-11-11).
- [x] Compare displayed numerical outcomes of the main tests and population fit with Stage 1.
- [x] Recheck the translated main test after wording-only fixes; master test passed (author-reported).

## Stage 4 — GitHub (publication candidate; not yet published)

- [x] Add MIT `LICENSE` with 2026 copyright holder; confirm authority to license all included source before publication.
- [x] Add software citation metadata in root `CITATION.cff` (no unverified DOI, ORCID, date, URL, or release version).
- [x] Add `AUTHORS.md` distinguishing software author from academic guidance.
- [x] Scan the ZIP for bundled libraries, credentials, generated binaries and personal absolute paths (no sensitive path discovered in source/config; remaining short-path examples are documentation only).
- [ ] Confirm copyright scope, other contributors' roles, and any included third-party source before publication.
- [ ] Confirm the author's full legal name as used in `LICENSE` and `CITATION.cff`.
- [ ] Create the public GitHub repository and verify rendered license and citation sidebar before tagging a version.
- [ ] Publish only the accepted source, tests, docs, license, attribution, citation, renv lockfile/bootstrap and checksum manifests.
- [ ] After first release, update `CITATION.cff` with a verified repository URL, version, release date and DOI (if assigned).

**Neither a full Monte Carlo robustness study nor a complete scientific article is included in this candidate.**
