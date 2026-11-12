# Report the actual R environment used to validate this repository.
# Run from the repository root, in a clean R session.
# This script does not install packages or modify validated model functions.

packages <- c("Rcpp", "gamlss", "gamlss.dist", "rmarkdown", "knitr", "renv")
cat("=== R SESSION ===\n")
print(utils::sessionInfo())
cat("\n=== PACKAGE VERSIONS ===\n")
for (pkg in packages) {
  version <- if (requireNamespace(pkg, quietly = TRUE))
    as.character(utils::packageVersion(pkg)) else "NOT INSTALLED"
  cat(sprintf("%s: %s\n", pkg, version))
}
cat("\n=== R BUILD DETAILS ===\n")
print(R.version[c("version.string", "platform", "arch", "os", "compiler")])
cat("\n=== C++ TOOLCHAIN ===\n")
print(tryCatch(system2(file.path(R.home("bin"), "R"),
                       c("CMD", "config", "CXX11"),
                       stdout = TRUE, stderr = TRUE),
               error = function(e) conditionMessage(e)))
cat("\nNumerical check default: Gauss-Legendre 64 nodes (GL64).\n")
