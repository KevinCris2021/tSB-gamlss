# Example validated in the research R Markdown, Section 5.146.
# Requires C++ and tsb_distribution.R to be loaded beforehand (via run_all_tests.R --fit).
# Uses weighted probabilities, NOT a random sample.

# ============================================================
# 1. Dependencies
# ============================================================

stopifnot(
  is.function(tSB),
  is.function(dtSB),
  is.function(tSB_moments_cpp),
  is.function(tSB_GL_rule)
)

regla_5146 <- tSB_GL_rule(64L)

# ============================================================
# 2. Define a known interior scenario
# ============================================================

parametros_5146 <- c(
  beta0  = -0.3,
  beta1  =  0.6,
  alpha0 = log(0.9),
  gamma0 = log(1.4),
  gamma1 = 0.9
)

niveles_x_5146 <- seq(
  -0.9, 0.9, length.out = 9L
)

m_5146 <- 35L

# Include every possible count at each level of x.
datos_5146 <- expand.grid(
  z = 0:m_5146,
  x = niveles_x_5146
)

datos_5146$m <- m_5146

datos_5146$mu_verdadera <- plogis(
  parametros_5146["beta0"] +
  parametros_5146["beta1"] * datos_5146$x
)

datos_5146$phi_verdadera <- exp(
  parametros_5146["alpha0"]
)

datos_5146$nu_verdadera <- exp(
  parametros_5146["gamma0"] +
  parametros_5146["gamma1"] * datos_5146$x
)

cat("========== TRUE PARAMETERS ==========\n")

print(parametros_5146, digits = 12)

cat("\nRange of nu:\n")
print(range(datos_5146$nu_verdadera))

# ============================================================
# 3. Construct the weighted discrete distribution
# ============================================================

datos_5146$probabilidad <- dtSB(
  x = datos_5146$z,
  mu = datos_5146$mu_verdadera,
  sigma = datos_5146$phi_verdadera,
  nu = datos_5146$nu_verdadera,
  bd = datos_5146$m,
  regla = regla_5146
)

# Normalize only the quadrature error
# within each covariate level.
datos_5146$grupo <- match(
  datos_5146$x,
  niveles_x_5146
)

masa_5146 <- as.numeric(rowsum(
  datos_5146$probabilidad,
  datos_5146$grupo,
  reorder = TRUE
))

error_masa_5146 <- max(abs(
  masa_5146 - 1
))

datos_5146$probabilidad <-
  datos_5146$probabilidad /
  masa_5146[datos_5146$grupo]

# Each level represents a weighted population
# with effective size 500.
datos_5146$w <- (
  500 * datos_5146$probabilidad
)

cat("\n========== WEIGHTED DISTRIBUTION ==========\n")

cat("Rows:", nrow(datos_5146), "\n")
cat("Total weight:", sum(datos_5146$w), "\n")
cat("Maximum mass error:", error_masa_5146, "\n")

stopifnot(
  error_masa_5146 < 1e-6,
  all(is.finite(datos_5146$w)),
  all(datos_5146$w > 0)
)

# ============================================================
# 4. Check the gradient at the known parameter values
# ============================================================

momentos_verdaderos_5146 <- tSB_moments_cpp(
  z = datos_5146$z,
  m = datos_5146$m,
  mu = datos_5146$mu_verdadera,
  phi = datos_5146$phi_verdadera,
  nu = datos_5146$nu_verdadera,
  nodes = regla_5146$nodes,
  weights = regla_5146$weights
)

X_5146 <- model.matrix(
  ~ x,
  data = datos_5146
)

w_5146 <- datos_5146$w

S_mu_5146 <- momentos_verdaderos_5146$score[, "mu"]
S_phi_5146 <- momentos_verdaderos_5146$score[, "phi"]
S_nu_5146 <- momentos_verdaderos_5146$score[, "nu"]

gradiente_verdadero_5146 <- c(
  as.numeric(crossprod(
    X_5146,
    w_5146 * S_mu_5146 *
      datos_5146$mu_verdadera *
      (1 - datos_5146$mu_verdadera)
  )),

  sum(
    w_5146 * S_phi_5146 *
      datos_5146$phi_verdadera
  ),

  as.numeric(crossprod(
    X_5146,
    w_5146 * S_nu_5146 *
      datos_5146$nu_verdadera
  ))
)

names(gradiente_verdadero_5146) <-
  names(parametros_5146)

cat("\n========== POPULATION GRADIENT ==========\n")

print(
  gradiente_verdadero_5146,
  digits = 12
)

cat(
  "Maximum absolute value:",
  max(abs(gradiente_verdadero_5146)),
  "\n"
)

# ============================================================
# 5. Fit GAMLSS from different interior starting values
# ============================================================

familia_5146 <- tSB(
  regla = regla_5146
)

error_ajuste_5146 <- NULL
advertencias_5146 <- character(0)

tiempo_5146 <- system.time({

  ajuste_5146 <- tryCatch(

    withCallingHandlers(

      gamlss::gamlss(
        formula = cbind(z, m - z) ~ x,
        sigma.formula = ~ 1,
        nu.formula = ~ x,

        family = familia_5146,
        data = datos_5146,
        weights = w,

        mu.start = plogis(
          -0.2 + 0.3 * datos_5146$x
        ),

        sigma.start = rep(
          0.75, nrow(datos_5146)
        ),

        nu.start = exp(
          log(1.3) + 0.2 * datos_5146$x
        ),

        control = gamlss::gamlss.control(
          n.cyc = 50,
          c.crit = 1e-6,
          trace = FALSE
        ),

        trace = FALSE
      ),

      warning = function(w) {
        advertencias_5146 <<- c(
          advertencias_5146,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),

    error = function(e) {
      error_ajuste_5146 <<- conditionMessage(e)
      NULL
    }
  )

})[["elapsed"]]

# ============================================================
# 6. Inspect the fitted model
# ============================================================

cat("\n========== GAMLSS FIT ==========\n")

cat("Elapsed time:", tiempo_5146, "seconds\n")

if (is.null(ajuste_5146)) {

  cat("Error:", error_ajuste_5146, "\n")

} else {

  cat(
    "Converged:",
    ajuste_5146$converged, "\n"
  )

  cat(
    "Iterations:",
    ajuste_5146$iter, "\n"
  )

  estimaciones_5146 <- c(
    as.numeric(coef(ajuste_5146, what = "mu")),
    as.numeric(coef(ajuste_5146, what = "sigma")),
    as.numeric(coef(ajuste_5146, what = "nu"))
  )

  names(estimaciones_5146) <-
    names(parametros_5146)

  tabla_estimaciones_5146 <- data.frame(
    Parameter = names(parametros_5146),
    True = unname(parametros_5146),
    Estimated = unname(estimaciones_5146),
    Difference = unname(
      estimaciones_5146 - parametros_5146
    )
  )

  cat("\n========== ESTIMATES ==========\n")

  print(
    tabla_estimaciones_5146,
    digits = 12,
    row.names = FALSE
  )

  # ==========================================================
  # 7. Check the log-likelihood and gradient
  # ==========================================================

  mu_hat <- as.numeric(
    fitted(ajuste_5146, what = "mu")
  )

  phi_hat <- as.numeric(
    fitted(ajuste_5146, what = "sigma")
  )

  nu_hat <- as.numeric(
    fitted(ajuste_5146, what = "nu")
  )

  momentos_hat <- tSB_moments_cpp(
    z = datos_5146$z,
    m = datos_5146$m,
    mu = mu_hat,
    phi = phi_hat,
    nu = nu_hat,
    nodes = regla_5146$nodes,
    weights = regla_5146$weights
  )

  ell_cpp <- sum(
    datos_5146$w * momentos_hat$log_prob
  )

  ell_gamlss <- -ajuste_5146$G.deviance / 2

  error_logLik <- abs(
    ell_cpp - ell_gamlss
  )

  gradiente_hat <- c(
    as.numeric(crossprod(
      X_5146,
      datos_5146$w *
        momentos_hat$score[, "mu"] *
        mu_hat * (1 - mu_hat)
    )),

    sum(
      datos_5146$w *
        momentos_hat$score[, "phi"] *
        phi_hat
    ),

    as.numeric(crossprod(
      X_5146,
      datos_5146$w *
        momentos_hat$score[, "nu"] *
        nu_hat
    ))
  )

  names(gradiente_hat) <-
    names(parametros_5146)

  cat("\n========== LOG-LIKELIHOOD ==========\n")

  cat("GAMLSS:", ell_gamlss, "\n")
  cat("C++:", ell_cpp, "\n")
  cat("Error:", error_logLik, "\n")

  cat("\n========== FINAL GRADIENT ==========\n")

  print(
    gradiente_hat,
    digits = 12
  )

  cat(
    "Maximum absolute gradient:",
    max(abs(gradiente_hat)), "\n"
  )

  # ==========================================================
  # 8. Final diagnostic
  # ==========================================================

  cat("\n========== DIAGNOSTIC ==========\n")

  if (
    isTRUE(ajuste_5146$converged) &&
    max(abs(gradiente_hat)) < 1e-2 &&
    error_logLik < 1e-5 &&
    max(abs(
      estimaciones_5146 - parametros_5146
    )) < 0.1
  ) {

    cat(
      "REGRESSION WITH VARIABLE NU VALIDATED ",
      "IN THE INTERIOR POPULATION SCENARIO\n"
    )

  } else {

    cat(
      "CHECK THE ESTIMATES OR THE ",
      "CONVERGENCE CRITERION\n"
    )
  }
}

if (length(advertencias_5146) > 0L) {

  cat("\n========== WARNINGS ==========\n")

  cat(
    paste(
      unique(advertencias_5146),
      collapse = "\n"
    ),
    "\n"
  )
}

cat("\nSECTION 5.146 COMPLETED\n")
