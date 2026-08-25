# =============================================================================
# Compare the R port's Monte Carlo results against the PM column of
# Chung, Chintagunta & Misra (2025) Table 1.
#
# Caveat on comparability: their Table 1 DGP carries a random coefficient on
# the first X (log sigma_theta1 = -0.916). The shipped Matlab -- and therefore
# this port -- is the homogeneous version, so we compare only the eight shared
# parameters. Their reported SDs are across 100 datasets at nd = 100, N = 1000.
# =============================================================================

res <- readRDS("R/out/replication_raw.rds")

PAR <- c("alpha0_2", "alpha0_3", "alpha0_4", "alpha0_5",
         "alpha1_1", "alpha1_2", "og_E", "log_c_bar")
TRUE_VAL <- c(-0.5, -0.8, 0.5, 0.8, 1.5, -1, 2.5, log(0.5))

# Chung et al. Table 1, PM row (est and SD across their 100 datasets)
PAPER_EST <- c(-0.500, -0.799, 0.500, 0.796, 1.495, -0.999, 2.524, -0.682)
PAPER_SD  <- c( 0.062,  0.068, 0.053, 0.051, 0.063,  0.066, 0.162,  0.052)
PAPER_RMSE <- 0.040          # over their 9 parameters
PAPER_MIN  <- 0.48           # minutes per estimation, Matlab

for (md in c("zK", "znext")) {
  d <- res[res$mode == md & res$conv == 0, ]
  cat("\n=========================================================\n")
  cat(sprintf("Pr(F) = %-6s   %d converged replications\n", md, nrow(d)))
  cat("=========================================================\n")
  cat(sprintf("%-11s %7s | %8s %7s %7s | %8s %7s\n",
              "param", "true", "R mean", "R sd", "R bias", "paper", "papSD"))
  est <- as.matrix(d[, PAR])
  for (k in seq_along(PAR)) {
    m <- mean(est[, k]); s <- sd(est[, k])
    cat(sprintf("%-11s %7.3f | %8.3f %7.3f %+7.3f | %8.3f %7.3f\n",
                PAR[k], TRUE_VAL[k], m, s, m - TRUE_VAL[k],
                PAPER_EST[k], PAPER_SD[k]))
  }
  rmse_i <- sqrt(rowMeans((est - matrix(TRUE_VAL, nrow(est), 8, byrow = TRUE))^2))
  cat(sprintf("\n  mean per-dataset RMSE (8 params) : %.4f   [paper: %.3f over 9]\n",
              mean(rmse_i), PAPER_RMSE))
  cat(sprintf("  mean seconds per estimation      : %.1f    [paper: %.1f s, Matlab]\n",
              mean(d$secs), PAPER_MIN * 60))
  cat(sprintf("  mean searches per consumer       : %.3f\n", mean(d$K_mean)))

  # is the mean estimate distinguishable from truth?
  tstat <- (colMeans(est) - TRUE_VAL) / (apply(est, 2, sd) / sqrt(nrow(est)))
  cat(sprintf("  |t| for bias, max over params    : %.2f  (%s)\n", max(abs(tstat)),
              paste(PAR[which.max(abs(tstat))])))
}

# ---- what the shipped loop code's Pr(F) does to the estimates --------------
a <- res[res$mode == "zK"    & res$conv == 0, ]
b <- res[res$mode == "znext" & res$conv == 0, ]
common <- intersect(a$rep, b$rep)
a <- a[match(common, a$rep), ]; b <- b[match(common, b$rep), ]
cat("\n=========================================================\n")
cat(sprintf("Effect of the shipped loop code's Pr(F) bug (%d paired reps)\n", length(common)))
cat("=========================================================\n")
cat(sprintf("%-11s %10s %10s %10s\n", "param", "mean(zK)", "mean(znext)", "shift"))
for (k in seq_along(PAR)) {
  cat(sprintf("%-11s %10.4f %10.4f %+10.4f\n", PAR[k],
              mean(a[[PAR[k]]]), mean(b[[PAR[k]]]),
              mean(b[[PAR[k]]] - a[[PAR[k]]])))
}
ra <- sqrt(rowMeans((as.matrix(a[, PAR]) - matrix(TRUE_VAL, nrow(a), 8, byrow = TRUE))^2))
rb <- sqrt(rowMeans((as.matrix(b[, PAR]) - matrix(TRUE_VAL, nrow(b), 8, byrow = TRUE))^2))
cat(sprintf("\n  mean RMSE  zK = %.4f   znext = %.4f   (%+.1f%%)\n",
            mean(ra), mean(rb), 100 * (mean(rb) / mean(ra) - 1)))
