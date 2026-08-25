# =============================================================================
# How does the PM likelihood scale? Needed to decide whether SIFT can be
# estimated inside the 12-hour budget.
#
# We time a single likelihood evaluation rather than a full fit, because the
# fit cost is (per-eval cost) x (number of evals), and the eval count is
# governed by the parameter count, not by N or J. Run solo -- contention
# roughly doubles everything.
# =============================================================================
source("R/chung_pm.R")
ginv <- make_ginv()

time_eval <- function(N, J, nd, reps = 5) {
  # sim_ssm's default alpha0 has J-1 = 4 entries; recycle it for other J
  dat  <- sim_ssm(N = N, J = J, seed = 1, ginv = ginv,
                  alpha0 = rep(c(-0.5, -0.8, 0.5, 0.8), length.out = J - 1L))
  info <- prep_data(dat); grp <- prep_groups(dat, info)
  dr   <- make_draws(N, J, nd = nd, seed = 7)
  p    <- dat$true_par
  invisible(nll_grouped(p, dat, info, dr, ginv, grp))          # warm up
  t0 <- Sys.time()
  for (i in seq_len(reps)) invisible(nll_grouped(p, dat, info, dr, ginv, grp))
  as.numeric(Sys.time() - t0, units = "secs") / reps
}

cat("=========================================================\n")
cat("A. Scaling in N (observations), J = 5, nd = 100\n")
cat("=========================================================\n")
cat(sprintf("  %8s %10s %12s\n", "N", "sec/eval", "sec per 1k"))
base <- NULL
for (N in c(500L, 1000L, 2000L, 4000L, 8000L)) {
  s <- time_eval(N, 5L, 100L)
  cat(sprintf("  %8d %10.3f %12.3f\n", N, s, 1000 * s / N))
}

cat("\n=========================================================\n")
cat("B. Scaling in J (alternatives per task), N = 1000, nd = 100\n")
cat("=========================================================\n")
cat(sprintf("  %8s %10s\n", "J", "sec/eval"))
for (J in c(3L, 5L, 8L, 12L, 16L)) {
  s <- time_eval(1000L, J, 100L)
  cat(sprintf("  %8d %10.3f\n", J, s))
}

cat("\n=========================================================\n")
cat("C. Scaling in nd (simulation draws), N = 1000, J = 5\n")
cat("=========================================================\n")
cat(sprintf("  %8s %10s\n", "nd", "sec/eval"))
for (nd in c(25L, 50L, 100L, 200L, 400L)) {
  s <- time_eval(1000L, 5L, nd)
  cat(sprintf("  %8d %10.3f\n", nd, s))
}
