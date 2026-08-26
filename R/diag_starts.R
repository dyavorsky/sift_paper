# Is the "bias" actually optimizer failure?
#
# BFGS from all-zeros on a weakly identified 8-parameter simulated likelihood
# may be landing in bad local optima. Per-replication estimates show huge spread
# and outright failures (beta1 = 0.47 in one 4800-task rep).
#
# Test: ONE dataset, many starting points. If starts disagree AND some reach a
# better nll than the cold start, the cold-start fits are not finding the MLE,
# and every "bias" number so far is contaminated.
source("R/sift_model.R"); ginv <- make_ginv()
PN <- c("beta1","beta2","kappa1","kappa2","log_sig_xi","a0","gamma0","gamma_pos")
d <- sim_sift(N = 150L, T = 8L, J = 6L, seed = 7001L, ginv = ginv)
info <- prep_sift(d); dr <- make_draws_sift(d, nd = 100L, seed = 8001L)
truth <- d$true_par

set.seed(11)
starts <- list(zeros = rep(0, 8))
starts$zeros[5] <- log(0.7)
starts$truth <- truth
for (j in 1:4) starts[[paste0("rand", j)]] <- truth + rnorm(8) * 0.6

cat(sprintf("dataset: %d tasks\n\n%-8s %10s %9s  %s\n", d$S, "start", "nll", "conv",
            paste(sprintf("%8s", PN), collapse = "")))
best <- Inf
for (nm in names(starts)) {
  f <- optim(starts[[nm]], function(p) nll_sift_ghk(p, d, info, dr, ginv, nq = 12L),
             method = "BFGS", control = list(maxit = 600, reltol = 1e-8))
  best <- min(best, f$value)
  cat(sprintf("%-8s %10.2f %9d  %s\n", nm, f$value, f$convergence,
              paste(sprintf("%8.3f", f$par), collapse = "")))
}
cat(sprintf("%-8s %10s %9s  %s\n", "TRUTH", sprintf("%.2f",
    nll_sift_ghk(truth, d, info, dr, ginv, nq = 12L)), "",
    paste(sprintf("%8.3f", truth), collapse = "")))
cat(sprintf("\nbest nll found = %.2f\n", best))
cat("If starts land at different nll, the surface is multimodal and cold-start\n")
cat("single fits are unreliable -> every bias number so far needs redoing.\n")
