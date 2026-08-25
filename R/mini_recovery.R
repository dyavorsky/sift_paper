# Small parameter-recovery check. Deliberately undersized -- this asks whether
# the estimator moves to the truth, not whether it is efficient.
source("R/sift_model.R")
ginv <- make_ginv()
d    <- sim_sift(N = 60, T = 8, J = 5, seed = 11, ginv = ginv)
info <- prep_sift(d)
dr   <- make_draws_sift(d, nd = 120, seed = 4)
p0   <- rep(0, length(d$true_par)); p0[length(p0) - 1] <- log(0.5)
t0 <- Sys.time()
fit <- optim(p0, function(p) nll_sift_ghk(p, d, info, dr, ginv),
             method = "BFGS", control = list(maxit = 200, reltol = 1e-6))
el <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("\ntasks=%d  nd=120  conv=%d  nll=%.2f  %.1f min\n",
            d$S, fit$convergence, fit$value, el))
cat(sprintf("%-12s %8s %8s %8s\n", "param", "true", "est", "err"))
for (k in seq_along(d$true_par))
  cat(sprintf("%-12s %8.3f %8.3f %+8.3f\n", d$par_names[k], d$true_par[k],
              fit$par[k], fit$par[k] - d$true_par[k]))
cat(sprintf("\nRMSE = %.4f\n", sqrt(mean((fit$par - d$true_par)^2))))
saveRDS(fit, "R/out/mini_recovery.rds")
