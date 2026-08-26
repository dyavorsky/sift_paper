# Misspecification or finite-sample bias?
#
# On a dataset large enough that sampling noise is negligible, the log-likelihood
# must be maximised AT THE TRUTH if the likelihood matches the DGP. If it prefers
# the biased point instead, the likelihood is misspecified relative to sim_sift()
# and no amount of data will fix it.
#
# nd is set high so simulation (Jensen) bias, measured at 0.019/task at nd=40,
# cannot drive the comparison.
source("R/sift_model.R"); ginv <- make_ginv()
d <- sim_sift(N = 1000L, T = 8L, J = 6L, seed = 4242L, ginv = ginv)
info <- prep_sift(d); dr <- make_draws_sift(d, nd = 400L, seed = 777L)
cat(sprintf("dataset: %d tasks, mean clicks %.2f\n\n", d$S, mean(rowSums(d$cons_set != 0))))

truth  <- d$true_par
biased <- c(1.012, -0.799, 1.218, -0.784, -0.007, 0.132, -0.514, 0.091)  # 24-rep mean
nt <- nll_sift_ghk(truth,  d, info, dr, ginv, nq = 16L)
nb <- nll_sift_ghk(biased, d, info, dr, ginv, nq = 16L)
cat(sprintf("nll at truth            = %10.2f\n", nt))
cat(sprintf("nll at 24-rep mean est  = %10.2f\n", nb))
cat(sprintf("difference              = %+10.2f  (%s)\n\n", nb - nt,
            if (nb < nt) "LIKELIHOOD PREFERS THE BIASED POINT -> misspecified"
            else "truth wins -> finite-sample bias, not misspecification"))

cat("profile in gamma0, all other parameters held at truth:\n")
cat(sprintf("  %9s %12s\n", "gamma0", "nll"))
best <- Inf; arg <- NA
for (g in seq(-1.6, -0.2, by = 0.2)) {
  p <- truth; p[7] <- g
  v <- nll_sift_ghk(p, d, info, dr, ginv, nq = 16L)
  if (v < best) { best <- v; arg <- g }
  cat(sprintf("  %9.2f %12.2f%s\n", g, v, if (abs(g - truth[7]) < 1e-9) "   <- truth" else ""))
}
cat(sprintf("\nprofile argmin = %.2f   truth = %.2f\n", arg, truth[7]))
