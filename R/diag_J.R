# The homogeneous MLE recovered cleanly at J=6 but was biased at J=5 in two
# earlier runs. Isolate J, holding everything else fixed.
source("R/sift_model.R"); ginv <- make_ginv()
PN <- c("beta1","beta2","kappa1","kappa2","log_sig_xi","a0","gamma0","gamma_pos")
for (J in c(5L, 6L, 7L)) {
  d <- sim_sift(N=150L, T=8L, J=J, seed=21L, ginv=ginv)
  info <- prep_sift(d); dr <- make_draws_sift(d, nd=100L, seed=5L)
  K <- rowSums(d$cons_set != 0)
  p0 <- rep(0, 8); p0[5] <- log(0.7)
  f <- optim(p0, function(p) nll_sift_ghk(p, d, info, dr, ginv, nq=12L),
             method="BFGS", control=list(maxit=300, reltol=1e-6))
  b <- f$par - d$true_par
  cat(sprintf("\n--- J = %d  (%d tasks, mean clicks %.2f, zero-click %.0f%%) ---\n",
              J, d$S, mean(K), 100*mean(K==0)))
  for (k in 1:8) cat(sprintf("  %-11s true %7.3f  est %7.3f  bias %+7.3f\n",
                             PN[k], d$true_par[k], f$par[k], b[k]))
  cat(sprintf("  mean|bias| = %.3f\n", mean(abs(b))))
}
