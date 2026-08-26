# Is the ridge broken by a STRONGER position shifter?
#
# @Yavorsky_2021's whole argument is that an exogenous, alternative-specific
# search-cost shifter is what separates search costs from the post-search
# spread. His shifter (driving distance) spans 0.01 to 102 miles. Ours is
# position, and at gamma_pos = 0.12 over 5 slots the cost gradient is only
# 1.13x to 1.82x -- possibly too weak to identify the block.
#
# Unlike a field study, the strength of our shifter is a DESIGN CHOICE.
source("R/sift_model.R"); ginv <- make_ginv()
PN <- c("beta1","beta2","kappa1","kappa2","log_sig_xi","a0","gamma0","gamma_pos")
res <- list()
for (gp in c(0.12, 0.30, 0.60, 1.00)) {
  d <- sim_sift(N = 100L, T = 8L, J = 6L, gpos = gp, seed = 21L, ginv = ginv)
  info <- prep_sift(d); dr <- make_draws_sift(d, nd = 100L, seed = 5L)
  p0 <- rep(0, 8); p0[5] <- log(0.7)
  t0 <- Sys.time()
  f <- optim(p0, function(p) nll_sift_ghk(p, d, info, dr, ginv, nq = 12L),
             method = "BFGS", control = list(maxit = 300, reltol = 1e-6))
  el <- as.numeric(Sys.time() - t0, units = "mins")
  b <- f$par - d$true_par
  grad <- exp(gp * d$J) / exp(gp)          # cost ratio, worst slot vs best
  cat(sprintf("\n--- gamma_pos = %.2f  (cost gradient %.1fx across %d slots, %.1f min) ---\n",
              gp, grad, d$J, el))
  for (k in 1:8) cat(sprintf("  %-11s true %7.3f  est %7.3f  bias %+7.3f\n",
                             PN[k], d$true_par[k], f$par[k], b[k]))
  cat(sprintf("  mean|bias| = %.3f   |bias| on (kappa,a0,gamma0) = %.3f\n",
              mean(abs(b)), mean(abs(b[c(3,4,6,7)]))))
  res[[as.character(gp)]] <- list(gpos = gp, gradient = grad, bias = b,
                                  est = f$par, truth = d$true_par, minutes = el)
  saveRDS(res, "R/out/diag_shifter.rds")
}
cat("\n===== does a stronger shifter break the ridge? =====\n")
for (nm in names(res)) cat(sprintf("  gamma_pos %.2f  gradient %5.1fx   mean|bias| %.3f   search-block |bias| %.3f\n",
    res[[nm]]$gpos, res[[nm]]$gradient, mean(abs(res[[nm]]$bias)),
    mean(abs(res[[nm]]$bias[c(3,4,6,7)]))))
