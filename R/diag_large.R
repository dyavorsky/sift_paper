# Consistency check: one LARGE dataset. If the estimator is consistent, a
# 9600-task fit must land near the truth. If it is still off by ~0.7 on gamma0,
# the model is misspecified relative to sim_sift().
source("R/sift_model.R"); ginv <- make_ginv()
PN <- c("beta1","beta2","kappa1","kappa2","log_sig_xi","a0","gamma0","gamma_pos")
d <- sim_sift(N = 1200L, T = 8L, J = 6L, seed = 31337L, ginv = ginv)
info <- prep_sift(d); dr <- make_draws_sift(d, nd = 60L, seed = 424L)
cat(sprintf("%d tasks, mean clicks %.2f\n", d$S, mean(rowSums(d$cons_set != 0))))
p0 <- rep(0, 8); p0[5] <- log(0.7)
t0 <- Sys.time()
f <- optim(p0, function(p) nll_sift_ghk(p, d, info, dr, ginv, nq = 12L),
           method = "BFGS", control = list(maxit = 400, reltol = 1e-8))
cat(sprintf("conv=%d  %.1f min\n\n", f$convergence, as.numeric(Sys.time()-t0, units="mins")))
for (k in 1:8) cat(sprintf("  %-11s true %7.3f  est %7.3f  bias %+7.3f\n",
                           PN[k], d$true_par[k], f$par[k], f$par[k]-d$true_par[k]))
cat(sprintf("\n  mean|bias| = %.3f   search-block |bias| = %.3f\n",
            mean(abs(f$par-d$true_par)), mean(abs((f$par-d$true_par)[3:7]))))
saveRDS(f, "R/out/diag_large.rds")
