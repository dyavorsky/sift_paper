# Bias or variance? Every recovery check so far has been a SINGLE dataset, and
# the gamma0 "bias" flips sign across runs (-0.03, -0.77, +0.64, +0.19). That is
# the signature of a flat ridge: the MLE lands at a sampling-noise-determined
# point along it. Only replication can separate bias from variance.
#
# Fits are independent, so parallelise over replications.
suppressPackageStartupMessages(library(parallel))
NREP <- 24L; NCORE <- 12L
one <- function(r) {
  source("R/sift_model.R", local = TRUE)
  ginv <- make_ginv()
  d <- sim_sift(N = 150L, T = 8L, J = 6L, seed = 5000L + r, ginv = ginv)
  info <- prep_sift(d); dr <- make_draws_sift(d, nd = 100L, seed = 9000L + r)
  p0 <- rep(0, 8); p0[5] <- log(0.7)
  f <- optim(p0, function(p) nll_sift_ghk(p, d, info, dr, ginv, nq = 12L),
             method = "BFGS", control = list(maxit = 300, reltol = 1e-6))
  c(f$par, conv = f$convergence, truth1 = d$true_par[1])
}
cl <- makeCluster(NCORE); on.exit(stopCluster(cl))
t0 <- Sys.time()
res <- do.call(rbind, parLapply(cl, seq_len(NREP), one))
cat(sprintf("%d reps on %d cores in %.1f min\n\n", NREP, NCORE,
            as.numeric(Sys.time() - t0, units = "mins")))
TRUTH <- c(1.0, -0.8, 0.9, -0.6, log(0.7), 0.4, log(0.35), 0.12)
PN <- c("beta1","beta2","kappa1","kappa2","log_sig_xi","a0","gamma0","gamma_pos")
est <- res[, 1:8]
cat(sprintf("%-11s %7s %9s %8s %8s %9s %8s\n",
            "param","truth","mean","bias","sd","t(bias)","|bias|/sd"))
for (k in 1:8) {
  m <- mean(est[, k]); s <- sd(est[, k]); b <- m - TRUTH[k]
  cat(sprintf("%-11s %7.3f %9.3f %+8.3f %8.3f %9.2f %8.2f\n",
              PN[k], TRUTH[k], m, b, s, b / (s / sqrt(NREP)), abs(b) / s))
}
cat(sprintf("\nconverged: %d/%d\n", sum(res[, "conv"] == 0), NREP))
saveRDS(list(est = est, truth = TRUTH, par_names = PN), "R/out/diag_reps.rds")
cat("\nt(bias) large  -> genuine bias.  sd large but t small -> pure variance.\n")
