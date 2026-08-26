# =============================================================================
# SECTION 6: minimum sample size for the search-cost block.
#
# Established: the likelihood is correct (enumeration + brute force at J=3 and
# J=4, sums to 1 across the parameter space) and the estimator is consistent --
# a single 9600-task fit cut search-block bias from 0.446 to 0.098. This
# replicates that at both ends, because single-dataset fits on this block are
# noise and misled this project repeatedly.
#
# beta is unbiased at every size tested (t < 1). The question is only how many
# tasks the (kappa, sigma_xi, a0, gamma0) block needs.
# =============================================================================
suppressPackageStartupMessages(library(parallel))
NREP <- 12L; NCORE <- 12L; ND <- 60L
TRUTH <- c(1.0,-0.8,0.9,-0.6,log(0.7),0.4,log(0.35),0.12)
PN <- c("beta1","beta2","kappa1","kappa2","log_sig_xi","a0","gamma0","gamma_pos")

one <- function(r, N) {
  source("R/sift_model.R", local = TRUE); ginv <- make_ginv()
  d <- sim_sift(N = N, T = 8L, J = 6L, seed = 60000L + r, ginv = ginv)
  info <- prep_sift(d); dr <- make_draws_sift(d, nd = ND, seed = 70000L + r)
  p0 <- rep(0, 8); p0[5] <- log(0.7)
  f <- optim(p0, function(p) nll_sift_ghk(p, d, info, dr, ginv, nq = 12L),
             method = "BFGS", control = list(maxit = 400, reltol = 1e-8))
  c(f$par, conv = f$convergence)
}
cl <- makeCluster(NCORE); on.exit(stopCluster(cl)); clusterExport(cl, "ND")
out <- list()
for (N in c(300L, 1200L)) {                      # 2400 and 9600 tasks
  tasks <- N * 8L; t0 <- Sys.time()
  res <- do.call(rbind, parLapply(cl, seq_len(NREP), one, N = N))
  el <- as.numeric(Sys.time() - t0, units = "mins")
  est <- res[, 1:8]
  cat(sprintf("\n===== %d tasks, %d reps, %.1f min, %d converged =====\n",
              tasks, NREP, el, sum(res[, "conv"] == 0)))
  cat(sprintf("%-11s %7s %9s %8s %8s %9s\n", "param", "truth", "mean", "bias", "sd", "t(bias)"))
  for (k in 1:8) {
    m <- mean(est[, k]); s <- sd(est[, k]); b <- m - TRUTH[k]
    cat(sprintf("%-11s %7.3f %9.3f %+8.3f %8.3f %9.2f\n", PN[k], TRUTH[k], m, b, s, b/(s/sqrt(NREP))))
  }
  bb <- colMeans(est) - TRUTH
  cat(sprintf("  beta |bias| = %.3f    search-block |bias| = %.3f\n",
              mean(abs(bb[1:2])), mean(abs(bb[3:7]))))
  out[[as.character(tasks)]] <- list(tasks = tasks, est = est, bias = bb,
                                     sd = apply(est, 2, sd), minutes = el)
  saveRDS(list(out = out, truth = TRUTH, par_names = PN, nrep = NREP, nd = ND),
          "R/out/simstudy_size.rds")
}
