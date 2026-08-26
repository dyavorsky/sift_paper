# =============================================================================
# Sample size x simulation draws, factorial, common seeds.
#
# Three prior results are not monotone in either dimension:
#   1200 tasks / nd=100  -> search-block |bias| 0.446
#   2400 tasks / nd=60   -> 0.070   (clean, all |t| < 2)
#   9600 tasks / nd=60   -> 0.149   (beta1 SD grows 4x)
#
# S and nd were varied together, so they cannot be separated from those runs.
# SML theory says the simulated-likelihood bias is O(1/nd) per observation and
# does NOT vanish as S grows at fixed nd, so nd must scale with S. This tests
# that directly, with the SAME data seeds across nd so the only thing changing
# within a row is the number of draws.
# =============================================================================
suppressPackageStartupMessages(library(parallel))
NREP <- 8L; NCORE <- 12L
TRUTH <- c(1.0,-0.8,0.9,-0.6,log(0.7),0.4,log(0.35),0.12)
PN <- c("beta1","beta2","kappa1","kappa2","log_sig_xi","a0","gamma0","gamma_pos")
GRID <- expand.grid(N = c(150L, 300L), nd = c(30L, 60L, 120L))

one <- function(r, N, nd) {
  source("R/sift_model.R", local = TRUE); ginv <- make_ginv()
  d <- sim_sift(N = N, T = 8L, J = 6L, seed = 90000L + r, ginv = ginv)   # seed fixed across nd
  info <- prep_sift(d); dr <- make_draws_sift(d, nd = nd, seed = 95000L + r)
  p0 <- rep(0, 8); p0[5] <- log(0.7)
  f <- optim(p0, function(p) nll_sift_ghk(p, d, info, dr, ginv, nq = 12L),
             method = "BFGS", control = list(maxit = 400, reltol = 1e-8))
  c(f$par, conv = f$convergence)
}
cl <- makeCluster(NCORE); on.exit(stopCluster(cl))
rows <- list()
for (i in seq_len(nrow(GRID))) {
  N <- GRID$N[i]; nd <- GRID$nd[i]; tasks <- N * 8L
  t0 <- Sys.time()
  res <- do.call(rbind, parLapply(cl, seq_len(NREP), one, N = N, nd = nd))
  el <- as.numeric(Sys.time() - t0, units = "mins")
  est <- res[, 1:8]; bb <- colMeans(est) - TRUTH
  rows[[i]] <- data.frame(tasks = tasks, nd = nd,
    beta_bias   = mean(abs(bb[1:2])),
    search_bias = mean(abs(bb[3:7])),
    gamma0_bias = bb[7],
    beta1_sd    = sd(est[, 1]),
    gamma0_sd   = sd(est[, 7]),
    conv = sum(res[, "conv"] == 0), minutes = el)
  cat(sprintf("tasks=%5d nd=%3d | beta|b|=%.3f search|b|=%.3f gamma0=%+.3f | sd(b1)=%.3f sd(g0)=%.3f | %d/%d conv | %.1f min\n",
              tasks, nd, rows[[i]]$beta_bias, rows[[i]]$search_bias, rows[[i]]$gamma0_bias,
              rows[[i]]$beta1_sd, rows[[i]]$gamma0_sd, rows[[i]]$conv, NREP, el))
  saveRDS(do.call(rbind, rows), "R/out/simstudy_factorial.rds")
}
cat("\n===== factorial =====\n"); print(do.call(rbind, rows), digits = 3, row.names = FALSE)
