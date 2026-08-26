# =============================================================================
# SECTION 6: how much data does each block of parameters need?
#
# Established by diagnostics (see R/DIAGNOSTICS.md):
#   * the likelihood is CORRECT -- on 8000 tasks the truth beats the biased
#     point by 29 nll units and the gamma0 profile peaks at the truth;
#   * what looked like a bug is finite-sample bias in the MLE, concentrated in
#     the (kappa, sigma_xi, a0, gamma0) block, which is weakly identified;
#   * beta is unbiased already at 1200 tasks (t = 0.8 and 0.1 over 24 reps).
#
# This measures the bias-vs-sample-size curve, which is the number a
# practitioner actually needs: how many tasks to field.
#
# Replication is not optional here. Single-dataset recovery checks on this block
# are meaningless -- the gamma0 error flips sign across draws.
# =============================================================================
suppressPackageStartupMessages(library(parallel))
SIZES <- list(c(N=150L,T=8L), c(N=300L,T=8L), c(N=600L,T=8L))   # 1200 / 2400 / 4800
NREP  <- 12L; NCORE <- 12L; ND <- 50L
TRUTH <- c(1.0,-0.8,0.9,-0.6,log(0.7),0.4,log(0.35),0.12)
PN <- c("beta1","beta2","kappa1","kappa2","log_sig_xi","a0","gamma0","gamma_pos")

one <- function(r, N, T) {
  source("R/sift_model.R", local = TRUE); ginv <- make_ginv()
  d <- sim_sift(N=N, T=T, J=6L, seed=7000L+r, ginv=ginv)
  info <- prep_sift(d); dr <- make_draws_sift(d, nd=ND, seed=8000L+r)
  p0 <- rep(0,8); p0[5] <- log(0.7)
  f <- optim(p0, function(p) nll_sift_ghk(p,d,info,dr,ginv,nq=12L),
             method="BFGS", control=list(maxit=300, reltol=1e-6))
  c(f$par, conv=f$convergence)
}
cl <- makeCluster(NCORE); on.exit(stopCluster(cl))
clusterExport(cl, c("ND"))
out <- list()
for (cfg in SIZES) {
  N <- cfg[["N"]]; T <- cfg[["T"]]; tasks <- N*T
  t0 <- Sys.time()
  res <- do.call(rbind, parLapply(cl, seq_len(NREP), one, N=N, T=T))
  el <- as.numeric(Sys.time()-t0, units="mins")
  est <- res[,1:8]
  cat(sprintf("\n===== %d tasks (N=%d, T=%d), %d reps, %.1f min =====\n", tasks, N, T, NREP, el))
  cat(sprintf("%-11s %7s %9s %8s %8s %9s\n","param","truth","mean","bias","sd","t(bias)"))
  tb <- numeric(8)
  for (k in 1:8) {
    m <- mean(est[,k]); s <- sd(est[,k]); b <- m-TRUTH[k]; tb[k] <- b/(s/sqrt(NREP))
    cat(sprintf("%-11s %7.3f %9.3f %+8.3f %8.3f %9.2f\n", PN[k], TRUTH[k], m, b, s, tb[k]))
  }
  out[[as.character(tasks)]] <- list(tasks=tasks, est=est, bias=colMeans(est)-TRUTH,
                                     sd=apply(est,2,sd), t=tb, minutes=el)
  saveRDS(list(out=out, truth=TRUTH, par_names=PN, nrep=NREP, nd=ND),
          "R/out/simstudy_scale.rds")
}
cat("\n===== bias vs sample size =====\n")
cat(sprintf("%8s %12s %14s %14s\n","tasks","|bias| beta","|bias| search","|bias| gamma0"))
for (nm in names(out)) {
  o <- out[[nm]]
  cat(sprintf("%8d %12.3f %14.3f %14.3f\n", o$tasks, mean(abs(o$bias[1:2])),
              mean(abs(o$bias[c(3,4,5,6,7)])), abs(o$bias[7])))
}
