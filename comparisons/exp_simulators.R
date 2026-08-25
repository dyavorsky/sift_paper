# =============================================================================
# COMPARISON: simulator choice -- Ursu, Seiler & Honka (2024) Table 2 on the
# SIFT design.
#
# Question. USH compare crude frequency, kernel-smoothed frequency, GHK and
# importance sampling for search models and report both accuracy and runtime.
# Their ranking is for a field design. Does it hold for SIFT, where position is
# randomized and the product page is designed?
#
# Design. One SIFT data set. Evaluate the log-likelihood at the TRUE parameters
# with each simulator, repeatedly, across independent draw sets. Two things
# matter and they are different: the SPREAD across draw sets (simulation noise,
# which is what the optimiser has to fight) and the LEVEL (bias -- crude
# frequency and kernel smoothing are not unbiased for the log-likelihood).
#
# Runtime is reported per evaluation because, as R/README.md documents, total
# fit time is (per-eval cost) x (eval count) and the eval count is governed by
# the parameter count rather than by the simulator.
# =============================================================================
source("comparisons/_common.R")

# ---- crude frequency: the unsmoothed accept-reject estimator ---------------
# Identical to nll_sift_ks but with a hard indicator instead of a logit kernel.
# It is unbiased for the PROBABILITY but discontinuous in the parameters, which
# is exactly why the literature smooths it.
nll_sift_cfs <- function(p, d, info, draws, ginv, nb = 2L, nk = 2L) {
  th <- unpack(p, nb, nk)
  sig_t <- sqrt(sum(th$kappa^2) + th$sig_xi^2)
  Xb <- as.vector(d$X %*% th$beta); mL <- as.vector(d$L %*% th$kappa)
  cst <- exp(th$gamma0 + th$gpos * d$pos); nd <- ncol(draws$H); u0 <- th$a0
  ll <- numeric(d$S)
  for (s in seq_len(d$S)) {
    a <- info[[s]]; ix <- a$ix; K <- a$K
    dl <- Xb[ix] + draws$H[ix, , drop = FALSE]
    zz <- resv(dl, cst[ix], sig_t, ginv)
    uu <- dl + mL[ix] + draws$E[ix, , drop = FALSE] * th$sig_xi
    ok <- rep(TRUE, nd)
    if (K == 0L) {
      for (j in seq_len(d$J)) ok <- ok & (u0 >= zz[j, ])
    } else {
      zS <- zz[a$S, , drop = FALSE]; uS <- uu[a$S, , drop = FALSE]; zK <- zS[K, ]
      if (K >= 2L) for (k in 1:(K - 1L)) ok <- ok & (zS[k, ] >= zS[k + 1L, ])
      if (length(a$U)) for (l in a$U) ok <- ok & (zK >= zz[l, ])
      if (K >= 2L) for (k in 1:(K - 1L)) ok <- ok & (uS[k, ] < zK)
      ok <- ok & (u0 < zK)
      znext <- if (length(a$U)) apply(zz[a$U, , drop = FALSE], 2, max) else -Inf
      ub <- if (a$ch == 0L) rep(u0, nd) else uS[a$ch_pos, ]
      if (length(a$U)) ok <- ok & (ub >= znext)
      for (k in seq_len(K)) if (k != a$ch_pos) ok <- ok & (ub >= uS[k, ])
      if (a$ch != 0L) ok <- ok & (ub >= u0)
    }
    ll[s] <- log(max(mean(ok), 1e-12))     # a zero cell must not be -Inf
  }
  -sum(ll)
}

# Two data sizes, and enough replications that the sd estimates are themselves
# usable. A single 6-replication sd has ~30% relative error, which is how an
# earlier draft of this project came to quote a variance ratio of 399x that did
# not survive re-measurement.
REPS <- 20
ND   <- 200
SIZES <- list(list(N = 25, T = 6, J = 5, seed = 31),
              list(N = 40, T = 8, J = 5, seed = 32))

run <- function(fn, label, d, info, p0) {
  v <- numeric(REPS); t0 <- Sys.time()
  for (r in seq_len(REPS))
    v[r] <- fn(p0, d, info, make_draws_sift(d, ND, seed = 700 + r), GINV)
  el <- as.numeric(Sys.time() - t0, units = "secs") / REPS
  # sd of an sd, for REPS normal draws, is about sd/sqrt(2(REPS-1))
  data.frame(simulator = label, mean_nll = mean(v), sd_nll = sd(v),
             sd_se = sd(v) / sqrt(2 * (REPS - 1)), sec_per_eval = el)
}

all <- list()
for (cfg in SIZES) {
  d    <- sim_sift(N = cfg$N, T = cfg$T, J = cfg$J, seed = cfg$seed, ginv = GINV)
  info <- prep_sift(d); p0 <- d$true_par
  r <- rbind(
    run(function(p, d, i, dr, g) nll_sift_cfs(p, d, i, dr, g), "crude frequency", d, info, p0),
    run(function(p, d, i, dr, g) nll_sift_ks (p, d, i, dr, g), "kernel-smoothed", d, info, p0),
    run(function(p, d, i, dr, g) nll_sift_ghk(p, d, i, dr, g), "GHK + analytic",  d, info, p0))
  r$tasks <- d$S
  r$var_ratio_vs_ghk <- round((r$sd_nll / r$sd_nll[3])^2, 1)
  # precision per second: matching GHK's sd needs var_ratio times as many draws,
  # and runtime is linear in draws, so this is the number that actually matters
  r$sec_for_ghk_precision <- round(r$sec_per_eval * r$var_ratio_vs_ghk, 3)
  all[[length(all) + 1L]] <- r
}
res <- do.call(rbind, all)
print(res, digits = 4)
save_result("simulators", list(res = res, nd = ND, reps = REPS,
  note = "log-likelihood at the true parameters; 20 independent draw sets per cell"))
