# =============================================================================
# COMPARISON: two instruments, one population
#
# Question. Does a conjoint SURVEY mismeasure preference, or does it measure
# preference correctly and then fail downstream when the estimates are carried
# into a market where consumers search?
#
# This distinction matters and is easy to blur. Our other conjoint comparison
# (exp_conjoint.R) fits a full-information choice model to SEARCH-GENERATED
# data, which is what one does with field clickstream. That is a
# misspecification result. It is NOT what a conjoint survey does: a conjoint
# survey constructs a full-information task, so the respondent genuinely does
# weigh the post-click attributes.
#
# Design. One population with fixed preferences (beta, kappa) and search costs.
# Field both instruments on it.
#
#   SIFT arm      respondents search, inspect, choose. Fit the search model.
#   Conjoint arm  the SAME attribute levels, all shown at once, no search.
#                 Respondent picks the best. Fit a conditional logit.
#
# Then ask two separate questions.
#
#   MEASUREMENT   does each instrument recover (beta, kappa)?
#   PREDICTION    does each reproduce market behaviour, where behaviour means
#                 the search-then-choose process the population actually follows?
#
# The prediction test is the sharp one, and position is the cleanest probe:
# conjoint has no position parameter at all, so it necessarily predicts that
# moving an alternative up the page does nothing.
# =============================================================================
source("comparisons/_common.R")

BETA  <- c(1.0, -0.8)
KAPPA <- c(0.9, -0.6)
SIGXI <- 0.7
A0    <- 0.4
G0    <- log(0.35)
GPOS  <- 0.12
N <- 150L; T <- 8L; J <- 6L

# ---------------------------------------------------------------------------
# Conjoint arm: same profiles, everything visible, no search.
# Both shocks are realized, because in a conjoint task the respondent knows
# their own tastes AND sees every attribute.
# ---------------------------------------------------------------------------
sim_conjoint_arm <- function(d, seed = 55) {
  set.seed(seed)
  n <- nrow(d$X)
  V <- as.vector(d$X %*% BETA) + as.vector(d$L %*% KAPPA)
  u <- V + rnorm(n) + rnorm(n) * SIGXI          # mu + eps, both realized
  ch <- integer(d$S)
  for (s in seq_len(d$S)) {
    ix <- ((s - 1L) * d$J + 1L):(s * d$J)
    cand <- c(A0 + rnorm(1), u[ix])
    ch[s] <- c(0L, seq_len(d$J))[which.max(cand)]
  }
  ch
}

# market behaviour: search-then-choose, given a parameter vector
simulate_market <- function(beta, kappa, sig_xi, a0, g0, gpos, d, nsim = 60L,
                            pos_override = NULL, seed = 7) {
  set.seed(seed)
  Xb <- as.vector(d$X %*% beta); mL <- as.vector(d$L %*% kappa)
  sig_t <- sqrt(sum(kappa^2) + sig_xi^2)
  pos <- if (is.null(pos_override)) d$pos else pos_override
  cst <- exp(g0 + gpos * pos)
  share_by_pos <- numeric(d$J); n_by_pos <- numeric(d$J)
  for (r in seq_len(nsim)) {
    delta <- Xb + rnorm(length(Xb))
    u <- delta + mL + rnorm(length(Xb)) * sig_xi
    z <- resv(delta, cst, sig_t, GINV)
    for (s in seq_len(d$S)) {
      ix <- ((s - 1L) * d$J + 1L):(s * d$J)
      ord <- order(z[ix], decreasing = TRUE); best <- a0; K <- 0L
      for (k in seq_len(d$J)) {
        if (z[ix][ord[k]] <= best) break
        K <- k; best <- max(best, u[ix][ord[k]])
      }
      cand <- c(a0, u[ix][ord[seq_len(K)]])
      pick <- c(0L, ord[seq_len(K)])[which.max(cand)]
      p <- pos[ix]
      for (j in seq_len(d$J)) n_by_pos[p[j]] <- n_by_pos[p[j]] + 1
      if (pick > 0L) share_by_pos[p[pick]] <- share_by_pos[p[pick]] + 1
    }
  }
  share_by_pos / n_by_pos
}

# conjoint's own prediction: full consideration, so position cannot enter
predict_conjoint <- function(par, d) {
  nb <- ncol(d$X); nk <- ncol(d$L)
  b <- par[1:nb]; k <- par[(nb + 1):(nb + nk)]; a0 <- par[nb + nk + 1]
  V <- as.vector(d$X %*% b) + as.vector(d$L %*% k)
  share_by_pos <- numeric(d$J); n_by_pos <- numeric(d$J)
  for (s in seq_len(d$S)) {
    ix <- ((s - 1L) * d$J + 1L):(s * d$J)
    v <- c(a0, V[ix]); p <- exp(v - max(v)); p <- p / sum(p)
    ps <- d$pos[ix]
    for (j in seq_len(d$J)) {
      share_by_pos[ps[j]] <- share_by_pos[ps[j]] + p[j + 1L]
      n_by_pos[ps[j]] <- n_by_pos[ps[j]] + 1
    }
  }
  share_by_pos / n_by_pos
}

# ---------------------------------------------------------------------------
cat("simulating the SIFT arm...\n")
d <- sim_sift(N = N, T = T, J = J, beta = BETA, kappa = KAPPA, sig_xi = SIGXI,
              a0 = A0, gamma0 = G0, gpos = GPOS, seed = 808L, ginv = GINV)
info <- prep_sift(d); dr <- make_draws_sift(d, nd = 100L, seed = 909L)
p0 <- rep(0, 8); p0[5] <- log(SIGXI)
t0 <- Sys.time()
f_sift <- optim(p0, function(p) nll_sift_ghk(p, d, info, dr, GINV, nq = 12L),
                method = "BFGS", control = list(maxit = 400, reltol = 1e-8))
cat(sprintf("  SIFT fit: conv=%d, %.1f min\n", f_sift$convergence,
            as.numeric(Sys.time() - t0, units = "mins")))

cat("simulating the conjoint arm...\n")
dc <- d; dc$choice <- sim_conjoint_arm(d)
dc$cons_set <- matrix(0L, d$S, d$J)               # no search in a conjoint task
f_conj <- fit_conjoint(dc)
cat(sprintf("  conjoint fit: conv=%d\n", f_conj$convergence))

# ---- MEASUREMENT: ratios, since logit and probit scales differ -------------
tr <- c(beta2_beta1 = BETA[2] / BETA[1],
        kappa1_beta1 = KAPPA[1] / BETA[1],
        kappa2_beta1 = KAPPA[2] / BETA[1])
rat <- function(b1, b2, k1, k2) c(b2 / b1, k1 / b1, k2 / b1)
r_sift <- rat(f_sift$par[1], f_sift$par[2], f_sift$par[3], f_sift$par[4])
r_conj <- rat(f_conj$par[1], f_conj$par[2], f_conj$par[3], f_conj$par[4])
meas <- data.frame(ratio = names(tr), truth = as.numeric(tr),
                   sift = r_sift, conjoint = r_conj,
                   sift_err = r_sift - tr, conjoint_err = r_conj - tr)
cat("\n--- MEASUREMENT (coefficient ratios) ---\n"); print(meas, digits = 3, row.names = FALSE)

# ---- PREDICTION: share by listing position --------------------------------
true_share <- simulate_market(BETA, KAPPA, SIGXI, A0, G0, GPOS, d)
sift_share <- simulate_market(f_sift$par[1:2], f_sift$par[3:4], exp(f_sift$par[5]),
                              f_sift$par[6], f_sift$par[7], f_sift$par[8], d)
conj_share <- predict_conjoint(f_conj$par, d)
pred <- data.frame(position = seq_len(J), truth = true_share,
                   sift = sift_share, conjoint = conj_share)
cat("\n--- PREDICTION (purchase share by listing position) ---\n")
print(pred, digits = 3, row.names = FALSE)
cat(sprintf("\n  true  share, position 1 vs %d : %.3f vs %.3f  (ratio %.2f)\n",
            J, true_share[1], true_share[J], true_share[1] / true_share[J]))
cat(sprintf("  SIFT  share, position 1 vs %d : %.3f vs %.3f  (ratio %.2f)\n",
            J, sift_share[1], sift_share[J], sift_share[1] / sift_share[J]))
cat(sprintf("  conj. share, position 1 vs %d : %.3f vs %.3f  (ratio %.2f)\n",
            J, conj_share[1], conj_share[J], conj_share[1] / conj_share[J]))

save_result("instruments", list(meas = meas, pred = pred, truth = tr,
  sift_par = f_sift$par, conj_par = f_conj$par, N = N, T = T, J = J,
  note = "same population, both instruments; measurement on ratios, prediction on share by position"))
