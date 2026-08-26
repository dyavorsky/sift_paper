# =============================================================================
# COMPARISON: Yavorsky, Honka & Chen (2021) vs SIFT
#
# Question. @Yavorsky_2021 shows the match-value standard deviation can be
# estimated rather than normalized to one, given an exogenous search-cost
# shifter, and that normalizing does real damage -- he estimates 8.16, and
# fixing it to one reverses a coefficient sign and understates a counterfactual
# by a factor of three to four.
#
# He has the shifter but no panel. @Morozov_2021 has the panel but normalizes.
# SIFT has both. Two things to establish:
#
#   (a) What does normalizing cost, as the true spread moves away from one?
#   (b) How much of the benefit comes from the PANEL as opposed to the shifter?
#       Holding total tasks fixed, T = 1 is the closest simulated analogue of
#       Yavorsky's design; T = 8 is SIFT.
#
# "Normalizing" is implemented the way a researcher actually does it: sigma_tilde
# is forced to 1 in the RESERVATION VALUE (the search rule), while kappa and
# sigma_xi still govern realized utility. See sig_fix in R/sift_model.R.
#
# The counterfactual is the managerially relevant quantity and mirrors
# Yavorsky's at-home-test-drive exercise: halve search costs, predict the change
# in clicks.
# =============================================================================
source("comparisons/_common.R")

TASKS <- 400L      # held fixed across the T arms so the panel is the only change
ND    <- 50L
NQ    <- 12L
BETA  <- c(1.0, -0.8)
KAPPA <- c(0.9, -0.6)
SIGXI <- 0.7

# forward-simulate search behaviour under a parameter vector
predict_clicks <- function(par, d, ginv, sig_fix = NULL, nsim = 40L, seed = 9) {
  th <- unpack(par)
  sig_t <- if (is.null(sig_fix)) sqrt(sum(th$kappa^2) + th$sig_xi^2) else sig_fix
  set.seed(seed)
  Xb <- as.vector(d$X %*% th$beta); mL <- as.vector(d$L %*% th$kappa)
  cst <- exp(th$gamma0 + th$gpos * d$pos)
  tot <- 0
  for (r in seq_len(nsim)) {
    delta <- Xb + rnorm(length(Xb))
    u     <- delta + mL + rnorm(length(Xb)) * th$sig_xi
    z     <- resv(delta, cst, sig_t, ginv)
    K <- 0
    for (s in seq_len(d$S)) {
      ix <- ((s - 1L) * d$J + 1L):(s * d$J)
      ord <- order(z[ix], decreasing = TRUE); best <- th$a0
      for (k in seq_len(d$J)) {
        if (z[ix][ord[k]] <= best) break
        K <- K + 1; best <- max(best, u[ix][ord[k]])
      }
    }
    tot <- tot + K / d$S
  }
  tot / nsim
}

rows <- list()
for (mult in c(1, 2, 4)) {
  kp <- KAPPA * mult; sx <- SIGXI * mult
  sig_true <- sqrt(sum(kp^2) + sx^2)
  for (T in c(1L, 8L)) {
    N <- as.integer(TASKS / T)
    d <- sim_sift(N = N, T = T, J = 5, beta = BETA, kappa = kp, sig_xi = sx,
                  seed = 3000L + 10L * mult + T, ginv = GINV)
    info <- prep_sift(d); dr <- make_draws_sift(d, nd = ND, seed = 77L)
    p0 <- rep(0, length(d$true_par)); p0[5] <- log(sx)

    fit <- function(sf) optim(p0, function(p)
      nll_sift_ghk(p, d, info, dr, GINV, nq = NQ, sig_fix = sf),
      method = "BFGS", control = list(maxit = 300, reltol = 1e-6))

    t0 <- Sys.time(); f_free <- fit(NULL); f_norm <- fit(1)
    el <- as.numeric(Sys.time() - t0, units = "mins")

    # counterfactual: halve search costs, predict mean clicks
    halve <- function(p) { p[7] <- p[7] - log(2); p }
    cf_true <- predict_clicks(halve(d$true_par), d, GINV)
    cf_free <- predict_clicks(halve(f_free$par), d, GINV)
    cf_norm <- predict_clicks(halve(f_norm$par), d, GINV, sig_fix = 1)
    base_true <- predict_clicks(d$true_par, d, GINV)

    g <- function(f) f$par[7]                          # gamma0
    gp <- function(f) f$par[8]                         # gamma_pos
    sig_hat <- sqrt(sum(f_free$par[3:4]^2) + exp(f_free$par[5])^2)

    rows[[length(rows) + 1L]] <- data.frame(
      sig_true = sig_true, T = T, N = N,
      sig_hat_free = sig_hat,
      gamma0_true = log(0.35), gamma0_free = g(f_free), gamma0_norm = g(f_norm),
      gpos_true = 0.12, gpos_free = gp(f_free), gpos_norm = gp(f_norm),
      clicks_base = base_true,
      cf_true = cf_true, cf_free = cf_free, cf_norm = cf_norm,
      minutes = el)
    cat(sprintf("sig_true=%.2f T=%d  sig_hat=%.2f  gamma0 %.2f/%.2f/%.2f  cf %.2f/%.2f/%.2f  (%.1f min)\n",
                sig_true, T, sig_hat, log(0.35), g(f_free), g(f_norm),
                cf_true, cf_free, cf_norm, el))
  }
}
res <- do.call(rbind, rows)
# counterfactual error, as a share of the true predicted change
res$cf_err_free <- (res$cf_free - res$cf_true) / (res$cf_true - res$clicks_base)
res$cf_err_norm <- (res$cf_norm - res$cf_true) / (res$cf_true - res$clicks_base)
print(res, digits = 3, row.names = FALSE)
save_result("mvsd", list(res = res, tasks = TASKS, nd = ND,
  note = "homogeneous fits; TASKS held fixed across T so the panel is the only change"))
