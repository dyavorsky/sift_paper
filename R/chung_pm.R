# =============================================================================
# R port of the sequential-search SMLE simulator of
#   Chung, Chintagunta & Misra (2025), QME 23:105-164.
#
# Ported from the authors' Matlab, `ssm_lik_for_loop.m` and `main_for_loop.m`
#   Dropbox/resources/Chung Chintagunta Misra - Matlab Code for Seq Search SMLE/
#
# Scope matches the shipped Matlab exactly: cross-sectional, homogeneous
# preferences, one homogeneous mean search cost, exponential search costs,
# conditional on at least one search (first search free, and it also reveals
# the outside option).
#
# Deliberate departures from the Matlab, all documented at the point of use:
#   1. g() uses an asymptotic expansion in the upper tail instead of `vpa`.
#   2. Probabilities are computed in log space with pnorm(log.p=) / expm1().
#   3. All randomness is pre-drawn once (it is parameter-free), instead of
#      re-seeding inside the likelihood.
#   4. Pr(F) defaults to using z_K, matching the paper and the authors' own
#      VECTORIZED code. The shipped LOOP code uses z_2 instead; reproduce that
#      with pr_F_use = "znext". See note at the bottom of nll_loop().
# =============================================================================


# ---- g: reservation-value offset -> standardised search cost ---------------
# g(eta) = phi(eta) - eta * (1 - Phi(eta)),  strictly decreasing, g: R -> (0, Inf)
#
# Port of z_to_c.m. The authors build their spline knots under Matlab's
# variable-precision arithmetic, evidently fearing cancellation in the upper
# tail. In R that is unnecessary: the two terms differ by a factor of roughly
# 1/eta^2, so the subtraction costs only ~1.8 digits at eta = 7.8, and
# pnorm(lower.tail = FALSE) is accurate to ~1e-15 relative deep into the tail.
#
# Checked against a 500-level Mills-ratio continued fraction (independent of
# pnorm): relative error stays at ~1e-14 out to eta = 12. A truncated
# asymptotic expansion was tried and is far worse (~1e-6 at eta = 7.8), being
# a divergent series; do not reintroduce it.
g_fn <- function(eta, ub = 7.8) {
  e <- pmin(eta, ub)
  dnorm(e) - e * pnorm(e, lower.tail = FALSE)
}


# ---- g inverse, as a cubic spline built once ------------------------------
# Matlab: std_z descending over [-7.8, 7.8] step 0.01 so that std_c ascends,
# then spline(std_c, std_z). R's splinefun(method = "fmm") is the same
# not-a-knot cubic Matlab's spline() uses.
#
# NOTE the argument is the STANDARDISED cost c / sigma. One spline therefore
# serves every value of sigma: eta = ginv(c / sigma), z = E + eta * sigma.
make_ginv <- function(min_eta = -7.8, max_eta = 7.8, step = 0.01) {
  eta <- rev(seq(min_eta, max_eta, by = step))   # descending
  cst <- g_fn(eta, ub = max_eta)                 # ascending
  stopifnot(!is.unsorted(cst, strictly = TRUE))
  f <- splinefun(cst, eta, method = "fmm")
  structure(f, c_range = range(cst), eta_range = c(min_eta, max_eta))
}


# ---- truncated exponential draws ------------------------------------------
# Ports of cpu_exprnd_trunc_below.m / cpu_exprnd_trunc_above.m, but taking
# pre-drawn uniforms rather than calling rand() internally.
#
# Below: memorylessness makes this a shift. This is the property that keeps
# the whole unsearched-alternative block linear, and is the real reason the
# authors use an exponential rather than a log-normal.
rexp_trunc_below <- function(mu, lb, u) {
  lb <- pmax(lb, 0)
  -mu * log(1 - u) + lb
}

rexp_trunc_above <- function(mu, ub, u) {
  -mu * log(1 - u * (-expm1(-ub / mu)))
}


# ---- log Pr(c > x) and log Pr(c < x) for c ~ Exp(mean = mu) ---------------
log_surv_exp <- function(x, mu) -x / mu                  # log Pr(c > x)
log_cdf_exp  <- function(x, mu) log(-expm1(-x / mu))     # log Pr(c < x)


# =============================================================================
# Data generating process -- port of main_for_loop.m lines 41-88
# =============================================================================
sim_ssm <- function(N = 1000, J = 5,
                    alpha0 = c(-0.5, -0.8, 0.5, 0.8),   # product 1 fixed at 0
                    alpha1 = c(1.5, -1),
                    og_E   = 2.5,
                    c_bar  = 0.5,
                    sig_u  = 1,
                    X1_mean = 1.5, X1_sd = 0.4,
                    X2_mean = -1,  X2_sd = 0.2,
                    seed = 1, ginv = NULL) {

  if (is.null(ginv)) ginv <- make_ginv()

  set.seed(seed * 1000L + 1L); X1 <- matrix(rnorm(N * J), N, J) * X1_sd + X1_mean
  set.seed(seed * 1000L + 2L); X2 <- matrix(rnorm(N * J), N, J) * X2_sd + X2_mean

  a0 <- c(0, alpha0)
  E_util <- sweep(X1 * alpha1[1] + X2 * alpha1[2], 2, a0, "+")
  E_util <- cbind(E_util, og_E)                       # column J+1 = outside option

  set.seed(seed * 1000L + 3L)
  util <- E_util + matrix(rnorm(N * (J + 1)), N, J + 1) * sig_u

  set.seed(seed * 1000L + 4L)
  s_cost <- matrix(rexp(N * J, rate = 1 / c_bar), N, J)

  z_mat <- E_util[, 1:J, drop = FALSE] + ginv(s_cost / sig_u) * sig_u

  # Weitzman: search in descending z; the first search is free and also
  # reveals the outside option; continue while z_(k+1) > max(u_0, u_(1..k)).
  # The difference z_(k+1) - cummax_u(k) is non-increasing, so "last k
  # satisfied" (the Matlab's test) equals "first violation minus one".
  cons_set   <- matrix(0L, N, J)
  choice_vec <- integer(N)

  for (i in seq_len(N)) {
    z_i <- z_mat[i, ]
    u_i <- util[i, ]
    ord <- order(z_i, decreasing = TRUE)

    max_u    <- pmax(u_i[ord], u_i[J + 1])
    cummax_u <- cummax(max_u)

    sat <- which(z_i[ord][-1] - cummax_u[-J] > 0)
    K   <- if (length(sat) == 0L) 1L else max(sat) + 1L

    cons_set[i, seq_len(K)] <- ord[seq_len(K)]
    choice_vec[i] <- which.max(c(u_i[ord[seq_len(K)]], u_i[J + 1]))
    choice_vec[i] <- c(ord[seq_len(K)], J + 1L)[choice_vec[i]]
  }

  list(N = N, J = J, X1 = X1, X2 = X2, sig_u = sig_u,
       cons_set = cons_set, choice_vec = choice_vec,
       true_par = c(alpha0, alpha1, og_E, log(c_bar)),
       par_names = c(paste0("alpha0_", 2:(J)), paste0("alpha1_", 1:2),
                     "og_E", "log_c_bar"))
}


# =============================================================================
# Per-consumer bookkeeping, done once
# =============================================================================
prep_data <- function(dat) {
  N <- dat$N; J <- dat$J
  info <- vector("list", N)
  for (i in seq_len(N)) {
    cons  <- dat$cons_set[i, ]
    cons  <- cons[cons != 0L]
    nview <- length(cons)
    consA <- c(cons, J + 1L)                       # searched, plus outside option
    ch    <- dat$choice_vec[i]
    # coerce throughout: J is stored as a double, so J + 1L promotes and the
    # indices silently become numeric
    info[[i]] <- list(
      cons       = as.integer(cons),
      nview      = nview,
      consA      = as.integer(consA),
      choice     = as.integer(ch),
      cons_other = as.integer(consA[consA != ch]),
      non_cons   = as.integer(setdiff(seq_len(J), consA)),
      last_is_ch = identical(as.integer(ch), as.integer(cons[nview]))
    )
  }
  info
}


# =============================================================================
# Pre-drawn randomness. All of it is parameter-free, so it is drawn once and
# held fixed across likelihood evaluations -- which is what the Matlab's
# repeated rng() calls accomplish, more cheaply and without re-seeding.
# =============================================================================
make_draws <- function(N, J, nd, seed = 99) {
  set.seed(seed)
  n_u <- 2L + 2L * max(J - 2L, 0L)                 # c_K for G and H, plus loops
  list(nd = nd,
       Z = matrix(rnorm(N * nd), N, nd),           # for u_star
       U = array(runif(N * nd * n_u), c(N, nd, n_u)))
}


# =============================================================================
# Negative log-likelihood -- port of ssm_lik_for_loop.m
# =============================================================================
nll_loop <- function(pars, dat, info, draws, ginv,
                     pr_F_use = c("zK", "znext"),
                     return_vec = FALSE) {

  pr_F_use <- match.arg(pr_F_use)
  N <- dat$N; J <- dat$J; sig_u <- dat$sig_u
  nd <- draws$nd; z_ub <- 7.8

  n_a0   <- J - 1L
  alpha0 <- c(0, pars[1:n_a0])
  alpha1 <- pars[(n_a0 + 1L):(n_a0 + 2L)]
  og_E   <- pars[n_a0 + 3L]
  c_bar  <- exp(pars[n_a0 + 4L])

  if (!is.finite(c_bar) || c_bar <= 0) return(1e10)

  E_util <- sweep(dat$X1 * alpha1[1] + dat$X2 * alpha1[2], 2, alpha0, "+")
  E_util <- cbind(E_util, og_E)

  ll <- numeric(N)

  for (i in seq_len(N)) {
    v      <- info[[i]]
    K      <- v$nview
    Ei     <- E_util[i, ]
    E_last <- Ei[v$cons[K]]                          # last SEARCHED product
    Eoth   <- Ei[v$cons_other]
    Enon   <- if (length(v$non_cons)) Ei[v$non_cons] else numeric(0)

    u_star <- draws$Z[i, ] * sig_u + Ei[v$choice]

    # Pr(A) = Pr(u* > u_j') for the other considered alternatives
    log_pr_A <- colSums(pnorm(outer(-Eoth, u_star, "+") / sig_u, log.p = TRUE))

    # Pr(B) = Pr(z_K > u*) ; Pr(D) = Pr(u* > z_K)
    # u* > z_K  <=>  eta_K < (u* - E_K)/sigma  <=>  c_K > g((u* - E_K)/sigma)
    c_K_up   <- g_fn((u_star - E_last) / sig_u, z_ub) * sig_u
    log_pr_D <- log_surv_exp(c_K_up, c_bar)
    log_pr_B <- log_cdf_exp(c_K_up, c_bar)

    # Pr(C) = Pr(u* > max_l z_l) over unsearched l -- one division apiece
    log_pr_C <- if (length(Enon)) {
      colSums(log_surv_exp(g_fn(outer(-Enon, u_star, "+") / sig_u, z_ub) * sig_u, c_bar))
    } else rep(0, nd)

    branch1 <- branch2 <- NULL

    if (K == 1L) {
      # ---- B branch: z_1 > u*
      branch1 <- log_pr_A + log_pr_B + log_pr_C
      # ---- D branch: u* > z_1, so the stopping rule binds at z_1 (Pr E)
      c_K  <- rexp_trunc_below(c_bar, c_K_up, draws$U[i, , 1])
      z_K  <- E_last + pmin(z_ub, ginv(c_K / sig_u)) * sig_u
      log_pr_E <- if (length(Enon)) {
        colSums(log_surv_exp(g_fn(outer(-Enon, z_K, "+") / sig_u, z_ub) * sig_u, c_bar))
      } else rep(0, nd)
      branch2 <- log_pr_A + log_pr_D + log_pr_E

    } else {
      # ---- G recursion: z_1 > z_2 > ... > z_K, given z_K > u*
      c_K   <- rexp_trunc_above(c_bar, c_K_up, draws$U[i, , 1])
      z_G   <- E_last + pmin(z_ub, ginv(c_K / sig_u)) * sig_u
      lp_G  <- rep(0, nd)
      for (k in seq_len(K - 1L)) {
        E_k    <- Ei[v$cons[K - k]]
        c_up_G <- g_fn((z_G - E_k) / sig_u, z_ub) * sig_u
        lp_G   <- lp_G + log_cdf_exp(c_up_G, c_bar)
        if (k < K - 1L) {
          uG  <- draws$U[i, , 2L + k]
          z_G <- ginv(rexp_trunc_above(c_bar, c_up_G, uG) / sig_u) * sig_u + E_k
        }
      }
      branch1 <- log_pr_A + log_pr_B + log_pr_C + lp_G

      if (v$last_is_ch) {
        # ---- H recursion: same ordering, but on the u* > z_K side
        c_K_H  <- rexp_trunc_below(c_bar, c_K_up, draws$U[i, , 2])
        z_K_H  <- E_last + pmin(z_ub, ginv(c_K_H / sig_u)) * sig_u
        z_H    <- z_K_H
        lp_H   <- rep(0, nd)
        for (k in seq_len(K - 1L)) {
          E_k    <- Ei[v$cons[K - k]]
          c_up_H <- g_fn((z_H - E_k) / sig_u, z_ub) * sig_u
          lp_H   <- lp_H + log_cdf_exp(c_up_H, c_bar)
          if (k < K - 1L) {
            uH  <- draws$U[i, , 2L + (J - 2L) + k]
            z_H <- ginv(rexp_trunc_above(c_bar, c_up_H, uH) / sig_u) * sig_u + E_k
          }
        }

        # Pr(E): stopping rule evaluated at z_K (not u*), on this branch
        log_pr_E <- if (length(Enon)) {
          colSums(log_surv_exp(g_fn(outer(-Enon, z_K_H, "+") / sig_u, z_ub) * sig_u, c_bar))
        } else rep(0, nd)

        # Pr(F): every other searched utility must sit below z_K.
        # The binding constraint is the LAST continuation condition, and z is
        # decreasing along the search order, so z_K is the tightest bound.
        # The shipped LOOP code uses the post-recursion z (= z_2 when K >= 3),
        # a looser bound; the authors' VECTORIZED code and the paper both use
        # z_K. We default to z_K and keep the other for comparison.
        z_for_F  <- if (pr_F_use == "zK") z_K_H else z_H
        log_pr_F <- colSums(pnorm(outer(-Eoth, z_for_F, "+") / sig_u, log.p = TRUE))

        branch2 <- log_pr_D + log_pr_E + log_pr_F + lp_H
      }
    }

    p <- exp(branch1)
    if (!is.null(branch2)) p <- p + exp(branch2)
    ll[i] <- log(mean(p))
  }

  if (return_vec) return(ll)
  val <- -sum(ll)
  if (!is.finite(val)) return(1e10)
  val
}


# =============================================================================
# Grouped/vectorised likelihood.
#
# Same model, same draws, same arithmetic as nll_loop() -- but consumers are
# batched by (K, whether the purchase was the last search). Every consumer in
# a batch follows an identical code path with identical loop lengths, so the
# per-consumer R loop collapses to one pass per batch (at most 2J of them).
#
# nll_loop() is the readable reference and the test oracle; this is what gets
# called in anger. test_equivalence.R checks they agree to machine precision.
# =============================================================================

# splinefun() drops dim(); restore it so ginv() can be applied to matrices
ginv_m <- function(ginv, x) {
  d <- dim(x); r <- ginv(as.vector(x)); if (!is.null(d)) dim(r) <- d; r
}

prep_groups <- function(dat, info) {
  N <- dat$N; J <- dat$J
  key <- vapply(info, function(v) sprintf("%d_%d", v$nview, as.integer(v$last_is_ch)),
                character(1))
  # NB: do NOT use t(vapply(...)). When the inner length is 1, vapply returns a
  # bare vector and t() makes it 1 x ng instead of ng x 1 -- which silently
  # corrupts K = 1 (cons, cons_other) and K = J-1 (non_cons), and blows up the
  # cost of gather() at the same time. rbind is unambiguous.
  stack <- function(lst, m, ng) {
    if (m == 0L) return(matrix(0L, ng, 0))
    matrix(as.integer(unlist(lst, use.names = FALSE)), nrow = ng, ncol = m, byrow = TRUE)
  }
  lapply(split(seq_len(N), key), function(idx) {
    v1 <- info[[idx[1]]]; K <- v1$nview; ng <- length(idx)
    list(idx        = idx,
         K          = K,
         last_is_ch = v1$last_is_ch,
         choice     = vapply(info[idx], function(v) v$choice, integer(1)),
         cons       = stack(lapply(info[idx], `[[`, "cons"),       K,     ng),
         cons_other = stack(lapply(info[idx], `[[`, "cons_other"), K,     ng),
         non_cons   = stack(lapply(info[idx], `[[`, "non_cons"),   J - K, ng))
  })
}

nll_grouped <- function(pars, dat, info, draws, ginv, groups,
                        pr_F_use = c("zK", "znext"), return_vec = FALSE) {

  pr_F_use <- match.arg(pr_F_use)
  N <- dat$N; J <- dat$J; sig_u <- dat$sig_u
  nd <- draws$nd; z_ub <- 7.8

  n_a0   <- J - 1L
  alpha0 <- c(0, pars[1:n_a0])
  alpha1 <- pars[(n_a0 + 1L):(n_a0 + 2L)]
  og_E   <- pars[n_a0 + 3L]
  c_bar  <- exp(pars[n_a0 + 4L])
  if (!is.finite(c_bar) || c_bar <= 0) return(1e10)

  E_util <- sweep(dat$X1 * alpha1[1] + dat$X2 * alpha1[2], 2, alpha0, "+")
  E_util <- cbind(E_util, og_E)

  ll <- numeric(N)
  Ug <- function(idx, s) matrix(draws$U[idx, , s], nrow = length(idx))

  # E_util[cbind(row, col)] gathered for a matrix of column indices
  gather <- function(idx, cols) {
    matrix(E_util[cbind(rep(idx, ncol(cols)), as.vector(cols))],
           nrow = length(idx))
  }

  for (g in groups) {
    idx <- g$idx; ng <- length(idx); K <- g$K
    E_cons <- gather(idx, g$cons)                    # ng x K, in search order
    E_last <- E_cons[, K]
    E_oth  <- gather(idx, g$cons_other)              # ng x K
    E_non  <- if (ncol(g$non_cons)) gather(idx, g$non_cons) else matrix(0, ng, 0)

    u_star <- matrix(draws$Z[idx, ], nrow = ng) * sig_u + E_util[cbind(idx, g$choice)]

    lp_A <- matrix(0, ng, nd)
    for (m in seq_len(ncol(E_oth))) lp_A <- lp_A + pnorm((u_star - E_oth[, m]) / sig_u, log.p = TRUE)

    c_K_up <- g_fn((u_star - E_last) / sig_u, z_ub) * sig_u
    lp_D   <- log_surv_exp(c_K_up, c_bar)
    lp_B   <- log_cdf_exp(c_K_up, c_bar)

    lp_C <- matrix(0, ng, nd)
    for (m in seq_len(ncol(E_non)))
      lp_C <- lp_C + log_surv_exp(g_fn((u_star - E_non[, m]) / sig_u, z_ub) * sig_u, c_bar)

    if (K == 1L) {
      br1  <- lp_A + lp_B + lp_C
      c_K  <- rexp_trunc_below(c_bar, c_K_up, Ug(idx, 1))
      z_K  <- E_last + pmin(z_ub, ginv_m(ginv, c_K / sig_u)) * sig_u
      lp_E <- matrix(0, ng, nd)
      for (m in seq_len(ncol(E_non)))
        lp_E <- lp_E + log_surv_exp(g_fn((z_K - E_non[, m]) / sig_u, z_ub) * sig_u, c_bar)
      br2 <- lp_A + lp_D + lp_E
      p   <- exp(br1) + exp(br2)

    } else {
      c_K  <- rexp_trunc_above(c_bar, c_K_up, Ug(idx, 1))
      z_G  <- E_last + pmin(z_ub, ginv_m(ginv, c_K / sig_u)) * sig_u
      lp_G <- matrix(0, ng, nd)
      for (k in seq_len(K - 1L)) {
        E_k    <- E_cons[, K - k]
        c_up_G <- g_fn((z_G - E_k) / sig_u, z_ub) * sig_u
        lp_G   <- lp_G + log_cdf_exp(c_up_G, c_bar)
        if (k < K - 1L)
          z_G <- ginv_m(ginv, rexp_trunc_above(c_bar, c_up_G, Ug(idx, 2L + k)) / sig_u) * sig_u + E_k
      }
      p <- exp(lp_A + lp_B + lp_C + lp_G)

      if (g$last_is_ch) {
        c_K_H <- rexp_trunc_below(c_bar, c_K_up, Ug(idx, 2))
        z_K_H <- E_last + pmin(z_ub, ginv_m(ginv, c_K_H / sig_u)) * sig_u
        z_H   <- z_K_H
        lp_H  <- matrix(0, ng, nd)
        for (k in seq_len(K - 1L)) {
          E_k    <- E_cons[, K - k]
          c_up_H <- g_fn((z_H - E_k) / sig_u, z_ub) * sig_u
          lp_H   <- lp_H + log_cdf_exp(c_up_H, c_bar)
          if (k < K - 1L)
            z_H <- ginv_m(ginv, rexp_trunc_above(c_bar, c_up_H, Ug(idx, 2L + (J - 2L) + k)) / sig_u) * sig_u + E_k
        }
        lp_E <- matrix(0, ng, nd)
        for (m in seq_len(ncol(E_non)))
          lp_E <- lp_E + log_surv_exp(g_fn((z_K_H - E_non[, m]) / sig_u, z_ub) * sig_u, c_bar)

        z_for_F <- if (pr_F_use == "zK") z_K_H else z_H
        lp_F <- matrix(0, ng, nd)
        for (m in seq_len(ncol(E_oth)))
          lp_F <- lp_F + pnorm((z_for_F - E_oth[, m]) / sig_u, log.p = TRUE)

        p <- p + exp(lp_D + lp_E + lp_F + lp_H)
      }
    }
    ll[idx] <- log(rowMeans(p))
  }

  if (return_vec) return(ll)
  val <- -sum(ll)
  if (!is.finite(val)) return(1e10)
  val
}
