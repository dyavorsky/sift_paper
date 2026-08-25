# =============================================================================
# SIFT: likelihood for the search-conjoint instrument.
#
# Builds on three sources:
#   * Yavorsky, Honka & Chen (2021, QME) -- utility decomposition, the
#     closed-form reservation utility of Kim et al. (2010), the kernel-smoothed
#     logit accept-reject simulator (their Appendix C), and the argument that an
#     *exogenous, product-specific* search-cost shifter identifies the
#     match-value standard deviation rather than forcing it to 1.
#   * Ursu, Seiler & Honka (2024) -- notation and the smoothing convention.
#   * GHK -- for the block of the likelihood where it legitimately applies.
#
# WHAT SIFT CHANGES relative to the field settings those papers study:
#   (1) A panel. Each respondent completes T tasks, so search spells are
#       observed repeatedly per person. This is what buys individual-level
#       heterogeneity rather than cross-sectional identification.
#   (2) Position on the landing page is RANDOMIZED, so the search-cost shifter
#       is exogenous by construction rather than by argument.
#   (3) The product page is DESIGNED. What a click reveals is L'kappa + xi,
#       where L is drawn by us from a known design distribution. From the
#       respondent's viewpoint the revealed component is therefore a mixture of
#       normals, not a single normal. We use the variance-matched normal
#       sigma_tilde^2 = kappa' Var(L) kappa + sigma_xi^2 in the reservation
#       value (see sections/model.qmd).
#   (4) A no-purchase outside option, and K = 0 (no clicks) is possible.
#
# NOTATION (task subscript t suppressed where unambiguous)
#   delta_ijt = X_ijt'beta + eta_ijt      known before clicking, eta ~ N(0,1)
#   u_ijt     = delta_ijt + L_ijt'kappa + xi_ijt   AFTER a click; xi ~ N(0,sigma_xi^2)
#   u_i0t     = a0                        outside option, known, never searched
#
#   The two views of the post-search component must be kept apart. BEFORE a
#   click the respondent does not know L, so from their viewpoint what a click
#   will reveal has the ex ante spread sigma_tilde -- that is what enters the
#   reservation value. AFTER a click, L is revealed to the respondent and is
#   observed by us (we designed it), so realized utility carries L'kappa in the
#   MEAN and only xi is random. Collapsing the two and using sigma_tilde for
#   realized utility leaves kappa entering solely through sigma_tilde, where it
#   trades off exactly against sigma_xi and is not identified.
#   c_ijt     = exp(gamma0 + gamma_pos * pos_ijt)      randomized position
#   z_ijt     = delta_ijt + zeta_ijt * sigma_tilde,  g(zeta) = c / sigma_tilde
# =============================================================================

source("R/chung_pm.R")     # g_fn(), make_ginv() -- shared, already validated

# ---------------------------------------------------------------------------
# Reservation values
# ---------------------------------------------------------------------------
resv <- function(delta, cost, sig, ginv) {
  r <- attr(ginv, "c_range")
  delta + ginv(pmin(pmax(cost / sig, r[1]), r[2])) * sig
}

# ---------------------------------------------------------------------------
# Data-generating process
# ---------------------------------------------------------------------------
sim_sift <- function(N = 200, T = 8, J = 6,
                     beta   = c(1.0, -0.8),      # landing-page attributes
                     kappa  = c(0.9, -0.6),      # product-page attributes
                     sig_xi = 0.7,
                     a0     = 0.4,               # outside-option utility
                     gamma0 = log(0.35),
                     gpos   = 0.12,              # position -> search cost
                     seed = 1, ginv = NULL) {
  if (is.null(ginv)) ginv <- make_ginv()
  set.seed(seed)
  n  <- N * T * J
  X  <- matrix(rnorm(n * length(beta)), n, length(beta))     # landing page
  L  <- matrix(rnorm(n * length(kappa)), n, length(kappa))   # product page
  id <- rep(seq_len(N), each = T * J)
  tk <- rep(rep(seq_len(T), each = J), times = N)
  pos <- as.vector(replicate(N * T, sample.int(J)))          # RANDOMIZED

  sig_t <- sqrt(sum(kappa^2) + sig_xi^2)     # variance-matched: Var(L) = I here
  delta <- as.vector(X %*% beta) + rnorm(n)
  eps   <- as.vector(L %*% kappa) + rnorm(n) * sig_xi
  cost  <- exp(gamma0 + gpos * pos)
  z     <- resv(delta, cost, sig_t, ginv)

  cons <- matrix(0L, N * T, J); choice <- integer(N * T)
  for (s in seq_len(N * T)) {
    ix <- ((s - 1) * J + 1):(s * J)
    zi <- z[ix]; ui <- delta[ix] + eps[ix]
    ord <- order(zi, decreasing = TRUE)
    best <- a0; K <- 0L
    for (k in seq_len(J)) {
      if (zi[ord[k]] <= best) break                 # stopping rule
      K <- k; best <- max(best, ui[ord[k]])         # search, update best
    }
    if (K > 0L) cons[s, seq_len(K)] <- ord[seq_len(K)]
    cand <- c(a0, ui[ord[seq_len(K)]])
    choice[s] <- c(0L, ord[seq_len(K)])[which.max(cand)]
  }
  list(N = N, T = T, J = J, S = N * T, X = X, L = L, pos = pos,
       id = id, task = tk, cons_set = cons, choice = choice,
       sig_t = sig_t,
       true_par = c(beta, kappa, log(sig_xi), a0, gamma0, gpos),
       par_names = c(paste0("beta", seq_along(beta)),
                     paste0("kappa", seq_along(kappa)),
                     "log_sig_xi", "a0", "gamma0", "gamma_pos"))
}

unpack <- function(p, nb = 2L, nk = 2L)
  list(beta = p[1:nb], kappa = p[(nb + 1):(nb + nk)],
       sig_xi = exp(p[nb + nk + 1]), a0 = p[nb + nk + 2],
       gamma0 = p[nb + nk + 3], gpos = p[nb + nk + 4])

# per-task bookkeeping, done once
prep_sift <- function(d) {
  K <- as.integer(rowSums(d$cons_set != 0))
  lapply(seq_len(d$S), function(s) {
    ix <- ((s - 1L) * d$J + 1L):(s * d$J)
    k  <- K[s]; S_ <- if (k > 0L) d$cons_set[s, seq_len(k)] else integer(0)
    list(ix = ix, K = k, S = S_, U = setdiff(seq_len(d$J), S_),
         ch = d$choice[s], ch_pos = if (d$choice[s] == 0L) 0L else match(d$choice[s], S_))
  })
}

# =============================================================================
# (A) Kernel-smoothed logit accept-reject simulator
#     Yavorsky et al. (2021) Appendix C / Ursu-Seiler-Honka smoothing.
#     Simple and general, but the smoothing induces bias that vanishes only as
#     the scale factors grow -- and large scale factors destroy the smoothness
#     the optimiser needs. That trade-off is the reason for (B).
# =============================================================================
nll_sift_ks <- function(p, d, info, draws, ginv,
                        lambda = c(15, 15, 15, 5), nb = 2L, nk = 2L,
                        return_vec = FALSE) {
  th <- unpack(p, nb, nk)
  sig_t <- sqrt(sum(th$kappa^2) + th$sig_xi^2)
  Xb  <- as.vector(d$X %*% th$beta)
  mL  <- as.vector(d$L %*% th$kappa)     # revealed on click, observed by us
  cst <- exp(th$gamma0 + th$gpos * d$pos)
  nd  <- ncol(draws$H)

  ll <- numeric(d$S)
  for (s in seq_len(d$S)) {
    a <- info[[s]]; ix <- a$ix; K <- a$K
    H <- draws$H[ix, , drop = FALSE]                 # J x nd  eta
    dl <- Xb[ix] + H                                 # delta
    zz <- resv(dl, cst[ix], sig_t, ginv)
    E  <- mL[ix] + draws$E[ix, , drop = FALSE] * th$sig_xi   # realized L, then xi
    uu <- dl + E
    u0 <- th$a0

    pen <- matrix(0, 1, nd)
    if (K == 0L) {                                   # no click: u0 >= all z
      for (j in seq_len(d$J)) pen <- pen + exp(-lambda[2] * (u0 - zz[j, ]))
    } else {
      zS <- zz[a$S, , drop = FALSE]; uS <- uu[a$S, , drop = FALSE]
      zK <- zS[K, ]
      # selection: z ordered along the observed search order
      if (K >= 2L) for (k in 1:(K - 1L))
        pen <- pen + exp(-lambda[2] * (zS[k, ] - zS[k + 1L, ]))
      if (length(a$U)) for (l in a$U)
        pen <- pen + exp(-lambda[2] * (zK - zz[l, ]))
      # continuation: binding constraint is z_K (z decreasing along the order)
      if (K >= 2L) for (k in 1:(K - 1L))
        pen <- pen + exp(-lambda[1] * (zK - uS[k, ]))
      pen <- pen + exp(-lambda[1] * (zK - u0))
      # stopping + choice
      znext <- if (length(a$U)) apply(zz[a$U, , drop = FALSE], 2, max) else -Inf
      ubest <- if (a$ch == 0L) rep(u0, nd) else uS[a$ch_pos, ]
      if (length(a$U)) pen <- pen + exp(-lambda[3] * (ubest - znext))
      for (k in seq_len(K)) if (k != a$ch_pos)
        pen <- pen + exp(-lambda[4] * (ubest - uS[k, ]))
      if (a$ch != 0L) pen <- pen + exp(-lambda[4] * (ubest - u0))
    }
    ll[s] <- log(mean(1 / (1 + pen)))
  }
  if (return_vec) return(ll)
  -sum(ll)
}

# =============================================================================
# (B) GHK on the eta block, analytic on the eps block.
#
# The structural fact that makes this work: CONDITIONAL ON THE OBSERVED CHOICE,
# every Weitzman condition is a conjunction of linear inequalities. The stopping
# rule max_h u_h >= max_l z_l is a union in general -- but the choice rule names
# which h attains the max, so it collapses to u_{choice} >= z_l for all l.
#
# Two further simplifications:
#   * Continuation requires u_h < z_k for every h < k <= K. Since z decreases
#     along the search order, the binding constraint is z_K alone.
#   * Conditional on eta AND on eps_{choice}, the remaining eps are mutually
#     independent one-sided constraints, so they integrate out in closed form as
#     a product of normal CDFs. Only a ONE-dimensional integral is left.
#
# So: GHK draws eta from the truncated selection region (smooth, low variance),
# and the eps block is evaluated by quadrature rather than simulation.
# =============================================================================

# constraint matrix A with A %*% eta >= b encoding selection + "search at all"
sel_system <- function(a, J) {
  K <- a$K
  if (K == 0L) {                       # u0 >= z_j for all j  ->  -eta_j >= ...
    return(list(A = -diag(J), kind = "none"))
  }
  A <- matrix(0, J, J); r <- 0L
  if (K >= 2L) for (k in 1:(K - 1L)) {
    r <- r + 1L; A[r, a$S[k]] <- 1; A[r, a$S[k + 1L]] <- -1
  }
  for (l in a$U) { r <- r + 1L; A[r, a$S[K]] <- 1; A[r, l] <- -1 }
  r <- r + 1L; A[r, a$S[K]] <- 1               # z_{s_K} > u0
  list(A = A, kind = "search")
}

# GHK: P(A eta >= b) with eta ~ N(0, I_J), returning draws from the region
ghk <- function(A, b, U) {
  J <- nrow(A); nd <- ncol(U)
  Sg <- A %*% t(A)
  Ch <- tryCatch(t(chol(Sg)), error = function(e) NULL)
  if (is.null(Ch)) return(NULL)
  W <- matrix(0, J, nd); w <- rep(1, nd)
  for (r in seq_len(J)) {
    mu <- if (r == 1L) rep(0, nd) else as.vector(Ch[r, 1:(r - 1), drop = FALSE] %*%
                                                W[1:(r - 1), , drop = FALSE])
    lo <- (b[r] - mu) / Ch[r, r]
    pr <- pnorm(lo, lower.tail = FALSE)
    pr <- pmin(pmax(pr, 1e-300), 1)
    w  <- w * pr
    # draw from N(0,1) truncated below at lo, by inverse CDF
    W[r, ] <- qnorm(pmin(pmax(pnorm(lo) + U[r, ] * pr, 1e-300), 1 - 1e-16))
  }
  # W holds the standardized draws v; the constrained vector is w = Ch %*% v,
  # and eta solves A eta = w.
  list(w = w, eta = solve(A, Ch %*% W))
}

# quadrature nodes on the probability scale for the one-dimensional eps integral
gq_nodes <- function(nq = 24L) (seq_len(nq) - 0.5) / nq

nll_sift_ghk <- function(p, d, info, draws, ginv, nb = 2L, nk = 2L,
                         nq = 24L, return_vec = FALSE) {
  th <- unpack(p, nb, nk)
  sig_t <- sqrt(sum(th$kappa^2) + th$sig_xi^2)
  Xb  <- as.vector(d$X %*% th$beta)
  mL  <- as.vector(d$L %*% th$kappa)     # revealed on click, observed by us
  cst <- exp(th$gamma0 + th$gpos * d$pos)
  qn  <- gq_nodes(nq)
  u0  <- th$a0
  ll  <- numeric(d$S)

  for (s in seq_len(d$S)) {
    a <- info[[s]]; ix <- a$ix; K <- a$K; J <- d$J
    m  <- Xb[ix] + ginv(cst[ix] / sig_t) * sig_t     # z = m + eta
    mu <- Xb[ix]                                     # delta = mu + eta
    sys <- sel_system(a, J)
    b <- if (K == 0L) -(u0 - m) else {
      bb <- numeric(J); r <- 0L
      if (K >= 2L) for (k in 1:(K - 1L)) { r <- r + 1L; bb[r] <- m[a$S[k + 1L]] - m[a$S[k]] }
      for (l in a$U) { r <- r + 1L; bb[r] <- m[l] - m[a$S[K]] }
      r <- r + 1L; bb[r] <- u0 - m[a$S[K]]
      bb
    }
    G <- ghk(sys$A, b, draws$U[[s]])
    if (is.null(G)) { ll[s] <- -1e6; next }

    if (K == 0L) { ll[s] <- log(mean(G$w)); next }   # no eps to integrate

    eta <- G$eta                                     # J x nd, inside the region
    zK  <- m[a$S[K]] + eta[a$S[K], ]
    # A click reveals the product page. Because WE designed it, the revealed
    # attributes L are observed for every searched alternative, so realized
    # utility is delta + L'kappa + xi with xi ~ N(0, sigma_xi^2). The ex ante
    # spread sigma_tilde belongs in the reservation value (computed BEFORE the
    # click); it must not be reused here. Folding L into the mean is what
    # identifies kappa -- without it kappa enters only through sigma_tilde and
    # trades off exactly against sigma_xi.
    dS  <- mu[a$S] + eta[a$S, , drop = FALSE] +
           matrix(mL[ix][a$S], K, ncol(eta))         # K x nd, delta + L'kappa
    sx  <- th$sig_xi
    znx <- if (length(a$U)) apply(m[a$U] + eta[a$U, , drop = FALSE], 2, max) else -Inf

    if (a$ch == 0L) {
      # no purchase: u0 beats every searched u, and u0 >= z_next
      ok <- if (length(a$U)) as.numeric(u0 >= znx) else 1
      ub <- matrix(u0, K, ncol(eta))
      if (K >= 2L) ub[1:(K - 1L), ] <- pmin(ub[1:(K - 1L), , drop = FALSE],
                                            matrix(zK, K - 1L, ncol(eta), byrow = TRUE))
      lp <- colSums(pnorm((ub - dS) / sx, log.p = TRUE))
      ll[s] <- log(mean(G$w * ok * exp(lp)))
    } else {
      cp <- a$ch_pos
      lo <- pmax(znx, u0) - dS[cp, ]                 # xi_choice >= lo
      # If the chosen alternative was NOT the last click, continuing past it
      # required u_choice < z_K. That caps xi_choice from above.
      hi <- if (cp < K) zK - dS[cp, ] else rep(Inf, ncol(eta))
      lo_p <- pnorm(lo / sx); hi_p <- pnorm(hi / sx)
      tail <- pmax(hi_p - lo_p, 0)
      acc <- numeric(ncol(eta))
      for (q in qn) {                                # integrate xi_choice
        e <- qnorm(pmin(pmax(lo_p + q * tail, 1e-300), 1 - 1e-16)) * sx
        uc <- dS[cp, ] + e
        others <- setdiff(seq_len(K), cp)
        if (length(others)) {
          ub <- matrix(uc, length(others), ncol(eta), byrow = TRUE)
          keep <- others < K
          if (any(keep)) ub[keep, ] <- pmin(ub[keep, , drop = FALSE],
                                            matrix(zK, sum(keep), ncol(eta), byrow = TRUE))
          lp <- colSums(pnorm((ub - dS[others, , drop = FALSE]) / sx, log.p = TRUE))
        } else lp <- 0
        acc <- acc + exp(lp)
      }
      ll[s] <- log(mean(G$w * tail * acc / nq))
    }
  }
  ll <- pmax(ll, -1e4)          # a zero-probability task must not return -Inf
  if (return_vec) return(ll)
  -sum(ll)
}

# ---------------------------------------------------------------------------
# Draws: all parameter-free randomness, made once
# ---------------------------------------------------------------------------
make_draws_sift <- function(d, nd, seed = 11) {
  set.seed(seed)
  list(H = matrix(rnorm(d$S * d$J * nd), d$S * d$J, nd),
       E = matrix(rnorm(d$S * d$J * nd), d$S * d$J, nd),
       U = lapply(seq_len(d$S), function(s) matrix(runif(d$J * nd), d$J, nd)))
}
