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

# GHK: P(A eta >= b) with eta ~ N(0, I_J), returning draws from the region.
#
# The weight is accumulated in LOG space. Accumulating it as a running product
# underflows to zero whenever the observed outcome is unlikely under the current
# parameters -- and a floored zero contributes ~ -691 to the log-likelihood,
# which swamps everything else and makes the optimiser flee the underflow region
# rather than fit the data. That single defect produced every spurious "bias"
# result in this project's diagnostics; see R/DIAGNOSTICS.md.
ghk <- function(A, b, U) {
  J <- nrow(A); nd <- ncol(U)
  Sg <- A %*% t(A)
  Ch <- tryCatch(t(chol(Sg)), error = function(e) NULL)
  if (is.null(Ch)) return(NULL)
  W <- matrix(0, J, nd); logw <- numeric(nd)
  for (r in seq_len(J)) {
    mu <- if (r == 1L) rep(0, nd) else as.vector(Ch[r, 1:(r - 1), drop = FALSE] %*%
                                                W[1:(r - 1), , drop = FALSE])
    lo  <- (b[r] - mu) / Ch[r, r]
    lpr <- pnorm(lo, lower.tail = FALSE, log.p = TRUE)      # log P(Z > lo)
    logw <- logw + lpr
    # draw Z | Z > lo entirely through the upper tail, so tiny probabilities
    # stay representable instead of collapsing to 0 or 1
    W[r, ] <- qnorm(lpr + log(pmax(U[r, ], 1e-300)),
                    lower.tail = FALSE, log.p = TRUE)
  }
  list(logw = logw, eta = solve(A, Ch %*% W))
}

# log(mean(exp(x))), overflow-safe
logmeanexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(mean(exp(x - m)))
}

# quadrature nodes on the probability scale for the one-dimensional eps integral
gq_nodes <- function(nq = 24L) (seq_len(nq) - 0.5) / nq

# Core: log-likelihood of ONE task. Kept separate from the loop so that the HB
# sampler can evaluate a single respondent's tasks at that respondent's own
# theta without touching anybody else's rows.
ll_task_ghk <- function(a, mu, m, mLs, u0, sig_t, sx, U, qn) {
  K <- a$K; J <- length(mu)
  sys <- sel_system(a, J)
  b <- if (K == 0L) -(u0 - m) else {
    bb <- numeric(J); r <- 0L
    if (K >= 2L) for (k in 1:(K - 1L)) { r <- r + 1L; bb[r] <- m[a$S[k + 1L]] - m[a$S[k]] }
    for (l in a$U) { r <- r + 1L; bb[r] <- m[l] - m[a$S[K]] }
    r <- r + 1L; bb[r] <- u0 - m[a$S[K]]
    bb
  }
  G <- ghk(sys$A, b, U)
  if (is.null(G)) return(-1e4)
  if (K == 0L) return(logmeanexp(G$logw))

  eta <- G$eta
  zK  <- m[a$S[K]] + eta[a$S[K], ]
  dS  <- mu[a$S] + eta[a$S, , drop = FALSE] + matrix(mLs[a$S], K, ncol(eta))
  znx <- if (length(a$U)) apply(m[a$U] + eta[a$U, , drop = FALSE], 2, max) else -Inf

  if (a$ch == 0L) {
    ok <- if (length(a$U)) as.numeric(u0 >= znx) else 1
    ub <- matrix(u0, K, ncol(eta))
    if (K >= 2L) ub[1:(K - 1L), ] <- pmin(ub[1:(K - 1L), , drop = FALSE],
                                          matrix(zK, K - 1L, ncol(eta), byrow = TRUE))
    lp <- colSums(pnorm((ub - dS) / sx, log.p = TRUE))
    return(logmeanexp(G$logw + log(ok) + lp))
  }
  cp <- a$ch_pos
  lo <- pmax(znx, u0) - dS[cp, ]
  hi <- if (cp < K) zK - dS[cp, ] else rep(Inf, ncol(eta))
  lo_p <- pnorm(lo / sx); hi_p <- pnorm(hi / sx)
  tail <- pmax(hi_p - lo_p, 0)
  acc  <- numeric(ncol(eta))
  others <- setdiff(seq_len(K), cp)
  keep <- others < K
  for (q in qn) {
    e  <- qnorm(pmin(pmax(lo_p + q * tail, 1e-300), 1 - 1e-16)) * sx
    uc <- dS[cp, ] + e
    if (length(others)) {
      ub <- matrix(uc, length(others), ncol(eta), byrow = TRUE)
      if (any(keep)) ub[keep, ] <- pmin(ub[keep, , drop = FALSE],
                                        matrix(zK, sum(keep), ncol(eta), byrow = TRUE))
      lp <- colSums(pnorm((ub - dS[others, , drop = FALSE]) / sx, log.p = TRUE))
    } else lp <- 0
    acc <- acc + exp(lp)
  }
  logmeanexp(G$logw + log(tail) + log(acc / length(qn)))
}

# Log-likelihood over a set of tasks. `beta`, `kappa`, `a0`, `gamma0` may vary
# by task (pass length-S vectors / S x p matrices) which is what the HB sampler
# needs; scalars are recycled.
ll_tasks_ghk <- function(tasks, d, info, draws, ginv, th, nq = 24L) {
  qn <- gq_nodes(nq)
  # th$sig_fix, when supplied, IMPOSES the post-search spread used in the
  # reservation value -- the literature's normalization -- while kappa and
  # sigma_xi continue to govern realized utility. That mismatch is exactly the
  # misspecification a normalizing analyst commits.
  sig_t <- if (!is.null(th$sig_fix)) th$sig_fix
           else sqrt(sum(th$kappa^2) + th$sig_xi^2)
  out <- numeric(length(tasks))
  for (n in seq_along(tasks)) {
    s <- tasks[n]; a <- info[[s]]; ix <- a$ix
    mu  <- as.vector(d$X[ix, , drop = FALSE] %*% th$beta)
    mLs <- as.vector(d$L[ix, , drop = FALSE] %*% th$kappa)
    cst <- exp(th$gamma0 + th$gpos * d$pos[ix])
    m   <- mu + ginv(pmin(pmax(cst / sig_t, attr(ginv, "c_range")[1]),
                          attr(ginv, "c_range")[2])) * sig_t
    out[n] <- ll_task_ghk(a, mu, m, mLs, th$a0, sig_t, th$sig_xi,
                          draws$U[[s]], qn)
  }
  out
}

nll_sift_ghk <- function(p, d, info, draws, ginv, nb = 2L, nk = 2L,
                         nq = 24L, return_vec = FALSE, tasks = NULL,
                         sig_fix = NULL) {
  th <- unpack(p, nb, nk); th$sig_fix <- sig_fix
  if (is.null(tasks)) tasks <- seq_len(d$S)
  ll <- pmax(ll_tasks_ghk(tasks, d, info, draws, ginv, th, nq), -1e4)
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
