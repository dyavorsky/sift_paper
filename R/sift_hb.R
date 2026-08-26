# =============================================================================
# Hierarchical Bayes for SIFT.
#
# Why HB rather than random-coefficients SMLE. A mixed-logit-style SMLE has to
# integrate over the heterogeneity distribution, so every likelihood evaluation
# costs (draws for the search integral) x (mixing draws). HB never does that
# integral: each respondent sits at their own theta_i, so the mixing-draw
# multiplier disappears entirely. That is what makes HB the CHEAPER route to
# heterogeneity here, not the fallback -- see sections/identification.qmd.
#
# Model
#   theta_i = (beta_i, kappa_i, a0_i, gamma0_i) ~ N(theta_bar, Sigma)
#   common across respondents: log sigma_xi, gamma_pos
#
# theta_i carries BOTH preferences (beta, kappa, a0) and the search-cost
# intercept (gamma0), which is the individual-level deliverable that makes SIFT
# a conjoint substitute rather than a search paper.
#
# Sampler
#   1. Metropolis random walk on each theta_i, using ONLY that respondent's
#      tasks. This is the step that must not touch other respondents' rows;
#      ll_tasks_ghk() in sift_model.R exists for exactly this.
#   2. Conjugate normal draw for theta_bar.
#   3. Conjugate inverse-Wishart draw for Sigma.
#   4. Metropolis on the common parameters, over all tasks.
#
# Cost per iteration is therefore about two passes over the data, independent of
# how many parameters are heterogeneous.
# =============================================================================
source("R/sift_model.R")

# respondent -> the task rows belonging to them
resp_tasks <- function(d) split(seq_len(d$S), rep(seq_len(d$N), each = d$T))

# assemble the parameter list ll_tasks_ghk() expects
mk_th <- function(ti, common) list(
  beta = ti[1:2], kappa = ti[3:4], a0 = ti[5], gamma0 = ti[6],
  sig_xi = exp(common[1]), gpos = common[2])

sift_hb <- function(d, info, draws, ginv,
                    niter = 2000L, burn = 1000L, thin = 1L,
                    nq = 12L,
                    init = NULL, common_init = NULL,
                    step0 = 0.25, step_common = 0.05,
                    A0 = 100, nu0 = 9L, V0scale = 0.5,
                    verbose = 200L, progress_file = NULL) {
  P  <- 6L
  RT <- resp_tasks(d)
  N  <- d$N

  th_i   <- if (is.null(init)) matrix(0, N, P) else init
  common <- if (is.null(common_init)) c(log(0.7), 0.1) else common_init
  th_bar <- colMeans(th_i)
  Sig    <- diag(P) * 0.5
  V0     <- diag(P) * V0scale * nu0

  step  <- rep(step0, N)
  acc_i <- numeric(N); acc_c <- 0

  # current per-respondent log-likelihood
  ll_i <- vapply(seq_len(N), function(i)
    sum(ll_tasks_ghk(RT[[i]], d, info, draws, ginv, mk_th(th_i[i, ], common), nq)), 0)

  keep <- seq(burn + 1L, niter, by = thin)
  out_bar <- matrix(NA_real_, length(keep), P)
  out_sig <- array(NA_real_, c(length(keep), P, P))
  out_com <- matrix(NA_real_, length(keep), 2)
  out_ll  <- numeric(length(keep))
  th_sum  <- matrix(0, N, P); th_n <- 0L
  k <- 0L
  t0 <- Sys.time()

  for (it in seq_len(niter)) {

    ## 1. theta_i -----------------------------------------------------------
    Sinv <- solve(Sig)
    for (i in seq_len(N)) {
      prop <- th_i[i, ] + rnorm(P) * step[i]
      llp  <- sum(ll_tasks_ghk(RT[[i]], d, info, draws, ginv, mk_th(prop, common), nq))
      dpr  <- -0.5 * (t(prop - th_bar) %*% Sinv %*% (prop - th_bar) -
                      t(th_i[i, ] - th_bar) %*% Sinv %*% (th_i[i, ] - th_bar))
      if (log(runif(1)) < (llp - ll_i[i] + dpr)) {
        th_i[i, ] <- prop; ll_i[i] <- llp; acc_i[i] <- acc_i[i] + 1
      }
    }
    # adapt the random-walk scale during burn-in toward ~30% acceptance
    if (it <= burn && it %% 50L == 0L) {
      rate <- acc_i / 50
      step <- step * exp(pmax(pmin(rate - 0.30, 0.5), -0.5))
      step <- pmin(pmax(step, 0.01), 3)
      acc_i[] <- 0
    }

    ## 2. theta_bar | Sigma -------------------------------------------------
    Vp <- solve(N * Sinv + diag(P) / A0)
    mp <- Vp %*% (Sinv %*% colSums(th_i))
    th_bar <- as.vector(mp + t(chol(Vp)) %*% rnorm(P))

    ## 3. Sigma | theta -----------------------------------------------------
    Dv <- sweep(th_i, 2, th_bar)
    S  <- crossprod(Dv)
    Sig <- solve(drop(rWishart(1, nu0 + N, solve(V0 + S))))

    ## 4. common parameters -------------------------------------------------
    prop <- common + rnorm(2) * step_common
    llp  <- sum(vapply(seq_len(N), function(i)
      sum(ll_tasks_ghk(RT[[i]], d, info, draws, ginv, mk_th(th_i[i, ], prop), nq)), 0))
    if (log(runif(1)) < (llp - sum(ll_i))) {
      common <- prop; acc_c <- acc_c + 1
      ll_i <- vapply(seq_len(N), function(i)
        sum(ll_tasks_ghk(RT[[i]], d, info, draws, ginv, mk_th(th_i[i, ], common), nq)), 0)
    }

    if (it %in% keep) {
      k <- k + 1L
      out_bar[k, ] <- th_bar; out_sig[k, , ] <- Sig
      out_com[k, ] <- common; out_ll[k] <- sum(ll_i)
      th_sum <- th_sum + th_i; th_n <- th_n + 1L
    }
    if (verbose && it %% verbose == 0L) {
      msg <- sprintf("iter %5d/%d  logL %10.1f  acc_c %.2f  step %.2f  %.1f min\n",
                     it, niter, sum(ll_i), acc_c / it, mean(step),
                     as.numeric(Sys.time() - t0, units = "mins"))
      cat(msg)
      if (!is.null(progress_file)) cat(msg, file = progress_file, append = TRUE)
    }
  }

  list(theta_bar = out_bar, Sigma = out_sig, common = out_com, logL = out_ll,
       theta_i_mean = th_sum / max(th_n, 1L),
       acc_common = acc_c / niter, step = step,
       minutes = as.numeric(Sys.time() - t0, units = "mins"),
       niter = niter, burn = burn, nq = nq, nd = ncol(draws$H))
}

# ---------------------------------------------------------------------------
# Heterogeneous DGP: theta_i drawn from N(theta_bar, Sigma)
# ---------------------------------------------------------------------------
sim_sift_het <- function(N = 100, T = 8, J = 6,
                         th_bar = c(1.0, -0.8, 0.9, -0.6, 0.4, log(0.35)),
                         sd_i   = c(0.5, 0.4, 0.45, 0.3, 0.3, 0.4),
                         sig_xi = 0.7, gpos = 0.12, seed = 1, ginv = NULL) {
  if (is.null(ginv)) ginv <- make_ginv()
  set.seed(seed)
  P <- length(th_bar)
  th_i <- matrix(rnorm(N * P), N, P) %*% diag(sd_i) +
          matrix(th_bar, N, P, byrow = TRUE)
  n <- N * T * J
  X <- matrix(rnorm(n * 2), n, 2); L <- matrix(rnorm(n * 2), n, 2)
  pos <- as.vector(replicate(N * T, sample.int(J)))
  cons <- matrix(0L, N * T, J); choice <- integer(N * T)

  for (i in seq_len(N)) {
    b <- th_i[i, 1:2]; kp <- th_i[i, 3:4]; a0 <- th_i[i, 5]; g0 <- th_i[i, 6]
    sig_t <- sqrt(sum(kp^2) + sig_xi^2)
    for (tt in seq_len(T)) {
      s  <- (i - 1L) * T + tt
      ix <- ((s - 1L) * J + 1L):(s * J)
      delta <- as.vector(X[ix, ] %*% b) + rnorm(J)
      eps   <- as.vector(L[ix, ] %*% kp) + rnorm(J) * sig_xi
      z <- resv(delta, exp(g0 + gpos * pos[ix]), sig_t, ginv)
      ord <- order(z, decreasing = TRUE); best <- a0; K <- 0L
      u <- delta + eps
      for (k in seq_len(J)) { if (z[ord[k]] <= best) break
                              K <- k; best <- max(best, u[ord[k]]) }
      if (K > 0L) cons[s, seq_len(K)] <- ord[seq_len(K)]
      cand <- c(a0, u[ord[seq_len(K)]])
      choice[s] <- c(0L, ord[seq_len(K)])[which.max(cand)]
    }
  }
  list(N = N, T = T, J = J, S = N * T, X = X, L = L, pos = pos,
       cons_set = cons, choice = choice, th_i = th_i,
       th_bar = th_bar, sd_i = sd_i, sig_xi = sig_xi, gpos = gpos,
       par_names = c("beta1", "beta2", "kappa1", "kappa2", "a0", "gamma0"))
}

# =============================================================================
# Parallel HB.
#
# The theta_i step is embarrassingly parallel: each respondent's Metropolis
# update touches only their own tasks. Measured single-core cost is ~1.8 ms per
# task per iteration, which puts a commercial-scale fit (800 respondents x 10
# tasks, 20k iterations) at roughly 80 hours. Parallelising this step is what
# brings it inside an overnight budget, so the runtime claim in the paper
# depends on this function rather than on sift_hb().
#
# Workers hold the data and rebuild the spline themselves; each iteration ships
# only theta_bar, Sigma, the common parameters, and the current theta_i block.
# =============================================================================
hb_block <- function(idx, th_blk, ll_blk, step_blk, th_bar, Sinv, common, nq) {
  P <- ncol(th_blk); acc <- numeric(length(idx))
  for (n in seq_along(idx)) {
    i <- idx[n]
    prop <- th_blk[n, ] + rnorm(P) * step_blk[n]
    llp  <- sum(ll_tasks_ghk(.RT[[i]], .d, .info, .draws, .ginv,
                             mk_th(prop, common), nq))
    dpr  <- -0.5 * (t(prop - th_bar) %*% Sinv %*% (prop - th_bar) -
                    t(th_blk[n, ] - th_bar) %*% Sinv %*% (th_blk[n, ] - th_bar))
    if (log(runif(1)) < (llp - ll_blk[n] + dpr)) {
      th_blk[n, ] <- prop; ll_blk[n] <- llp; acc[n] <- 1
    }
  }
  list(th = th_blk, ll = ll_blk, acc = acc)
}

hb_ll_block <- function(idx, th_blk, common, nq)
  vapply(seq_along(idx), function(n)
    sum(ll_tasks_ghk(.RT[[idx[n]]], .d, .info, .draws, .ginv,
                     mk_th(th_blk[n, ], common), nq)), 0)

sift_hb_par <- function(d, info, draws, ginv, ncore = 6L,
                        niter = 2000L, burn = 1000L, thin = 1L, nq = 12L,
                        init = NULL, common_init = NULL,
                        step0 = 0.25, step_common = 0.05,
                        A0 = 100, nu0 = 9L, V0scale = 0.5,
                        verbose = 200L, seed = 1L) {
  suppressPackageStartupMessages(library(parallel))
  P <- 6L; N <- d$N; RT <- resp_tasks(d)
  blocks <- split(seq_len(N), cut(seq_len(N), ncore, labels = FALSE))

  cl <- makeCluster(ncore)
  on.exit(stopCluster(cl), add = TRUE)
  clusterSetRNGStream(cl, seed)
  clusterExport(cl, c("d", "info", "draws", "RT"), envir = environment())
  clusterEvalQ(cl, {
    source("R/sift_model.R")
    .d <- d; .info <- info; .draws <- draws; .RT <- RT
    .ginv <- make_ginv()
    mk_th <- function(ti, common) list(
      beta = ti[1:2], kappa = ti[3:4], a0 = ti[5], gamma0 = ti[6],
      sig_xi = exp(common[1]), gpos = common[2])
    NULL
  })
  clusterExport(cl, c("hb_block", "hb_ll_block"), envir = environment())

  th_i   <- if (is.null(init)) matrix(0, N, P) else init
  common <- if (is.null(common_init)) c(log(0.7), 0.1) else common_init
  th_bar <- colMeans(th_i); Sig <- diag(P) * 0.5
  V0 <- diag(P) * V0scale * nu0
  step <- rep(step0, N); acc_i <- numeric(N); acc_c <- 0

  ll_i <- unlist(clusterMap(cl, hb_ll_block, blocks,
                            lapply(blocks, function(b) th_i[b, , drop = FALSE]),
                            MoreArgs = list(common = common, nq = nq)))

  keep <- seq(burn + 1L, niter, by = thin)
  out_bar <- matrix(NA_real_, length(keep), P)
  out_sig <- array(NA_real_, c(length(keep), P, P))
  out_com <- matrix(NA_real_, length(keep), 2); out_ll <- numeric(length(keep))
  th_sum <- matrix(0, N, P); th_n <- 0L; k <- 0L
  t0 <- Sys.time()

  for (it in seq_len(niter)) {
    Sinv <- solve(Sig)
    r <- clusterMap(cl, hb_block, blocks,
                    lapply(blocks, function(b) th_i[b, , drop = FALSE]),
                    lapply(blocks, function(b) ll_i[b]),
                    lapply(blocks, function(b) step[b]),
                    MoreArgs = list(th_bar = th_bar, Sinv = Sinv,
                                    common = common, nq = nq))
    for (j in seq_along(blocks)) {
      b <- blocks[[j]]
      th_i[b, ] <- r[[j]]$th; ll_i[b] <- r[[j]]$ll; acc_i[b] <- acc_i[b] + r[[j]]$acc
    }
    if (it <= burn && it %% 50L == 0L) {
      step <- pmin(pmax(step * exp(pmax(pmin(acc_i / 50 - 0.30, 0.5), -0.5)), 0.01), 3)
      acc_i[] <- 0
    }

    Vp <- solve(N * Sinv + diag(P) / A0)
    th_bar <- as.vector(Vp %*% (Sinv %*% colSums(th_i)) + t(chol(Vp)) %*% rnorm(P))
    Sig <- solve(drop(rWishart(1, nu0 + N, solve(V0 + crossprod(sweep(th_i, 2, th_bar))))))

    prop <- common + rnorm(2) * step_common
    llp  <- unlist(clusterMap(cl, hb_ll_block, blocks,
                              lapply(blocks, function(b) th_i[b, , drop = FALSE]),
                              MoreArgs = list(common = prop, nq = nq)))
    if (log(runif(1)) < (sum(llp) - sum(ll_i))) {
      common <- prop; ll_i <- llp; acc_c <- acc_c + 1
    }

    if (it %in% keep) {
      k <- k + 1L
      out_bar[k, ] <- th_bar; out_sig[k, , ] <- Sig
      out_com[k, ] <- common; out_ll[k] <- sum(ll_i)
      th_sum <- th_sum + th_i; th_n <- th_n + 1L
    }
    if (verbose && it %% verbose == 0L)
      cat(sprintf("iter %5d/%d  logL %10.1f  acc_c %.2f  %.1f min\n",
                  it, niter, sum(ll_i), acc_c / it,
                  as.numeric(Sys.time() - t0, units = "mins")))
  }
  list(theta_bar = out_bar, Sigma = out_sig, common = out_com, logL = out_ll,
       theta_i_mean = th_sum / max(th_n, 1L), acc_common = acc_c / niter,
       step = step, minutes = as.numeric(Sys.time() - t0, units = "mins"),
       niter = niter, burn = burn, nq = nq, nd = ncol(draws$H), ncore = ncore)
}
