# =============================================================================
# COMPARISON: Ursu (2018) position analysis vs SIFT
#
# Question. Every position-based search-cost shifter rests on the assumption
# that position affects the COST of inspecting an alternative and nothing else.
# @Ursu_2018 §4.3 enumerates the ways position could instead enter the model and
# rules them out with observational tests. Can SIFT run those tests natively,
# and do they actually discriminate?
#
# Design. Three data-generating processes, identical except for where position
# acts:
#   "cost"     position raises the search cost           (the maintained assumption)
#   "utility"  position lowers pre-search utility        (Ursu's Cases 1,4,5,7)
#   "variance" position raises the post-search spread    (Ursu's Case 2)
#
# Two tests, both computable from raw SIFT data with no model:
#   T1  Among CLICKED alternatives, does position predict purchase?
#       Position-in-cost implies NO: position governed whether you looked, not
#       how much you liked what you saw.
#   T2  Is the first click near the top of the page?
#       Position-in-variance implies the opposite -- a wider post-search spread
#       raises the reservation value, so bad positions would be opened FIRST.
#
# The point of the comparison: Ursu must condition on covariates because
# position is chosen by the platform. We randomize it within every task, so both
# tests are experimental contrasts.
# =============================================================================
source("comparisons/_common.R")

sim_pos <- function(where, N = 400, T = 8, J = 6, seed = 1,
                    beta = c(1.0, -0.8), kappa = c(0.9, -0.6),
                    sig_xi = 0.7, a0 = 0.4, gamma0 = log(0.35),
                    effect = 0.35) {
  set.seed(seed)
  n <- N * T * J
  X <- matrix(rnorm(n * 2), n, 2); L <- matrix(rnorm(n * 2), n, 2)
  pos <- as.vector(replicate(N * T, sample.int(J)))          # RANDOMIZED
  sig_t <- sqrt(sum(kappa^2) + sig_xi^2)

  delta <- as.vector(X %*% beta) + rnorm(n)
  eps   <- as.vector(L %*% kappa) + rnorm(n) * sig_xi
  cost  <- rep(exp(gamma0), n)
  sg    <- rep(sig_t, n)

  if (where == "cost")     cost  <- exp(gamma0 + effect * pos)
  if (where == "utility")  delta <- delta - effect * pos
  if (where == "variance") sg    <- sig_t * exp(0.18 * pos)

  z <- resv(delta, cost, sg, GINV)
  cons <- matrix(0L, N * T, J); choice <- integer(N * T); firstpos <- integer(N * T)
  for (s in seq_len(N * T)) {
    ix <- ((s - 1) * J + 1):(s * J)
    zi <- z[ix]; ui <- delta[ix] + eps[ix]
    ord <- order(zi, decreasing = TRUE); best <- a0; K <- 0L
    for (k in seq_len(J)) { if (zi[ord[k]] <= best) break
                            K <- k; best <- max(best, ui[ord[k]]) }
    if (K > 0L) { cons[s, seq_len(K)] <- ord[seq_len(K)]; firstpos[s] <- pos[ix][ord[1]] }
    cand <- c(a0, ui[ord[seq_len(K)]])
    choice[s] <- c(0L, ord[seq_len(K)])[which.max(cand)]
  }
  list(J = J, S = N * T, pos = pos, cons_set = cons, choice = choice,
       firstpos = firstpos)
}

# --- the two tests, computed from raw data only ----------------------------
tests <- function(d) {
  # T1: among clicked alternatives, does position predict purchase?
  rec <- list()
  for (s in seq_len(d$S)) {
    k <- d$cons_set[s, ]; k <- k[k > 0]
    if (!length(k)) next
    ix <- ((s - 1L) * d$J + 1L):(s * d$J)
    rec[[length(rec) + 1L]] <- data.frame(pos = d$pos[ix][k],
                                          bought = as.integer(k == d$choice[s]))
  }
  cl <- do.call(rbind, rec)
  m  <- glm(bought ~ pos, family = binomial(), data = cl)
  ci <- suppressMessages(confint(m))
  fp <- d$firstpos[d$firstpos > 0]
  list(t1_slope = unname(coef(m)["pos"]),
       t1_lo = unname(ci["pos", 1]), t1_hi = unname(ci["pos", 2]),
       t1_p  = unname(summary(m)$coefficients["pos", 4]),
       n_clicked = nrow(cl),
       t2_first_pos = mean(fp), t2_uniform = (d$J + 1) / 2,
       t2_p = t.test(fp, mu = (d$J + 1) / 2)$p.value,
       mean_clicks = mean(rowSums(d$cons_set != 0)))
}

arms <- c("cost", "utility", "variance")
res <- do.call(rbind, lapply(seq_along(arms), function(i) {
  d <- sim_pos(arms[i], seed = 900 + i)
  cbind(data.frame(dgp = arms[i]), as.data.frame(tests(d)))
}))
print(res, digits = 3)
save_result("position", list(res = res,
  note = "N=400 x T=8, J=6; position randomized within task in all three arms"))
