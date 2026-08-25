# =============================================================================
# The decisive correctness test for the likelihood.
#
# For one consumer with fixed pre-search utilities there are finitely many
# possible outcomes (ordered search set, purchase). The model assigns each a
# probability, and they must sum to one. Two checks:
#
#   A. Do the simulated likelihoods over all enumerated outcomes sum to 1?
#   B. Do they match brute-force frequencies from simulating the DGP directly?
#
# (B) is the cheap version of the authors' own SSE check (their §5.3), which
# they run against 100 million simulated consumers.
#
# IMPORTANT: the simulated likelihood is an unbiased but NOISY estimate of each
# probability, and at these settings its Monte Carlo error is several times the
# brute-force error. Deviations must therefore be scored against BOTH sources,
# which is why the model side is replicated below rather than computed once.
#
# This file deliberately re-implements the search loop rather than calling
# sim_ssm(), so the DGP used for the brute force is written independently of
# the one under test.
# =============================================================================
source("R/chung_pm.R")

J     <- 3L                # 48 outcomes at J = 3
SIG_U <- 1
C_BAR <- 0.5
ND    <- 4000L             # draws per model replicate
NREP  <- 12L               # independent model replicates
NSIM  <- 2e6               # brute-force consumers

ginv <- make_ginv()
E1   <- c(0.30, -0.20, 0.55, 0.40)      # J products + outside option
stopifnot(length(E1) == J + 1L)

# ---- enumerate every reachable (ordered search set, purchase) --------------
perms <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i)
    lapply(perms(v[-i]), function(p) c(v[i], p))))

outcomes <- list()
for (K in 1:J)
  for (S in unique(lapply(perms(seq_len(J)), function(p) p[seq_len(K)])))
    for (ch in c(S, J + 1L))
      outcomes[[length(outcomes) + 1L]] <- list(cons = S, choice = ch)
key      <- vapply(outcomes, function(o) paste0(paste(o$cons, collapse = "-"), "|", o$choice), "")
outcomes <- outcomes[!duplicated(key)]; key <- key[!duplicated(key)]
M <- length(outcomes)
cat(sprintf("J = %d  ->  %d enumerated outcomes\n\n", J, M))

# ---- pack as a synthetic dataset, all rows sharing the same E_util ---------
X1 <- matrix(0, M, J); for (j in 1:J) X1[, j] <- E1[j]   # alpha1 = c(1,0) => E = X1
cons_set <- matrix(0L, M, J)
for (m in 1:M) cons_set[m, seq_along(outcomes[[m]]$cons)] <- outcomes[[m]]$cons
dat <- list(N = M, J = J, X1 = X1, X2 = matrix(0, M, J), sig_u = SIG_U,
            cons_set = cons_set,
            choice_vec = vapply(outcomes, function(o) as.integer(o$choice), 1L),
            true_par = c(0, 0, 1, 0, E1[J + 1L], log(C_BAR)))
info <- prep_data(dat); grp <- prep_groups(dat, info)

reps <- vapply(seq_len(NREP), function(s) {
  dr <- make_draws(M, J, nd = ND, seed = 3000L + 17L * s)
  exp(nll_grouped(dat$true_par, dat, info, dr, ginv, grp, return_vec = TRUE))
}, numeric(M))
p_model  <- rowMeans(reps)
se_model <- apply(reps, 1, sd) / sqrt(NREP)

cat("=========================================================\n")
cat("A. Do the model probabilities sum to 1?\n")
cat("=========================================================\n")
cat(sprintf("  sum over %d outcomes  : %.6f\n", M, sum(p_model)))
cat(sprintf("  MC s.e. of the sum    : %.6f\n", sqrt(sum(se_model^2))))
cat(sprintf("  deviation in s.e.     : %.2f\n", abs(sum(p_model) - 1) / sqrt(sum(se_model^2))))
cat(sprintf("  any negative / NaN    : %s\n\n", any(!is.finite(p_model) | p_model < 0)))

# ---- brute force ----------------------------------------------------------
cat("=========================================================\n")
cat("B. Model probabilities vs brute-force DGP frequencies\n")
cat("=========================================================\n")
set.seed(4242)
U  <- matrix(rnorm(NSIM * (J + 1)), NSIM, J + 1) * SIG_U +
      matrix(E1, NSIM, J + 1, byrow = TRUE)
Z  <- matrix(E1[1:J], NSIM, J, byrow = TRUE) +
      ginv(matrix(rexp(NSIM * J, rate = 1 / C_BAR), NSIM, J) / SIG_U) * SIG_U

bkey <- character(NSIM)
for (i in seq_len(NSIM)) {
  ord <- order(Z[i, ], decreasing = TRUE)
  cm  <- cummax(pmax(U[i, ord], U[i, J + 1L]))
  sat <- which(Z[i, ord][-1] - cm[-J] > 0)
  K   <- if (length(sat) == 0L) 1L else max(sat) + 1L
  S   <- ord[seq_len(K)]
  bkey[i] <- paste0(paste(S, collapse = "-"), "|",
                    c(S, J + 1L)[which.max(c(U[i, S], U[i, J + 1L]))])
}
tab   <- table(bkey)
p_sim <- as.numeric(tab[key]); p_sim[is.na(p_sim)] <- 0; p_sim <- p_sim / NSIM
se_sim <- sqrt(pmax(p_sim, 1e-12) * (1 - p_sim) / NSIM)
se_tot <- sqrt(se_model^2 + se_sim^2)
z      <- (p_model - p_sim) / se_tot

ord <- order(-p_sim)
cat(sprintf("  %-10s %9s %9s %10s %9s %6s\n", "outcome", "model", "brute", "diff", "se(both)", "z"))
for (m in head(ord, 12))
  cat(sprintf("  %-10s %9.6f %9.6f %10.2e %9.2e %6.2f\n",
              key[m], p_model[m], p_sim[m], p_model[m] - p_sim[m], se_tot[m], z[m]))
cat(sprintf("\n  sum |model - brute|  : %.5f\n", sum(abs(p_model - p_sim))))
cat(sprintf("  max |z|              : %.2f   (expect ~2-3 across %d outcomes)\n", max(abs(z)), M))
cat(sprintf("  RMS z                : %.2f   (expect ~1 if both are correct)\n", sqrt(mean(z^2))))
cat(sprintf("  mean z               : %+.2f  (expect ~0; a systematic sign is a bug)\n", mean(z)))
