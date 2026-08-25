# Mini-tests for R/sift_model.R. Deliberately tiny -- these check that the
# machinery is wired correctly, not that the estimator is efficient.
source("R/sift_model.R")
ginv <- make_ginv()

cat("=========================================================\n")
cat("1. DGP produces sensible search behaviour\n")
cat("=========================================================\n")
d <- sim_sift(N = 60, T = 8, J = 6, seed = 3, ginv = ginv)
K <- rowSums(d$cons_set != 0)
cat(sprintf("  tasks                    : %d  (%d respondents x %d tasks)\n", d$S, d$N, d$T))
cat(sprintf("  sigma_tilde (variance-matched) : %.4f\n", d$sig_t))
cat("  clicks per task:\n"); print(round(prop.table(table(K)), 3))
cat(sprintf("  share no-purchase        : %.3f\n", mean(d$choice == 0)))
cat(sprintf("  share zero-click         : %.3f\n", mean(K == 0)))

# position was randomized -> it should predict click order but not preferences
first <- d$cons_set[, 1]; ok <- first > 0
pos1 <- d$pos[((which(ok) - 1) * d$J) + first[ok]]
cat(sprintf("  mean position of 1st click : %.3f  (uniform would be %.2f)\n",
            mean(pos1), (d$J + 1) / 2))

cat("\n=========================================================\n")
cat("2. Likelihood is finite, and a proper log-likelihood\n")
cat("=========================================================\n")
d2   <- sim_sift(N = 25, T = 6, J = 5, seed = 7, ginv = ginv)
info <- prep_sift(d2)
dr   <- make_draws_sift(d2, nd = 200, seed = 21)
p0   <- d2$true_par

t0 <- Sys.time(); nll_k <- nll_sift_ks (p0, d2, info, dr, ginv); tk <- Sys.time() - t0
t0 <- Sys.time(); nll_g <- nll_sift_ghk(p0, d2, info, dr, ginv); tg <- Sys.time() - t0
cat(sprintf("  kernel-smoothed  nll = %10.3f   (%.2f s)\n", nll_k, as.numeric(tk, units = "secs")))
cat(sprintf("  GHK + analytic   nll = %10.3f   (%.2f s)\n", nll_g, as.numeric(tg, units = "secs")))
cat(sprintf("  tasks = %d, so mean log-lik per task: KS %.3f, GHK %.3f\n",
            d2$S, -nll_k / d2$S, -nll_g / d2$S))
lv <- nll_sift_ghk(p0, d2, info, dr, ginv, return_vec = TRUE)
cat(sprintf("  all per-task log-liks finite and <= 0 : %s\n",
            all(is.finite(lv)) && all(lv <= 1e-8)))

cat("\n=========================================================\n")
cat("3. Simulation error: GHK should be far steadier than KS\n")
cat("=========================================================\n")
reps <- 6
vk <- vapply(seq_len(reps), function(r)
  nll_sift_ks (p0, d2, info, make_draws_sift(d2, 200, seed = 100 + r), ginv), 0)
vg <- vapply(seq_len(reps), function(r)
  nll_sift_ghk(p0, d2, info, make_draws_sift(d2, 200, seed = 100 + r), ginv), 0)
cat(sprintf("  KS  across %d draw sets : mean %9.3f   sd %.4f\n", reps, mean(vk), sd(vk)))
cat(sprintf("  GHK across %d draw sets : mean %9.3f   sd %.4f\n", reps, mean(vg), sd(vg)))
cat(sprintf("  variance ratio KS/GHK   : %.1fx\n", (sd(vk) / sd(vg))^2))

cat("\n=========================================================\n")
cat("4. Truth beats perturbations (both likelihoods)\n")
cat("=========================================================\n")
cat(sprintf("  %-14s %12s %12s\n", "parameter", "KS", "GHK"))
cat(sprintf("  %-14s %12.3f %12.3f   <- should be lowest\n", "at truth", nll_k, nll_g))
set.seed(5)
for (k in c(1, 3, 5, 8)) {
  pp <- p0; pp[k] <- pp[k] + 0.6
  cat(sprintf("  %-14s %12.3f %12.3f\n", paste0("+0.6 on ", d2$par_names[k]),
              nll_sift_ks(pp, d2, info, dr, ginv),
              nll_sift_ghk(pp, d2, info, dr, ginv)))
}

cat("\n=========================================================\n")
cat("5. GHK selection block integrates to a probability\n")
cat("=========================================================\n")
# Summing the likelihood over ALL possible (ordered search set, choice)
# outcomes for one task must give 1. Enumerate at J = 3.
J3 <- 3L
dd <- sim_sift(N = 1, T = 1, J = J3, seed = 2, ginv = ginv)
perms <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i)
    lapply(perms(v[-i]), function(q) c(v[i], q))))
outs <- list(list(S = integer(0), ch = 0L))
for (K in 1:J3) for (S in unique(lapply(perms(seq_len(J3)), function(q) q[seq_len(K)])))
  for (ch in c(0L, S)) outs[[length(outs) + 1L]] <- list(S = S, ch = ch)
tot <- 0
for (o in outs) {
  e <- dd; e$cons_set <- matrix(0L, 1, J3)
  if (length(o$S)) e$cons_set[1, seq_along(o$S)] <- o$S
  e$choice <- o$ch
  inf <- prep_sift(e)
  tot <- tot + exp(nll_sift_ghk(dd$true_par, e, inf,
                                make_draws_sift(e, 4000, seed = 99), ginv,
                                return_vec = TRUE))
}
cat(sprintf("  outcomes enumerated : %d\n", length(outs)))
cat(sprintf("  probabilities sum to: %.4f   (target 1.000)\n", tot))
