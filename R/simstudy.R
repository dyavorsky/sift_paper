# =============================================================================
# SIMULATION STUDY (section 6)
#
# Three questions, one design.
#
#  Q1  Does HB recover the population parameters -- both preferences and the
#      search-cost intercept -- at commercial task counts?
#  Q2  Does it recover INDIVIDUAL-level parameters? This is the claim that makes
#      SIFT a conjoint substitute, and it is the claim @Morozov_2021 is careful
#      about: his nonparametric identification result takes sessions per
#      consumer to infinity, and with a mean of 5 he reports hyperparameters
#      only. We should find out where the hierarchical prior stops doing most of
#      the work.
#  Q3  How many tasks does a practitioner need to field?
#
# Design: sweep T in {4, 8, 16} holding N fixed. Report population recovery,
# individual-level recovery (correlation and RMSE of posterior means against the
# truth), and the share of individual variance actually recovered.
#
# Run from the repo root:  Rscript R/simstudy.R
# =============================================================================
source("R/sift_hb.R")
ginv <- make_ginv()
dir.create("R/out", showWarnings = FALSE, recursive = TRUE)

N      <- 150L
J      <- 5L
TGRID  <- c(4L, 8L, 16L)
NITER  <- 5000L
BURN   <- 2500L
ND     <- 40L
NQ     <- 12L
NCORE  <- 12L

TH_BAR <- c(1.0, -0.8, 0.9, -0.6, 0.4, log(0.35))
SD_I   <- c(0.5, 0.4, 0.45, 0.3, 0.3, 0.4)
PN     <- c("beta1", "beta2", "kappa1", "kappa2", "a0", "gamma0")

out <- list()
for (T in TGRID) {
  cat(sprintf("\n===== T = %d  (%d tasks) =====\n", T, N * T))
  d    <- sim_sift_het(N = N, T = T, J = J, th_bar = TH_BAR, sd_i = SD_I,
                       seed = 100L + T, ginv = ginv)
  info <- prep_sift(d)
  dr   <- make_draws_sift(d, nd = ND, seed = 500L + T)

  fit <- sift_hb_par(d, info, dr, ginv, ncore = NCORE, niter = NITER,
                     burn = BURN, nq = NQ, verbose = 500L, seed = T)

  # population level
  pm  <- colMeans(fit$theta_bar)
  psd <- apply(fit$theta_bar, 2, sd)
  q   <- t(apply(fit$theta_bar, 2, quantile, c(0.025, 0.975)))
  pop <- data.frame(T = T, param = PN, truth = TH_BAR, post_mean = pm,
                    post_sd = psd, lo = q[, 1], hi = q[, 2],
                    covered = TH_BAR >= q[, 1] & TH_BAR <= q[, 2])

  # heterogeneity: posterior mean of the diagonal of Sigma vs the truth
  sd_hat <- sqrt(apply(fit$Sigma, 2, function(z) mean(z))[seq(1, 36, by = 7)])
  het <- data.frame(T = T, param = PN, sd_true = SD_I, sd_post = sd_hat)

  # individual level
  ind <- data.frame(T = T, param = PN,
    cor  = vapply(1:6, function(k) cor(fit$theta_i_mean[, k], d$th_i[, k]), 0),
    rmse = vapply(1:6, function(k) sqrt(mean((fit$theta_i_mean[, k] - d$th_i[, k])^2)), 0),
    sd_of_truth = apply(d$th_i, 2, sd))
  # share of individual variance recovered: 1 - MSE/Var, the shrinkage-aware R^2
  ind$var_explained <- 1 - (ind$rmse^2 / ind$sd_of_truth^2)

  cat("\npopulation:\n"); print(pop, digits = 3, row.names = FALSE)
  cat("\nindividual:\n"); print(ind, digits = 3, row.names = FALSE)

  out[[as.character(T)]] <- list(pop = pop, het = het, ind = ind,
                                 minutes = fit$minutes, logL = fit$logL,
                                 acc_common = fit$acc_common,
                                 common_post = colMeans(fit$common),
                                 tasks = d$S)
  saveRDS(out, "R/out/simstudy.rds")     # checkpoint after every T
}

cat("\n===== summary =====\n")
print(do.call(rbind, lapply(out, function(o)
  data.frame(T = o$pop$T[1], tasks = o$tasks, minutes = round(o$minutes, 1),
             pop_rmse = round(sqrt(mean((o$pop$post_mean - o$pop$truth)^2)), 3),
             pop_coverage = mean(o$pop$covered),
             ind_cor = round(mean(o$ind$cor), 3),
             ind_var_expl = round(mean(o$ind$var_explained), 3)))),
  row.names = FALSE)
