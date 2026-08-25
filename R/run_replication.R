# =============================================================================
# Monte Carlo replication of the PM column of Chung, Chintagunta & Misra (2025)
# Table 1, using the R port in chung_pm.R.
#
# The DGP and estimator match the authors' shipped `main_for_loop.m` exactly.
# NOTE their Table 1 additionally puts a random coefficient on the first X
# (log sigma_theta1 = -0.916); the shipped code is the homogeneous version and
# so is this. We therefore compare on the eight shared parameters only.
#
# nd = 100 matches the paper (their §5.2.1); the shipped driver uses 500.
# RNG streams differ between R and Matlab, so individual datasets cannot be
# matched -- the comparison is distributional.
#
# Two things learned the hard way and encoded here:
#   * A two-stage fit (coarse nd, then refine) is far cheaper than optimising
#     at full nd from zeros. Starting values do not change the target.
#   * reltol must not be tighter than the simulation noise floor, or BFGS
#     grinds to maxit chasing numerical wiggle. 1e-8 is plenty.
#   * Results are written per replication, so a partial run is still usable.
# =============================================================================

suppressPackageStartupMessages(library(parallel))

NREPS <- 100
N     <- 1000
ND    <- 100      # final draw count, matching the paper
ND0   <- 25       # coarse draw count for stage 1
NCORE <- 14       # physical cores; 18 on 14 cores doubled per-fit time
OUT   <- "R/out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PROG <- file.path(OUT, "progress.log")
PART <- file.path(OUT, "parts"); dir.create(PART, showWarnings = FALSE)

one_rep <- function(r) {
  f <- file.path(PART, sprintf("rep_%03d.rds", r))
  if (file.exists(f)) return(NULL)
  source("R/chung_pm.R", local = TRUE)
  ginv <- make_ginv()
  dat  <- sim_ssm(N = N, seed = r, ginv = ginv)
  info <- prep_data(dat)
  grp  <- prep_groups(dat, info)
  d0   <- make_draws(dat$N, dat$J, nd = ND0, seed = 20000L + r)
  d1   <- make_draws(dat$N, dat$J, nd = ND,  seed = 10000L + r)
  p0   <- rep(0, length(dat$true_par))

  fit <- function(mode) {
    t0 <- Sys.time()
    s1 <- optim(p0, function(p) nll_grouped(p, dat, info, d0, ginv, grp, pr_F_use = mode),
                method = "BFGS", control = list(maxit = 200, reltol = 1e-6))
    s2 <- optim(s1$par, function(p) nll_grouped(p, dat, info, d1, ginv, grp, pr_F_use = mode),
                method = "BFGS", control = list(maxit = 200, reltol = 1e-8))
    data.frame(rep = r, mode = mode, conv = s2$convergence, nll = s2$value,
               secs = as.numeric(Sys.time() - t0, units = "secs"),
               t(setNames(s2$par, dat$par_names)),
               K_mean = mean(rowSums(dat$cons_set != 0)),
               stringsAsFactors = FALSE)
  }
  out <- rbind(fit("zK"), fit("znext"))
  saveRDS(out, f)
  cat(sprintf("rep %d done (%.0f s)\n", r, sum(out$secs)), file = PROG, append = TRUE)
  NULL
}

cat(sprintf("Running %d replications on %d cores (N=%d, nd=%d->%d)...\n",
            NREPS, NCORE, N, ND0, ND))
t0 <- Sys.time()
cl <- makeCluster(NCORE)
clusterExport(cl, c("N", "ND", "ND0", "PROG", "PART"))
invisible(parLapplyLB(cl, seq_len(NREPS), one_rep))
stopCluster(cl)

res <- do.call(rbind, lapply(list.files(PART, full.names = TRUE), readRDS))
saveRDS(res, file.path(OUT, "replication_raw.rds"))
write.csv(res, file.path(OUT, "replication_raw.csv"), row.names = FALSE)
cat(sprintf("done in %.1f min; %d rows written\n",
            as.numeric(Sys.time() - t0, units = "mins"), nrow(res)))
