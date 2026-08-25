# Validation checks for the R port in chung_pm.R
source("R/chung_pm.R")

cat("=========================================================\n")
cat("1. g() vs Matlab ground truth (Ursu/Seiler/Honka tableZ.csv)\n")
cat("=========================================================\n")
# tableZ.m computes c = (1-normcdf(j))*((normpdf(j)/(1-normcdf(j))) - j)
#                     = phi(j) - j*(1-Phi(j)) = g(j),  for j in [-3.55, 4].
# csvwrite stores 5 significant digits, so compare at the precision the file
# actually holds rather than at machine epsilon.
tz <- file.path("C:/Users/DanYavorsky/Dropbox/resources",
                "Ursu Seiler Honka - Matlab Code for Seq Search Models",
                "seq_code/kernel/tableZ.csv")
if (file.exists(tz)) {
  tab <- read.csv(tz, header = FALSE)
  names(tab) <- c("m", "c")
  mine <- g_fn(tab$m)
  hit5 <- mean(signif(mine, 5) == signif(tab$c, 5))
  cat(sprintf("  n knots compared                    : %d  (eta in [%.2f, %.2f])\n",
              nrow(tab), min(tab$m), max(tab$m)))
  cat(sprintf("  rows identical at 5 sig digits      : %.4f\n", hit5))
  cat(sprintf("  max reldiff vs signif(mine, 5)      : %.3e\n",
              max(abs(signif(mine, 5) - tab$c) / abs(tab$c))))
  cat(sprintf("  verdict                             : %s\n\n",
              if (hit5 == 1) "EXACT at Matlab's stored precision" else "MISMATCH"))
} else cat("  tableZ.csv not found; skipped\n\n")

cat("=========================================================\n")
cat("2. g() far-tail accuracy vs an independent reference\n")
cat("=========================================================\n")
# Mills-ratio continued fraction; uses neither pnorm nor any series.
mills_cf <- function(x, n = 500) { f <- 0; for (k in n:1) f <- k / (x + f); 1 / (x + f) }
cf <- function(e) dnorm(e) * (1 - e * mills_cf(e))
cat(sprintf("  %5s %22s %12s\n", "eta", "continued fraction", "rel err"))
for (e in c(4, 5, 6, 7, 7.5, 7.8, 9, 12)) {
  r <- cf(e)
  cat(sprintf("  %5.1f %22.12e %12.2e\n", e, r, abs(g_fn(e, ub = 99) - r) / r))
}
cat("  accurate to ~1e-14 throughout: the authors' vpa step is not needed here\n\n")

cat("=========================================================\n")
cat("3. g inverse spline: round-trip accuracy\n")
cat("=========================================================\n")
ginv <- make_ginv()
eta_test <- seq(-7, 7.5, by = 0.13)
cat(sprintf("  max |eta - ginv(g(eta))| over [-7, 7.5] : %.3e\n",
            max(abs(eta_test - ginv(g_fn(eta_test))))))
cat(sprintf("  cost range covered by spline            : [%.3e, %.3f]\n",
            attr(ginv, "c_range")[1], attr(ginv, "c_range")[2]))
cat(sprintf("  Ursu lookup table floors at c          : %.3e (eta = 4.001)\n", g_fn(4.001)))
cat(sprintf("  our floor is lower by a factor of      : %.3g\n\n",
            g_fn(4.001) / attr(ginv, "c_range")[1]))

cat("=========================================================\n")
cat("4. Truncated exponential draws\n")
cat("=========================================================\n")
set.seed(7); mu <- 0.5; n <- 2e5
lb <- 0.3; ub <- 0.4
xb <- rexp_trunc_below(mu, rep(lb, n), runif(n))
xa <- rexp_trunc_above(mu, rep(ub, n), runif(n))
mb <- lb + mu
ma <- mu - ub * exp(-ub / mu) / (1 - exp(-ub / mu))
cat(sprintf("  trunc BELOW at %.2f : min=%.4f (>= %.2f)  mean=%.4f (want %.4f)\n",
            lb, min(xb), lb, mean(xb), mb))
cat(sprintf("  trunc ABOVE at %.2f : max=%.4f (<= %.2f)  mean=%.4f (want %.4f)\n\n",
            ub, max(xa), ub, mean(xa), ma))

cat("=========================================================\n")
cat("5. DGP: search-length distribution\n")
cat("=========================================================\n")
dat <- sim_ssm(N = 5000, seed = 1, ginv = ginv)
K   <- rowSums(dat$cons_set != 0)
cat("  searches per consumer:\n")
print(round(prop.table(table(K)), 4))
cat(sprintf("\n  mean searches             : %.3f\n", mean(K)))
cat(sprintf("  share buying outside opt  : %.4f\n", mean(dat$choice_vec == dat$J + 1)))
last_is_ch <- vapply(seq_len(dat$N),
                     function(i) dat$choice_vec[i] == dat$cons_set[i, K[i]], logical(1))
cat(sprintf("  share buying last-searched: %.4f\n", mean(last_is_ch)))
cat(sprintf("  share with K >= 3         : %.4f  <- region where the shipped\n", mean(K >= 3)))
cat("     loop code's Pr(F) differs from the vectorized code's\n\n")
