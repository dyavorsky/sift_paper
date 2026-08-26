# Decisive diagnostic: start the sampler AT THE TRUTH.
#   stays put  -> the posterior is fine; the cold-start chain simply had not
#                 mixed, and the fix is longer burn-in / better proposals
#   drifts off -> the posterior genuinely peaks elsewhere, which is either
#                 misspecification or the sigma_tilde / gamma0 ridge biting
source("R/sift_hb.R"); ginv <- make_ginv()
TH_BAR <- c(1.0,-0.8,0.9,-0.6,0.4,log(0.35)); SD_I <- c(0.5,0.4,0.45,0.3,0.3,0.4)
PN <- c("beta1","beta2","kappa1","kappa2","a0","gamma0")
d <- sim_sift_het(N=150L, T=8L, J=5L, th_bar=TH_BAR, sd_i=SD_I, seed=108L, ginv=ginv)
info <- prep_sift(d); dr <- make_draws_sift(d, nd=100L, seed=508L)

runs <- list(
  cold  = list(init=NULL,      common=NULL),
  warm  = list(init=d$th_i,    common=c(log(d$sig_xi), d$gpos)))   # start at truth

out <- list()
for (nm in names(runs)) {
  cat(sprintf("\n===== %s start =====\n", nm))
  r <- runs[[nm]]
  f <- sift_hb_par(d, info, dr, ginv, ncore=12L, niter=6000L, burn=3000L, nq=12L,
                   init=r$init, common_init=r$common, verbose=1500L, seed=11L)
  pm <- colMeans(f$theta_bar)
  cat(sprintf("%-8s %7s %9s %8s\n","param","truth","post","bias"))
  for (k in 1:6) cat(sprintf("%-8s %7.3f %9.3f %+8.3f\n", PN[k], TH_BAR[k], pm[k], pm[k]-TH_BAR[k]))
  cat(sprintf("mean|bias| = %.3f   mean logL = %.1f\n", mean(abs(pm-TH_BAR)), mean(f$logL)))
  out[[nm]] <- list(post=pm, bias=pm-TH_BAR, logL=f$logL,
                    ind_cor=vapply(1:6,function(k) cor(f$theta_i_mean[,k], d$th_i[,k]),0))
  saveRDS(out, "R/out/diag_start.rds")
}
cat("\n===== verdict =====\n")
cat(sprintf("mean|bias|  cold %.3f   warm %.3f\n",
            mean(abs(out$cold$bias)), mean(abs(out$warm$bias))))
cat(sprintf("mean logL   cold %.1f   warm %.1f  (higher = better fit)\n",
            mean(out$cold$logL), mean(out$warm$logL)))
cat("If warm keeps a much better logL AND lower bias, the cold chain never mixed.\n")
cat("If both land in the same place, the posterior itself is off.\n")
