# Diagnostic: is the HB bias simulation bias (Jensen, O(1/nd)) or bad mixing?
# log of an unbiased probability estimate is biased DOWNWARD, and that bias
# shrinks with nd, NOT with sample size -- which matches a bias that persists
# (and in gamma0's case worsens) as T grows. Hold everything fixed but nd.
source("R/sift_hb.R"); ginv <- make_ginv()
TH_BAR <- c(1.0,-0.8,0.9,-0.6,0.4,log(0.35)); SD_I <- c(0.5,0.4,0.45,0.3,0.3,0.4)
PN <- c("beta1","beta2","kappa1","kappa2","a0","gamma0")
d <- sim_sift_het(N=150L, T=8L, J=5L, th_bar=TH_BAR, sd_i=SD_I, seed=108L, ginv=ginv)
info <- prep_sift(d)
out <- list()
for (nd in c(40L, 200L)) {
  cat(sprintf("\n===== nd = %d =====\n", nd))
  dr <- make_draws_sift(d, nd=nd, seed=508L)
  f <- sift_hb_par(d, info, dr, ginv, ncore=12L, niter=5000L, burn=2500L,
                   nq=12L, verbose=1000L, seed=8L)
  pm <- colMeans(f$theta_bar)
  cat(sprintf("%-8s %7s %9s %8s\n","param","truth","post_mean","bias"))
  for (k in 1:6) cat(sprintf("%-8s %7.3f %9.3f %+8.3f\n", PN[k], TH_BAR[k], pm[k], pm[k]-TH_BAR[k]))
  cat(sprintf("mean |bias| = %.3f   minutes = %.1f\n", mean(abs(pm-TH_BAR)), f$minutes))
  out[[as.character(nd)]] <- list(post=pm, bias=pm-TH_BAR, minutes=f$minutes,
                                  ind_cor=vapply(1:6,function(k) cor(f$theta_i_mean[,k], d$th_i[,k]),0))
  saveRDS(out, "R/out/diag_nd.rds")
}
cat("\n===== does the bias shrink with nd? =====\n")
for (k in 1:6) cat(sprintf("%-8s nd=40 %+7.3f   nd=200 %+7.3f\n", PN[k],
                           out[["40"]]$bias[k], out[["200"]]$bias[k]))
cat(sprintf("mean|bias|  nd=40 %.3f -> nd=200 %.3f\n",
            mean(abs(out[["40"]]$bias)), mean(abs(out[["200"]]$bias))))
