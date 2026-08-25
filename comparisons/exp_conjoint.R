# =============================================================================
# COMPARISON: conjoint analysis vs SIFT
#
# Question. SIFT nests conjoint: as the search cost goes to zero every
# alternative is inspected and the model collapses to choice over fully realized
# utilities. So what does a conjoint analyst get wrong when search costs are NOT
# zero, and does the error have a predictable direction?
#
# Design. Sweep the search-cost intercept from nearly free to expensive. At each
# level simulate SIFT data and fit the model a conjoint analyst would fit: a
# conditional logit over the FULL choice set, using every attribute, with no
# notion of search.
#
# Comparison is on RATIOS of coefficients, not levels, because logit and probit
# scales differ. Ratios are scale-free and are what a conjoint analyst actually
# interprets (relative importance, willingness-to-pay).
#
# Prediction worth testing: attributes that are only visible AFTER a click
# (kappa) should be attenuated MORE than attributes visible on the landing page
# (beta), because conjoint assumes everyone saw everything while in fact the
# post-click attributes only influenced the choices of people who looked.
# =============================================================================
source("comparisons/_common.R")

BETA  <- c(1.0, -0.8)
KAPPA <- c(0.9, -0.6)
TRUE_RATIO <- c(beta2_beta1  = BETA[2] / BETA[1],
                kappa1_beta1 = KAPPA[1] / BETA[1],
                kappa2_beta1 = KAPPA[2] / BETA[1])

grid <- c(0.02, 0.06, 0.15, 0.35, 0.80, 2.00)   # mean search cost
rows <- list()

for (i in seq_along(grid)) {
  g <- grid[i]
  d <- sim_sift(N = 250, T = 8, J = 6, beta = BETA, kappa = KAPPA,
                gamma0 = log(g), seed = 400 + i, ginv = GINV)
  K <- rowSums(d$cons_set != 0)
  fit <- fit_conjoint(d)
  nb <- ncol(d$X)
  b <- fit$par[1:nb]; k <- fit$par[(nb + 1):(nb + ncol(d$L))]
  rows[[i]] <- data.frame(
    cost        = g,
    mean_clicks = mean(K),
    frac_seen   = mean(K) / d$J,          # share of the choice set inspected
    frac_nobuy  = mean(d$choice == 0),
    beta2_beta1  = b[2] / b[1],
    kappa1_beta1 = k[1] / b[1],
    kappa2_beta1 = k[2] / b[1])
}
res <- do.call(rbind, rows)

# error in each ratio, and the key contrast: post-click vs pre-click attenuation
res$err_beta2  <- res$beta2_beta1  - TRUE_RATIO["beta2_beta1"]
res$err_kappa1 <- res$kappa1_beta1 - TRUE_RATIO["kappa1_beta1"]
res$err_kappa2 <- res$kappa2_beta1 - TRUE_RATIO["kappa2_beta1"]
# attenuation factor: estimated ratio divided by truth, averaged over the two
# coefficients of each type. 1 = unbiased, <1 = attenuated toward zero.
res$atten_pre  <- (res$beta2_beta1 / TRUE_RATIO["beta2_beta1"])
res$atten_post <- rowMeans(cbind(res$kappa1_beta1 / TRUE_RATIO["kappa1_beta1"],
                                 res$kappa2_beta1 / TRUE_RATIO["kappa2_beta1"]))

print(round(res, 3))
save_result("conjoint", list(res = res, truth = TRUE_RATIO,
                             note = "conditional logit on the full choice set, N=250 x T=8, J=6"))
