# Shared setup for the comparison studies.
# Every experiment script is run FROM THE REPO ROOT and writes its result to
# comparisons/out/<name>.rds. The notebooks then load that file and display it,
# so the notebooks render in seconds and always show numbers that were actually
# produced rather than numbers that were typed.
source("R/sift_model.R")
GINV <- make_ginv()
OUT  <- "comparisons/out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

save_result <- function(name, obj) {
  obj$generated <- format(Sys.time(), "%Y-%m-%d %H:%M")
  saveRDS(obj, file.path(OUT, paste0(name, ".rds")))
  cat(sprintf("wrote %s/%s.rds\n", OUT, name))
  invisible(obj)
}

# ---------------------------------------------------------------------------
# A conjoint analyst's model: choice among ALL alternatives plus the outside
# option, with no notion of search. Conditional logit over the full attribute
# vector -- what you would fit if you ran a conjoint study on these products.
# ---------------------------------------------------------------------------
nll_conjoint <- function(p, d, searched_only = FALSE) {
  nb <- ncol(d$X); nk <- ncol(d$L)
  beta <- p[1:nb]; kappa <- p[(nb + 1):(nb + nk)]; a0 <- p[nb + nk + 1]
  V <- as.vector(d$X %*% beta) + as.vector(d$L %*% kappa)
  ll <- 0
  for (s in seq_len(d$S)) {
    ix <- ((s - 1L) * d$J + 1L):(s * d$J)
    keep <- if (searched_only) {
      k <- d$cons_set[s, ]; k <- k[k > 0]
      if (!length(k)) integer(0) else k
    } else seq_len(d$J)
    v <- c(a0, V[ix][keep])                       # outside option first
    ch <- d$choice[s]
    j <- if (ch == 0L) 1L else match(ch, keep) + 1L
    if (is.na(j)) next
    ll <- ll + v[j] - log(sum(exp(v - max(v)))) - max(v)
  }
  -ll
}

fit_conjoint <- function(d, searched_only = FALSE) {
  p0 <- rep(0, ncol(d$X) + ncol(d$L) + 1)
  optim(p0, nll_conjoint, d = d, searched_only = searched_only,
        method = "BFGS", control = list(maxit = 400, reltol = 1e-9))
}
