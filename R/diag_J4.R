# Enumeration + brute force at J = 4 (261 outcomes). The J = 3 test may be too
# small to expose an error in the selection system, which grows with the number
# of UNSEARCHED alternatives -- and the failing configs are J = 5 and J = 6.
source("R/sift_model.R"); g <- make_ginv()
J <- 4L
d0 <- sim_sift(N=1,T=1,J=J,seed=2,ginv=g)
th <- unpack(d0$true_par); sig_t <- d0$sig_t
Xb <- as.vector(d0$X %*% th$beta); mL <- as.vector(d0$L %*% th$kappa)
cst <- exp(th$gamma0 + th$gpos * d0$pos); u0 <- th$a0

perms <- function(v) if (length(v)<=1) list(v) else do.call(c, lapply(seq_along(v),
  function(i) lapply(perms(v[-i]), function(q) c(v[i], q))))
outs <- list(list(S=integer(0), ch=0L))
for (K in 1:J) for (S in unique(lapply(perms(seq_len(J)), function(q) q[seq_len(K)])))
  for (ch in c(0L,S)) outs[[length(outs)+1L]] <- list(S=S, ch=ch)
cat(sprintf("J=%d -> %d outcomes\n", J, length(outs)))

set.seed(1); NS <- 1500000; key <- character(NS)
for (i in seq_len(NS)) {
  eta <- rnorm(J); xi <- rnorm(J)*th$sig_xi
  dl <- Xb + eta; zz <- resv(dl, cst, sig_t, g); uu <- dl + mL + xi
  ord <- order(zz, decreasing=TRUE); best <- u0; K <- 0L
  for (k in 1:J) { if (zz[ord[k]] <= best) break; K <- k; best <- max(best, uu[ord[k]]) }
  S <- ord[seq_len(K)]
  key[i] <- paste0(paste(S, collapse="-"), "|", c(0L,S)[which.max(c(u0, uu[S]))])
}
tab <- table(key)/NS
tot <- 0; worst <- NULL
cat(sprintf("\n%-12s %10s %10s %8s %8s\n","outcome","model","brute","ratio","z"))
for (o in outs) {
  e <- d0; e$cons_set <- matrix(0L,1,J)
  if (length(o$S)) e$cons_set[1, seq_along(o$S)] <- o$S
  e$choice <- o$ch; inf <- prep_sift(e)
  pk <- exp(nll_sift_ghk(d0$true_par, e, inf, make_draws_sift(e,20000,seed=99), g, return_vec=TRUE))
  k <- paste0(paste(o$S, collapse="-"), "|", o$ch)
  pb <- if (k %in% names(tab)) as.numeric(tab[k]) else 0
  tot <- tot + pk
  se <- sqrt(max(pb,1e-9)*(1-pb)/NS)
  z  <- (pk-pb)/se
  if (pb > 0.002 || pk > 0.002) cat(sprintf("%-12s %10.5f %10.5f %8.2f %8.1f\n", k, pk, pb, if(pb>0) pk/pb else NA, z))
  if (pb > 0.001 && (is.null(worst) || abs(z) > abs(worst$z))) worst <- list(k=k, z=z, pk=pk, pb=pb)
}
cat(sprintf("\nmodel total = %.4f   brute total = %.4f\n", tot, sum(tab)))
cat(sprintf("worst z among outcomes with p>0.001: %s  z=%.1f (model %.5f vs brute %.5f)\n",
            worst$k, worst$z, worst$pk, worst$pb))
