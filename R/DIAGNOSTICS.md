# Diagnosing the search-cost block

A record of a long debugging session, kept because the wrong turns are as
instructive as the fix, and because anyone who sees odd estimates from a
sequential-search likelihood should check the same things in the same order.

## The symptom

Estimates of the preference parameters `beta` were clean at every sample size.
Estimates of `(kappa, sigma_xi, a0, gamma0)` — the post-search spread and the
search cost — were badly off, and the error flipped sign across datasets:
gamma0 errors of −0.03, −0.77, +0.64, +0.19 on four runs of identical code.

## What it was NOT

Each of these was proposed, tested, and rejected. Four of the five were
proposed *before* being tested, which is the actual lesson of this document.

| hypothesis | test | verdict |
|---|---|---|
| MCMC not mixing | start the HB sampler at the truth | rejected — cold and warm chains land identically (gamma0 bias +0.613 vs +0.612) |
| flat `sigma_tilde`/`gamma0` ridge | sweep shifter strength | rejected — bias flat through a 20× cost gradient |
| too few alternatives (J) | J ∈ {5,6,7} at fixed everything else | rejected — J=6 was clean at 800 tasks and biased at 1200 |
| model misspecified vs the DGP | enumerate all outcomes at J=3 (49) and J=4 (261) against brute force | rejected — matches within Monte Carlo error, sums to 1 across the parameter space |
| simulation (Jensen) bias | log-likelihood vs `nd` at fixed parameters | rejected as the cause — 0.019 nll/task at nd=40 |

## What it was, in part: GHK underflow

`ghk()` accumulated its importance weight as a running **product** of up to J
probabilities. For an outcome unlikely under the current parameters that product
underflows to zero, and a floored zero contributes `-log(1e-300) ≈ -691` to the
log-likelihood.

One such task in 1200 was **21.8% of the total nll**. The optimiser was not
fitting data in those regions — it was fleeing underflow.

**Fix:** accumulate in log space (`pnorm(..., log.p = TRUE)`), draw the truncated
normals through the upper tail so tiny probabilities stay representable, and
combine with `logmeanexp`. See `ghk()` in `sift_model.R`.

**Effect:** nll at the truth improved by 691 units, floor incidence went to zero,
and outright failed fits (beta1 = 0.47) disappeared.

**But this did not remove the bias.** It is a real bug that had to be fixed; it
was not the explanation.

## What it actually was: finite-sample bias

The estimator is consistent. Search-block bias falls roughly like 1/N:

| tasks | search-block \|bias\| | gamma0 bias |
|---|---|---|
| 1200 | 0.446 | +0.709 |
| 9600 | 0.098 | +0.124 |

`beta` is unbiased at every size (t < 1). The `(kappa, sigma_xi, a0, gamma0)`
block is weakly identified and needs substantially more data — which is a
sample-size recommendation for the paper, not a defect.

## Rules this session earned

1. **Never diagnose bias from a single dataset.** On this block the error flips
   sign across draws. The replication harness should be the first thing built,
   not the fifth.
2. **Validate away from the truth.** Every test in the original suite ran at or
   near the true parameters, where underflow is rare. The bug only bit where the
   optimiser roams. `test_sift_model.R` now checks floor incidence.
3. **Report convergence codes and per-replication estimates**, not just means.
   One failed fit in twelve masqueraded as a sample-size effect.
4. **A standard deviation that grows with sample size is an optimiser problem**,
   never a statistical one.
