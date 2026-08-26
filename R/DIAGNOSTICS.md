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

## What it actually was: finite-sample bias, quantified

The estimator is not broken. `beta` is essentially unbiased everywhere
(|bias| 0.007-0.06). The `(kappa, sigma_xi, a0, gamma0)` block carries a
persistent POSITIVE bias in `gamma0` that declines slowly in both sample size
and simulation draws.

Factorial, 8 replications per cell, data seeds held fixed across `nd` so only
the draw count varies within a row:

| tasks | nd=30 | nd=60 | nd=120 |
|---|---|---|---|
| 1200 | +0.686 | +0.401 | +0.597 |
| 2400 | +0.330 | +0.245 | +0.209 |

Sample size roughly halves the bias at every draw count. Draws help monotonically
at 2400 tasks. The non-monotonicity at 1200 is within one standard error — with
8 reps and gamma0 SD ~0.3 the SE per cell is ~0.1, which is exactly why single
cells looked erratic and produced several false diagnoses above.

Even at the best cell measured, gamma0 bias is +0.209 against a true −1.05:
**20% of the parameter value.** Reaching a negligible bias needs both a large
sample and many draws.

This is a known property of the literature rather than a SIFT defect.
@Yavorsky_2021 writes of his own simulations that "some of the average parameter
estimates are not within two standard errors of their true values. This is not
uncommon for search models estimated using SMLE," citing Honka (2014),
Ursu (2018) and Ursu et al. (2020). What is new here is the quantification.

### Earlier cells that looked anomalous were noise

A 12-replication cell at 2400 tasks / nd=60 gave 0.070 and looked clean; the
same cell with different seeds gives 0.149, as does 9600 tasks / nd=60. Nothing
was anomalous. The apparent non-monotonicity in sample size that prompted three
separate investigations was run-to-run variation being read as signal.

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
