# R port of the Chung–Chintagunta–Misra sequential-search SMLE

This directory contains an R port of the "PM" (probability-mapping) estimator from
Chung, Chintagunta & Misra (2025), together with the tests that establish it is
correct. The port exists to answer one question for SIFT: **can a Weitzman search
model be estimated inside a 12-hour budget?** The PM estimator was the leading
candidate; this is the check that it works before we build on it.

Source material: `C:\Users\DanYavorsky\Dropbox\resources\Chung Chintagunta Misra -
Matlab Code for Seq Search SMLE\`, downloaded from the authors' site. The port
follows `ssm_lik_for_loop.m` / `main_for_loop.m` — the readable loop version.

## Files

| file | purpose |
|---|---|
| `chung_pm.R` | the port: `g_fn`, `make_ginv`, truncated-exponential draws, `sim_ssm` (DGP), `prep_data`, `make_draws`, `nll_loop` (readable oracle), `prep_groups` + `nll_grouped` (batched, the one to use) |
| `test_chung_pm.R` | five component checks (see below) |
| `test_probability.R` | the decisive check: enumerate all 48 outcomes at J=3, verify they sum to 1 and match a brute-force DGP simulation |
| `run_replication.R` | 100-replication Monte Carlo against Table 1, PM column |
| `analyze_replication.R` | the comparison table |

## What the port does differently, and why

**No `vpa`.** The authors compute `g(η) = φ(η) − η(1−Φ(η))` in variable-precision
arithmetic because the two terms cancel catastrophically in the tail. In R the
direct formula is accurate to ~1e-14 out to η = 12, verified against a 500-level
Mills-ratio continued fraction that uses neither `pnorm` nor any series. An
asymptotic-expansion branch was tried and **removed**: the Mills series is
divergent, and it was *less* accurate than the direct formula (rel. error 1e-6 at
η = 7.8). Do not reintroduce it.

**The g-inverse spline is built in the standardised argument `c/σ`,** so one spline
serves every σ rather than being rebuilt per parameter draw. It covers costs down
to 3.85e-16, against the 7.11e-06 floor of the Ursu/Seiler/Honka lookup table —
lower by a factor of 1.85e10, which matters because the likelihood evaluates `g⁻¹`
at whatever the optimiser proposes, not at economically sensible values.

**Exponential search costs are exploited directly.** Memorylessness makes
truncation-from-below a shift, and makes every unsearched-alternative term a single
division rather than a CDF evaluation.

## Validation

All five pass. Run `Rscript R/test_chung_pm.R` and `Rscript R/test_probability.R`.

1. **`g()` against Matlab ground truth.** Compared to `tableZ.csv` from the
   Ursu/Seiler/Honka archive: **100.00% of 7551 rows identical at 5 significant
   digits**, max reldiff 0.000e+00. (The apparent 4.99e-05 mismatch on first
   inspection was `csvwrite`'s 5-digit storage, not a discrepancy.)
2. **Far-tail accuracy** ~1e-14 through η = 12, per the continued fraction above.
3. **Spline round-trip** `max |η − g⁻¹(g(η))| = 5.3e-15` over [−7, 7.5].
4. **`nll_loop` ≡ `nll_grouped`** to machine precision (|diff| ≤ 4.55e-13). The
   batched version is what the optimiser calls; the loop version is the readable
   statement of the same thing, kept as an oracle.
5. **The likelihood is an unbiased estimator of the true outcome probabilities.**
   Enumerating all 48 outcomes at J = 3, the simulated probabilities sum to 1 and
   converge to brute-force frequencies at the 1/√nd rate (SD fell 7.6× over a 64×
   increase in draws). At nd = 32000 the model gives 0.131583 against brute force's
   0.131625.

   *Methodological note:* an earlier version of this test scored deviations against
   brute-force noise alone and produced a spurious z = −12.58. The simulated
   likelihood's own Monte Carlo error is ~7× larger; the test now replicates the
   model side and combines both error sources.

## Two bugs in the authors' shipped code

**1. `ssm_lik_vec.m` cannot run as shipped.** It calls
`new_exprnd_trunc_above` / `new_exprnd_trunc_below`, which are not in the archive.
(The loop version calls `cpu_exprnd_trunc_*`.) `main_vectorize.m` therefore fails
immediately. Documented shims are added in the reference directory; see
`README_dy.md` there.

**2. The shipped loop code computes Pr(F) with the wrong reservation value.**
`ssm_lik_for_loop.m` uses `z_next_H` — which equals z₂ once K ≥ 3 — where it should
use `z_K_H`. The vectorized code and the paper both use z_K. z_K is correct: the
binding constraint is the *last* continuation condition, and z is decreasing along
the search order. This affects the **12.76%** of consumers with K ≥ 3. At the true
parameters the buggy version reports a lower nll (3033.57 vs 3101.27), i.e. it
overstates the likelihood.

Both `nll_loop` and `nll_grouped` take `pr_F_use = c("zK", "znext")` so the
replication can quantify what the bug does to the estimates. **Use `"zK"`.**

## Estimation notes

Single dataset (N = 1000, J = 5), solo on this machine: **77 s at nd = 100**, 451 s
at nd = 500, with near-identical estimates — residual error is finite-sample, not
simulation noise, so nd = 100 is sufficient.

`run_replication.R` fits in two stages (nd = 25, then nd = 100). This does not
speed up a solo fit but it converges cleanly in both stages. Two settings matter:
`reltol` must not be tighter than the simulation noise floor or BFGS grinds to
`maxit` chasing numerical wiggle, and results are written per replication so a
partial run is still usable and the job is resumable.

Profiling: `pnorm` 51% + `dnorm` 24% = **76% of runtime**, essentially all of it
inside `g_fn` (65% of total). The R code is near its floor without Rcpp — which is
the obvious next lever if SIFT needs more speed.

---

## Replication result (100 datasets, N = 1000, J = 5, nd = 100)

`Rscript R/run_replication.R && Rscript R/analyze_replication.R`. 34.6 min on 14
cores; **100/100 converged** in both Pr(F) modes.

| param | true | Chung Table 1 (PM) | this port |
|---|---|---|---|
| α₀₂ | −0.500 | −0.500 (0.062) | −0.494 (0.059) |
| α₀₃ | −0.800 | −0.799 (0.068) | −0.808 (0.068) |
| α₀₄ | 0.500 | 0.500 (0.053) | 0.497 (0.051) |
| α₀₅ | 0.800 | 0.796 (0.051) | 0.799 (0.055) |
| α₁₁ | 1.500 | 1.495 (0.063) | 1.498 (0.054) |
| α₁₂ | −1.000 | −0.999 (0.066) | −1.008 (0.085) |
| μ₀ | 2.500 | 2.524 (0.162) | 2.522 (0.150) |
| log c̄ | −0.693 | −0.682 (0.052) | −0.687 (0.049) |

Centres and sampling spreads both match. Largest |t| for bias across the eight
parameters is **1.50** — no detectable bias anywhere.

**On their headline RMSE of 0.040.** We get 0.0729. This is a definitional gap, not
a failure: their *own* reported SDs imply a per-dataset RMSE of
√(mean SD²) = **0.080**, i.e. twice the number in their table. Our SDs imply 0.078,
and our measured 0.0729 sits just below that (Jensen). So the port is marginally
*tighter* than the paper at the parameter level; we simply cannot reproduce that one
summary statistic without knowing how they formed it.

**What the Pr(F) bug costs.** Estimating with the shipped loop code's `z_next`:

| param | zK (correct) | znext (shipped) | shift |
|---|---|---|---|
| log c̄ | −0.687 | −0.569 | **+0.118** |
| μ₀ | 2.522 | 2.617 | +0.095 |
| α₁₁ | 1.498 | 1.545 | +0.046 |

Mean RMSE rises 35% (0.0729 → 0.0986), and the bias in log c̄ is **26.6 standard
errors** — it inflates the estimated mean search cost by 12% (0.503 → 0.566).
Search cost is the object SIFT exists to measure, so this is exactly the parameter
we could least afford to get wrong.

## Scaling (`R/bench_scaling.R`, solo)

Per likelihood evaluation:

* **N: exactly linear** — 0.26 s per 1000 observations (J = 5, nd = 100), flat from
  N = 500 to N = 8000.
* **nd: exactly linear** — 0.24 s per 100 draws (N = 1000, J = 5).
* **J: sublinear and saturating** — 0.147 s at J = 3, 0.475 at J = 8, 0.696 at
  J = 16. It flattens because search stops early: the number of *searched*
  alternatives saturates even as the choice set grows.

So `sec/eval ≈ 0.26 · (N/1000) · (nd/100) · f(J)`, with f(5) = 1, f(8) ≈ 1.9,
f(12) ≈ 2.7.

Evaluation counts, backed out from the fits: `evals ≈ n_fn + p · n_grad` with
n_fn ≈ 50 and n_grad ≈ 25 per stage, i.e. **≈ (50 + 25p) per stage**, two stages.
The p term dominates — parameter count, not sample size, drives the eval count.

### What this means for SIFT's 12-hour budget

Assume 800 respondents × 10 tasks = **8000 task observations**, J = 8 per landing
page, p = 20 (attribute levels plus search cost), nd = 100 with an nd = 25 first
stage. Then sec/eval ≈ 4.0 (stage 2) and ≈ 1.0 (stage 1), and ≈ 550 evals/stage:

| specification | fit time | basis | verdict |
|---|---|---|---|
| homogeneous, SMLE | ~45 min | extrapolated | comfortable |
| random coefficients, SMLE, 100 mixing draws | ~75 h | extrapolated | **infeasible in R** |
| hierarchical Bayes, 20k iterations, single-threaded | **~80 h** | **measured** | **infeasible** |
| hierarchical Bayes, 20k iterations, 12 cores | **~11 h** | **measured** | feasible |

The HB rows are measured rather than extrapolated, and they corrected an earlier
estimate. `sift_hb()` costs **1.8 ms per task per iteration** — linear in tasks,
and only weakly increasing in `nd` because per-task overhead dominates the draw
count. An earlier extrapolation of ~11 h single-threaded was wrong by a factor
of seven: it assumed one pass over the data per iteration, where the sampler
makes two and carries per-respondent overhead on top.

Parallelising the `theta_i` step (`sift_hb_par()`) measures a **7.5× speedup on
12 cores** at 640 tasks, which is conservative — at that size communication is a
material share of each iteration, and the share falls as the problem grows. That
is what brings a commercial-scale fit inside an overnight budget. **HB meets the
12-hour target only in parallel.**

Two consequences worth carrying into the paper. First, the homogeneous model is not
the constraint — heterogeneity is, and it is the whole point, since recovering
individual-level preferences is what makes SIFT a conjoint substitute. Second, on
these numbers **HB is the cheaper route to heterogeneity than SMLE with mixing
draws**, because HB never integrates over the heterogeneity distribution: each
respondent is evaluated at their own θ, so the mixing-draw multiplier disappears.
That inverts the usual ordering and supports doing HB properly rather than as a
robustness check.

*These are extrapolations from measured per-eval costs and measured eval counts, not
timed fits of the SIFT model, which does not exist yet.*
