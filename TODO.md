# Working notes
Lifted out of the manuscript so the `.qmd` files read as submitted prose.
Nothing here is part of the paper.

## sections\analysis.qmd

```
<!--
TODO: awaiting guidance and data.

Proposed scope — the split-sample pilot described in the project brief:

  - Category, sample, and fielding partner.
  - The design: a SIFT block fielded alongside a conventional conjoint in a
    split sample, in a category with a genuine list-screen.
  - Descriptive evidence first: how much do respondents actually search,
    where do they click, does position move clicks.
  - Convergent validity — do the preference parameters from the search block
    agree with the conjoint part-worths. The most persuasive result
    available, and the one to lead with.
  - Is the estimated search cost distinguishable from zero.
  - Holdout choice prediction, search model versus conjoint.
  - Does position randomization move consideration in a way conjoint cannot
    capture.
  - Task-order effects: learning and fatigue as drift in the search cost.
-->
```

## sections\appendix-lit.qmd

```
<!-- TODO: cite the Sawtooth technical paper -->
```

```
<!-- TODO: vendor scan, then decide whether this belongs in the paper at all -->
```

```
<!--
TODO: prose. The four claims to be stated here as hypotheses the paper
tests, not as established gaps:

  1. A fieldable, productized survey instrument that yields Weitzman
     search-cost parameters the way conjoint yields part-worths. Lab
     analogues exist (@Gabaix_2006); a survey method with design norms
     and reporting conventions does not.
  2. Individual-level search costs identified from within-subject panel
     data rather than cross-sectional heterogeneity.
  3. Designed exogenous variation in position and information architecture
     in every task, rather than obtained once from a platform (@Ursu_2018).
  4. Segmentation on search behavior — how people shop — rather than on
     preference. No prior art known.

Each of these is a claim to disprove during the verification pass, not a
result. Do not write them as results.
-->
```

## sections\appendix.qmd

```
<!--
TODO:
  B. Derivation of the mixture-of-normals reservation value.
  C. Proof of the zero-search-cost nesting result (@sec-model-nesting).
  D. Sampler details and convergence diagnostics.
  E. Validation against the reference implementations (see @sec-estimation).
  F. Instrument screenshots and the full attribute list.
  G. Additional simulation results.
-->
```

## sections\conclusion.qmd

```
<!--
TODO: awaiting guidance. Draft last.

The closing argument, in one line: conjoint tells you what people prefer once
they are looking; SIFT tells you what they look at, what it costs them to
look, and what they never see at all. Because conjoint is the zero-search-cost
special case, adopting the instrument costs nothing when search frictions are
absent and buys a second set of parameters when they are not.
-->
```

## sections\discussion.qmd

```
<!--
TODO before circulating:
  - @sec-disc-estimation cites Honka_2014 and Ursu_2018 as reported by Yavorsky.
    Verify those are the papers he cites for SMLE bias rather than paraphrasing
    from memory.
  - @sec-disc-expectations promises a learning diagnostic in @sec-analysis, which
    does not exist yet. Keep the promise or cut the sentence.
  - the claim that ratios are more robust than levels to finite-sample bias is a
    conjecture. It is cheap to test in @sec-simstudy and should be tested before
    it is asserted here or in @sec-managerial-friction.
-->
```

## sections\identification.qmd

```
<!--
TODO: verify against Ursu (2018) §4.3 directly before asserting what SHE
concludes from conversion-on-rank. Our simulation establishes what the tests do
in OUR model; her data-generating environment differs (platform-assigned rather
than randomized position, one click per impression), and the selection channel
that drives our T1 result may operate differently there. Do not attribute the
sign pattern above to her.
-->
```

```
<!--
TODO once the simulation study exists (@sec-simstudy):
  - timings at commercial task counts with random coefficients, replacing the
    extrapolations in R/README.md with measured HB fits;
  - RMSE of the hybrid estimator against the kernel-smoothed one on identical
    data, which is the apples-to-apples number reviewers will want;
  - the recall-rate separation moment (Morozov Appendix D) computed on simulated
    data, to show it moves with the intercept and not with search cost.
-->
```

## sections\instrument.qmd

```
<!-- TODO: fix the partition for the empirical
application in @sec-analysis, and state it explicitly there. -->
```

```
<!-- TODO: confirm ten
against the fielding partner's length constraints; the identification
consequences are in @sec-estimation. -->
```

```
<!--
TODO: making a click cost something.

In a survey a click is nearly free, which is the central threat to
interpreting estimated search-cost levels (see @sec-discussion). Design
levers to consider, none yet chosen:
  - a short enforced latency on each product-page load, so a click costs
    real time rather than a mouse movement;
  - an explicit click budget per task;
  - incentive alignment in the sense of @Ding_2005 and @Ding_2007, so that
    time spent searching has an opportunity cost denominated in something
    the respondent wants;
  - denominating the cost in time directly and modeling it against
    observed dwell (@Ursu_2020).
Whatever is chosen has to be described here and defended in @sec-discussion.
Note that the enforced-latency option interacts with the timing data in
@sec-instrument-data: an imposed delay must be separable from voluntary
reading time.
-->
```

```
<!-- TODO: a figure showing one task's listing page and one product page, plus
a schematic of the event stream for a single respondent-task. Deferred until
the instrument is built. -->
```

## sections\intro.qmd

```
<!--
TODO before circulating:
  - the 64% attenuation figure and the sample-size claims point at @sec-simstudy,
    which is still being written from the factorial in R/out/. Reconcile the
    numbers once that section exists.
  - the empirical application does not exist yet; soften or cut the roadmap
    sentence on @sec-analysis if it is still absent at submission.
  - the behavioural-evidence paragraph cites no one. Pull the specific
    experimental papers from sections/litreview.qmd once their citations are
    checked, rather than asserting the pattern.
-->
```

## sections\model.qmd

```
<!--
DECIDED: carry L'kappa; estimate the variance-matched normal, report the exact
mixture as the generalization. Show the collapse into a scalar match value (giving
up on pricing product-page attributes altogether) as a misspecification case in
@sec-simstudy.

DECIDED: randomized attribute profiles, freshly drawn each task -- NOT a fixed set of
named alternatives recurring across tasks. Reasons to state in @sec-instrument:
  - independence across tasks is true by construction rather than assumed; a fresh
    profile really is a new object, whereas a recurring alternative would carry its
    match value forward and the model would have to absorb that as learning
    (cf. @Morozov_2021 §6.4, who estimate an MA(1) carryover of rho = 0.3 across
    sessions that are weeks apart -- ours are minutes apart, so carryover would be
    much larger);
  - heterogeneity then lives in beta_i over attributes, which is the conjoint HB
    object clients want, rather than in alternative-specific intercepts;
  - the nesting result needs attribute part-worths. Brand intercepts are not conjoint.
Note the generalization: a fixed recurring alternative set is the case
@Morozov_2021 and @Chung_2025 estimate, and their identification arguments apply to
it directly; ours are the attribute-level analogue.

Note for @sec-lit: @Ursu_2018 has exactly this list-page / product-page architecture
and needs no mixture, because she does not observe what the product page contains
(her §3.2.1) and so models it as an unobserved scalar match value. The difference is
not the architecture, it is that we design the page. Say so explicitly rather than
letting a referee find it.
-->
```

```
<!--
TODO: Proposition 1 is stated with a sketch rather than a proof. Either write the
formal argument into @sec-appendix and restore a pointer to it, or leave the sketch
and drop the claim to a formal proposition.
-->
```

## sections\simstudy.qmd

```
<!--
TODO: awaiting guidance. Run this before any client sees the instrument.

Proposed scope:

  - Parameter recovery at commercial scale: how many tasks per respondent,
    how many alternatives per list-screen, how many respondents.
  - Recovery of individual-level search costs specifically — the most
    optimistic claim in the paper, and the one most likely to fail.
  - Sensitivity to search-cost heterogeneity and to the strength of the
    randomization.
  - Cost of misspecification in both directions: fitting SIFT data with a
    conjoint model (search costs forced to zero), and fitting conjoint-like
    data with the search model.
  - Power to reject a zero search cost.
  - Behavioral robustness: data generated by satisficing or by rules of
    thumb, then fit with the Weitzman likelihood.
-->
```
