# bifactory 0.6.0

Fixes and additions prompted by the review of the Teacher's Corner
manuscript. The minor version is bumped because reported values change:
SRMR now follows the Mplus definition, composite reliabilities use absolute
loadings, and lavaan >= 0.7-2 is required. Point estimates and fit of the
previously correct paths are unchanged. Prose now says "ordinal" instead of
"ordered-categorical"; function and argument names (`esem_ordered()`,
`ordered = TRUE`) are unchanged.

## Dependencies

* lavaan >= 0.7-2 is required. From that version the polychoric
  correlations equal the exact two-step maximum-likelihood values, which
  multi-group comparisons with Mplus depend on.

## Documentation

* Roxygen markdown is enabled, so the reference manual renders headings,
  bold text, bullet lists and tables. A test fails if raw markdown reaches
  the Rd files.

## Multi-group reference group

* `esem_invariance()` passes the sorted group levels to lavaan as
  `group.label`, so the reference group is the first sorted level in R and
  in the Mplus syntax written by `run_mplus_besem_invariance()` alike.
  Loadings and factor covariances are therefore oriented to the same group
  in both programs.

## Ordinal configural identification

* The ordinal configural model (ESEM and B-ESEM, theta parameterization)
  uses the standard identification in every group: residual variances 1,
  factor means 0, all thresholds free. Same df, equivalent model, faster
  convergence, and the generated Mplus syntax follows. Configural fit
  indices change by at most 0.001.

## ESEM-within-CFA

* `fit_ewc(var_fixed = FALSE)` re-expresses the source ESEM exactly and is
  honoured for a B-ESEM source (referent primary loadings fixed, factor
  variances free, covariances 0).
* The ordinal plain EWC is refitted under the source ESEM's
  parameterization, so freed factor variances read 1.00.
* The ordinal B-ESEM-within-CFA reports its own fit (df = B-ESEM df +
  k(k-1)/2); `print()` and `compare_ewc()` agree with each other and with
  Mplus.

## Reporting under MLR

* `print()` of ESEM, B-ESEM and EWC fits, `fit_indices()`, `compare_ewc()`,
  the `parameters()` header and `esem_compare()` report the robust
  chi-square, CFI, TLI and RMSEA, matching the comparison table.
* `esem_compare()` fits its CFA with the ESEM's estimator, robust test,
  standard errors, missing-data setting and, for ordinal items,
  parameterization, so both columns are like for like.

## Fit indices

* SRMR is computed as Mplus defines it (Asparouhov & Muthén, 2018) for
  every lavaan-based fit: correlation, variance and standardized-mean
  residuals for continuous items (denominator `p(p+1)/2 + p`), polychoric
  and category-probability residuals for ordinal items (denominator
  `p(p-1)/2 +` total categories), pooled over groups by sample size.
  Multi-group values agree with Mplus to three decimals at every invariance
  level of the validation runs; single-group values are unchanged at three
  decimals. The custom DWLS B-ESEM path uses the same denominator.

## Reliability

* `compute_indices()` sums absolute standardized loadings in every
  composite (Morin, Arens & Marsh, 2016), so reverse-scoring an item or
  reflecting a factor changes nothing (`test-reliability.R`). For CFA and
  ESEM the residual variance comes from the model-implied item variance
  (one minus the diagonal of Lambda Phi Lambda'), which accounts for the
  factor correlations.
* The subscale omega footnotes and the `compute_indices()` help describe
  the other factor's variance as an additional source of error rather than
  as partialled out (Morin et al., 2020).
* `print()` of the reliability indices follows the console width (80 to 100
  columns) instead of a fixed 100, and `alignment_check()` wraps its
  recommendation line.

## Measurement invariance

* The "Fitting <level> ..." line of `esem_invariance()` is printed in bold
  blue where the console supports ANSI colour (RStudio, or a terminal whose
  `TERM` is not "dumb"); `options(bifactory.ansi_colour = FALSE)` turns it
  off, `TRUE` forces it. Same detection as `parameters()`.
* `esem_invariance()` catches lavaan's estimation warnings per level, on the
  master and on the worker processes, and reports them reworded in
  `inv$notes` and under "Estimation notes" after the table. Admissibility is
  read from the reported solution, so a note names the group and the
  parameter, e.g. "a negative variance estimate (Heywood case): 2 (G latent
  variance -0.120)". Only the reported fit's notes are kept; package warnings
  are unaffected.
* A closing block ends every run and `print()`: levels fitted and
  considered, that the table can be interpreted, the most constrained
  admissible level supported by dCFI >= -0.010 among all fitted levels, any
  note at that level, a tip toward `through = "varcov"` / `"means"` or
  `partial_invariance()` when strict is supported but inadmissible, and
  "Remember: print(inv)". `factor_scores(level = "auto")` uses the same
  rule, so it never selects an inadmissible level and can select varcov or
  means when they were fitted.
* `esem_invariance()` gains `through = c("strict", "varcov", "means")`: the
  latent variance/covariance and latent mean invariance levels of Morin's
  sequence, for ESEM and B-ESEM under MLR and WLSMV. The default sequence
  is unchanged; `run_mplus_besem_invariance()` mirrors the fitted levels.
* `partial_invariance()` follows the partial-invariance rules for ESEM and
  B-ESEM. Only `level = "strong"` (intercepts or thresholds) and
  `level = "strict"` (residual variances) are accepted; `level = "weak"` is
  refused with an explanation, because the rotated loadings are constrained
  as a block. Constraints are released one at a time by score test; an item
  is released only if every factor it belongs to keeps at least two items
  with invariant intercepts/thresholds (Byrne, Shavelson & Muthén, 1989;
  van de Schoot, Lugtig & Hox, 2012); releasing more than 20% of the
  constraints added at a level warns (Dimitrov, 2010); every released
  parameter is listed. On ordinal items at least one threshold per item
  stays invariant (Millsap & Yun-Tein, 2004); an item's last threshold can
  be released by hand through `esem_invariance(..., group.partial =
  "item|t1")`. Released thresholds are honoured on the ordinal path, the
  continuous refit is retried from simple starting values, and B-ESEM is
  supported. A recovery test plants intercept shifts on two items in one of
  two simulated groups and checks that exactly those items are released, as
  intercepts on the continuous path and as thresholds on the ordinal path.
* The ordinal multi-group configural model is fitted from the per-group
  optima: each group is fitted alone from two starts and the joint fit
  starts from those values, passed as a parameter table. Every other level
  is fitted from lavaan's default and its `start = "simple"` starting
  values; the alternative is kept when its fit function is lower by more
  than 0.1%. The new `cores` argument (default 2) runs the two starts on
  two worker processes at the same time with identical results, so a level
  takes as long as its slower start; `cores = 1` fits sequentially.
* `align_loadings()` compares against another fit, or against a supplied
  solution (`lambda`, `psi`, `theta` per group), by orientation transfer:
  the unstandardized loadings are rotated to the reference and the factor
  covariances are transformed with the same rotation before standardizing.
  This is exact in every group, including groups with free factor
  covariances, where Procrustes alignment of standardized loadings is not.
  An `esem_invariance` target is read at the requested `level`.
* `run_mplus_besem_invariance()` writes `CONVERGENCE = 0.000001` for its
  first Mplus attempt; if Mplus does not converge there, the retries start
  at Mplus's default and then loosen (1e-4 to 5e-2). When a level's Mplus
  run prints no fit (singular information matrix, no DIFFTEST derivatives
  file), `DIFFTEST` is omitted from the next level's input.
  `extract_mplus_loadings()` returns numeric matrices for such outputs.

# bifactory 0.5.2

This release contains defensive fixes; estimation results are unchanged.

* `parameters()` on a Heywood-corrected ESEM fit now reports `NA` for standard
  errors, z-values, and p-values when the rotation changed after those
  quantities were computed, while retaining the corrected standardized
  loadings.
* Reliability calculations floor residual variances at 1e-6, preventing a
  Heywood case from producing a negative residual variance and corrupting the
  H index.
* `factor_scores()` now raises an informative error when the polychoric
  correlation or score-information matrix is singular, instead of exposing a
  bare `solve()` failure.
* Heywood retries now recognize the `NA`/0 target convention, so primary cells
  in targets built that way are freed correctly.
* B-ESEM sign correction now follows the target column order before orienting
  loadings, rather than assuming that the general factor is the first column.
* `factor_scores()` now warns when one pooled polychoric matrix is applied
  across multiple groups.
* `esem_ordered()` and `besem_ordered()` now stop with a clear message when
  supplied item names collide case-insensitively, instead of silently
  mismatching them.
* The ECV index now returns `NA` rather than `NaN` for a degenerate all-zero
  loading matrix.
* The RMSEA confidence-interval upper bound retries with a wider bracket only
  when the original search fails, leaving previously published values
  unchanged when the original search succeeds.

# bifactory 0.5.1

Changes in response to CRAN review (no changes to estimation results):

* DESCRIPTION: all acronyms are now spelled out on first use (WLSMV, MLR,
  DWLS, ICM-CFA); removed quotes around `efa()`.
* Console output outside `print()`/`summary()` methods now uses `message()`
  (or `warning()` where appropriate) instead of `cat()`/`print()`, so it can
  be suppressed with `suppressMessages()`. Progress output keeps its
  `verbose` gating.
* `parameters()` now returns a classed data frame (`bifactory_parameters`)
  that prints the same formatted table via its own `print()` method.
  Assigning the result no longer prints as a side effect.
* Random-number state is no longer saved to and restored from `.GlobalEnv`;
  internal fixed-seed sections use `withr::with_seed()` (new Import: withr).
  Draws are unchanged, so results are identical to 0.5.0.
* Examples rewritten to be self-contained and executable
  (`lavaan::HolzingerSwineford1939`); `\dontrun{}` is retained only for
  examples that require a licensed Mplus installation or open an interactive
  editor. Longer-running examples use `\donttest{}`. No runnable example
  uses a package from Suggests.

# bifactory 0.5.0

* Initial CRAN submission.

## Features

* Bifactor ESEM (B-ESEM), standard ESEM, and CFA for continuous and
  ordered-categorical data.
* Continuous models via 'lavaan' native `efa()` blocks (MLR).
* Ordered-categorical ESEM via 'lavaan' WLSMV; ordered B-ESEM via a custom
  diagonally weighted least squares (DWLS) path with polychoric correlations
  from 'psych', rotation-delta standard errors via 'numDeriv', and
  mean-and-variance-adjusted chi-square.
* Target, geomin, and oblimin rotations via 'GPArotation'; the bifactor ESEM
  approach follows Morin, Arens and Marsh (2016).
* Multi-group measurement invariance (configural through strict, with partial
  invariance), ESEM-within-CFA conversion, McDonald's omega reliability suite,
  and the Mehrvarz and Rouder (2026) alignment-ratio check for ICM-CFA
  misspecification.
* Optional side-by-side comparison with Mplus output via 'MplusAutomation'.
