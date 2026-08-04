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
