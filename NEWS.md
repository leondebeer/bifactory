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
