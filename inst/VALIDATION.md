# Validation summary (bifactory)

This note summarizes numerical checks against Mplus. Full reproduction scripts
live in the package source repository under `cran/validation/`; that directory is
**not** included in the CRAN tarball (see `.Rbuildignore`).

## Continuous and ordered single-group fits

| Dataset | N | Items | Factors | Checks |
|---------|---|-------|---------|--------|
| BAT3 | 495 | 23 | 4 (ESEM) | Loadings, SEs, CFI/TLI/RMSEA/SRMR (CFA, ESEM, B-ESEM; WLSMV and MLR) |
| neil.csv | 197 | 18 | 3 | Same |

Scripts: `validation/bat3_lavaan_vs_mplus.R`, `validation/neil_lavaan_vs_mplus.R`.

Typical tolerances after lavaan 0.6-22.2568: loadings max |Δ| ≤ 0.001; SEs ≤ 0.008
(max across test set); fit indices to 3 decimal places.

## Ordered B-ESEM

Default production path: `besem_ordered(method = "rotation")` (custom DWLS/WLSMV +
orthogonal target rotation). Validated against Mplus `ESTIMATOR = WLSMV;
ROTATION = TARGET (orthogonal)`.

`besem_ordered(method = "set-esem")` is a restricted lavaan CFA and is **not**
expected to match Mplus B-ESEM loadings.

## Multi-group measurement invariance (lavaan path)

| Model | Data | Groups validated | Notes |
|-------|------|------------------|-------|
| B-ESEM ordered | BFI (25 items, 5 specific + G) | G ∈ {2, 3, 4} | df match Mplus; max \|ΔCFI\| ≤ 0.001; max \|ΔSRMR\| ≤ 0.002 |
| ESEM ordered | Same | G ∈ {2, 3, 4} | Theta-parameterization identification fixes |

**Caution:** B-ESEM ordered invariance at **G ≥ 5** with ≥ 18 items is outside the
validated range; convergence or SE instability may occur---verify against external
software before reporting.

WLSMV χ² difference tests from `lavTestLRT()` may differ numerically from Mplus
`DIFFTEST`.

## Reproducing checks locally

```r
devtools::load_all("path/to/cran")
source("validation/bat3_lavaan_vs_mplus.R")  # requires Mplus output files
```

GitHub: <https://github.com/leondebeer/bifactory>
