# Validation summary (bifactory)

This note summarizes numerical checks against Mplus. Full reproduction scripts
live in the package source repository under `cran/validation/`; that directory is
**not** included in the CRAN tarball (see `.Rbuildignore`).

## Continuous and ordinal single-group fits

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
| B-ESEM ordinal | BFI (25 items, 5 specific + G) | G ∈ {2, 3, 4, 5} | df match Mplus at all six levels; max \|ΔCFI\| ≤ 0.001; max \|ΔSRMR\| ≤ 0.001; standardized loadings within 0.001 at all six levels and in every group once the R solution is carried into Mplus's rotation orientation (model-implied polychoric correlations within 0.001 as well); group-wise Procrustes alignment of standardized loadings understates the agreement in non-reference groups from the weak level on (0.01 to 0.02), because their free factor covariances make the standardization rotation-dependent |
| B-ESEM ordinal | BAT 2023 (23 items, 4 specific + G) | G ∈ {2, 4, 8} | df match; see the paragraph below for what agrees and what does not |
| ESEM ordinal | BFI | G ∈ {2, 3, 4} | Theta-parameterization identification fixes |

**Caution:** with many groups the R side is slow (bfi: 19 / 40 / 100 minutes at
G = 3 / 4 / 5 on 25 items with the two starts of each level fitted one after the
other; eight BAT countries about three and a half hours; the default `cores = 2`
fits the two starts at the same time with identical results, 13 instead of 18
minutes at G = 3), and the
levels whose non-reference residual variances and factor covariances are empirically
unidentified on your data (a non-positive-definite factor covariance matrix in some
group) have no unique solution in any program---verify against external software
before reporting them.

WLSMV χ² difference tests from `lavTestLRT()` may differ numerically from Mplus
`DIFFTEST`. With identical df the WLSMV chi-square itself sits about 0.1 to 0.2
percent below Mplus's (2 to 4 units on bfi): both programs reach the same minimum
of the fit function and use the same scaling and shift, but lavaan weights each
group's discrepancy by n_g − 1 where Mplus uses n_g.

**What agreement with Mplus requires (BAT 2023 check, 23 items, up to 8 countries,
0.6.0).** lavaan ≥ 0.7-2 (earlier versions put sparse, highly correlated polychoric
pairs up to 0.03 too high); Mplus run at a convergence criterion of 1e-6 or tighter
(its default stopped 0.1 to 1 percent above its own optimum; the generated inputs now
use 1e-6); and the configural fit started from the per-group optima, which
`esem_invariance()` now does (the joint optimizer's path depends on the number of
groups and left one country in a local optimum 9 percent above Mplus's; the starts are
passed as a parameter table because lavaan caps a numeric `start` at 1000 entries and
four bfi tiers or five BAT countries already need more). With those
three in place, standardized loadings agree to about 0.01 after Procrustes alignment at
every level whose fit function has one optimum. Two things stay apart and are properties
of the data: levels with a non-positive-definite factor covariance matrix in some group
(both programs report it, each settles in its own point), and from six groups on the
latent-means model, which has two or three equally deep optima 0.11 to 0.16 apart;
which one a program reaches depends on its starting values (at eight and at five groups
lavaan's default start and Mplus's reach the same one, 0.011 and 0.006 apart, and
`esem_invariance()` keeps the default start's solution unless a second start is deeper
by more than 0.1 percent).

## Reproducing checks locally

```r
devtools::load_all("path/to/cran")
source("validation/bat3_lavaan_vs_mplus.R")  # requires Mplus output files
```

GitHub: <https://github.com/leondebeer/bifactory>
