# bifactory

**Mplus-aligned (bifactor) ESEM, CFA, and measurement invariance for continuous and ordered data in R.**

`bifactory` fits Exploratory Structural Equation Models (ESEM), bifactor ESEM
(B-ESEM), and Confirmatory Factor Analysis (CFA), and runs multi-group
measurement invariance — with numerical output aligned to Mplus conventions. It
combines EFA-style rotation (via **GPArotation** / **psych**) with structural
estimation (via **lavaan**), and handles both continuous data (MLR) and
ordered/Likert data (WLSMV).

The distinctive piece is the **ordered bifactor ESEM** path: a custom diagonally
weighted least squares (DWLS) routine over polychoric correlations, with an
orthogonal target rotation, a mean-and-variance-adjusted scaled chi-square, and
rotation-aware delta-method standard errors — the combination needed to
reproduce Mplus's `ROTATION = TARGET (ORTHOGONAL)` bifactor solution under WLSMV.

## Installation

```r
# install.packages("remotes")
remotes::install_github("leondebeer/bifactory")
```

Requires R (>= 4.1.0) and **lavaan (>= 0.6-21)** (0.6-22.2568+ recommended for
ordered polychoric parity with Mplus).

## Quick start

```r
library(bifactory)

# Continuous ESEM on the classic Holzinger-Swineford data
data("HolzingerSwineford1939", package = "lavaan")
d <- HolzingerSwineford1939[, paste0("x", 1:9)]

fit <- esem(d, nfactors = 3, factor_names = c("Visual", "Textual", "Speed"))

std_loadings(fit)                                   # standardized loadings (with cross-loadings)
factor_correlations(fit)                            # factor correlation matrix
fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr"))  # fit indices
```

## Bifactor ESEM for ordered (Likert) items

A general factor **G** loading on all items, plus orthogonal domain-specific
factors. For ordered indicators this uses the custom DWLS/WLSMV path:

```r
fit_b <- besem_ordered(
  data = my_items,
  specific_factors = list(
    EX = c("ex1", "ex2", "ex3", "ex4"),
    MD = c("md1", "md2", "md3"),
    CI = c("ci1", "ci2", "ci3")
  )
)
summary(fit_b, fit.measures = TRUE, standardized = TRUE)
compute_indices(fit_b)        # McDonald's omega suite (omega_H, ECV, PUC, H)
```

## Measurement invariance

```r
spec <- specify_model(
  EX = ex_items, MD = md_items, CI = ci_items,
  data = df, group = "sex", ordered = TRUE
)

inv <- esem_invariance(spec)  # configural -> weak -> strong -> strict
print(inv)                    # comparison table with delta-CFI / delta-RMSEA
```

Multi-group B-ESEM invariance (`model = "besem"`) is validated against Mplus for
small numbers of groups (the function warns at larger group counts, where the
optimiser can stall on the bifactor target rotation).

## What's in the box

| Function | Purpose |
|---|---|
| `specify_model()` | Define a model (factors, data, group, ordered) |
| `esem()` / `besem()` | ESEM / bifactor ESEM for continuous data (MLR) |
| `esem_ordered()` / `besem_ordered()` | ESEM / bifactor ESEM for ordered data (WLSMV) |
| `run_comparison()` | Fit CFA, ESEM, and B-ESEM together and compare |
| `esem_invariance()` | Multi-group configural / weak / strong / strict invariance |
| `compute_indices()` | McDonald's omega reliability suite |
| `make_target()` / `make_bifactor_target()` | Target rotation matrices |
| `alignment_check()` | Pre-fit ICM-CFA misspecification diagnostic |
| `std_loadings()`, `factor_correlations()`, `parameters()` | Tidy extractors |

## Optional: side-by-side comparison with Mplus

If you have Mplus and **MplusAutomation** installed, `run_comparison()` can
generate the `.inp` files, run them, and report R-vs-Mplus deltas. This is
entirely optional — every estimator runs in R alone — and engages only when you
supply an Mplus output folder.

```r
results <- run_comparison(spec, mplus_folder = "path/to/mplus_output")
```

## Citation

```r
citation("bifactory")
```

## Contributing

Issues and pull requests are welcome. Please open an issue to discuss larger
changes before submitting a PR.

## License

[AGPL-3](LICENSE). You may use `bifactory` (including commercially), but modified
or network-served versions must make their source available under the same
license.
