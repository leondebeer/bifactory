# =============================================================================
#  bifactory -- Analysis Template
#  Open this file, copy it to your project, edit the CONFIG block, then run.
#
#    file.copy(system.file("templates", "template.R", package = "bifactory"),
#              "my_analysis.R")
#
#  Demo dataset: psych::bfi  (Big Five Inventory, 25 Likert items, N = 2800)
# =============================================================================
library(bifactory)

# -- 1. Load data -------------------------------------------------------------
data(bfi, package = "psych")
# psych::bfi has 5 factors x 5 items on a 6-point Likert scale + demographics

# -- 2. CONFIG -----------------------------------------------------------------
# Edit the factor-to-item mapping for your own data.
factors_list <- list(
  A = c("A1", "A2", "A3", "A4", "A5"),   # Agreeableness
  C = c("C1", "C2", "C3", "C4", "C5"),   # Conscientiousness
  E = c("E1", "E2", "E3", "E4", "E5"),   # Extraversion
  N = c("N1", "N2", "N3", "N4", "N5"),   # Neuroticism
  O = c("O1", "O2", "O3", "O4", "O5")    # Openness
)

use_ordered   <- TRUE           # TRUE -> WLSMV + polychoric; FALSE -> MLR
use_mplus     <- TRUE          # TRUE requires MplusAutomation + Mplus binary
output_folder <- "bfi_results"  # where save_results() will write

# Missing-data handling.  Applied uniformly to CFA, ESEM, and B-ESEM so all
# three fit indices are computed on the same sample.
#   NULL         -> "pairwise" (ordinal) or "listwise" (continuous) default
#   "listwise"   -> complete-case; always safe, drops rows with any NA
#   "pairwise"   -> ordered only; all available pairs for each polychoric
#   "fiml"/"ml"  -> continuous only; uses all rows under MAR. With >= 4
#                   specific factors, the FIML + GPArotation::targetT path
#                   used for B-ESEM can settle into a local optimum (G
#                   absorbs specific-factor variance); a warning fires so
#                   you can verify against other software such as Mplus, or
#                   fall back to "listwise".
missing_data  <- NULL


# -- 3. Alignment / misspecification check ------------------------------------
# Mehrvarz & Rouder (2026) ICM-CFA misspecification diagnostic -- run BEFORE
# fitting any model.  A high alignment-ratio dispersion indicates that the
# factor-to-item mapping is under-specifying cross-loadings, i.e. a standard
# CFA is inappropriate and ESEM / B-ESEM is warranted.
alignment <- alignment_check(bfi, factors_list)
print(alignment)
# plot(alignment)                              # per-pair visualisation


# -- 4. Specify model ----------------------------------------------------------
spec <- do.call(specify_model, c(
  factors_list,
  list(data    = bfi,
       ordered = use_ordered,
       missing = missing_data,
       label   = "BFI 5-factor")
))
print(spec)


# -- 5. Run the three-model comparison (CFA / ESEM / B-ESEM) -------------------
# run_comparison() fits all three models and builds a side-by-side fit table.
# Pass mplus_folder = "mplus_output/" to additionally generate .inp files.
# run_alignment = FALSE -- alignment was already computed in step 3.
results <- run_comparison(
  spec,
  mplus_folder  = if (use_mplus) "mplus_output" else NULL,
  run_alignment = FALSE
)
print(results)


# -- 6. Reliability indices (omega, omega_H, ECV, PUC, H) ----------------------
# Composite sums use absolute standardized loadings (Morin, Arens & Marsh
# 2016), so reverse-worded items need no reverse-scoring beforehand.
indices <- compute_indices(results)
print(indices)


# -- 7. Inspect individual fits ------------------------------------------------
# Standardised loadings with primary (target) items highlighted in blue.
# Works on fit_cfa / fit_esem / fit_besem, or on the whole pipeline.
parameters(results$fit_esem)                  # ESEM, all loadings
# parameters(results$fit_esem, suppress = 1.96)   # hide |z| < 1.96
# parameters(results)                         # all three models side by side

# Raw loading matrix (if you need it as a plain numeric matrix)
# round(std_loadings(results$fit_esem), 2)

# Fit statistics per model.
#   results$comparison_table -- the canonical, Mplus-matched view (recommended)
#
# `fit_indices()` returns CFI / TLI / RMSEA [90% CI] / SRMR as a formatted
# row, using the Mplus-matched values (.scaled for MLR/WLSMV, $wlsmv_stats
# for B-ESEM WLSMV).  Stack rows with rbind() + noquote() to compare models.
noquote(rbind(
  ESEM  = fit_indices(results$fit_esem),
  BESEM = fit_indices(results$fit_besem)
))

# A CFA of your own against the ESEM on the same data, like for like
# (same estimator, robust test, missing-data setting):
# esem_compare(results$fit_esem, "A =~ A1 + A2 + A3 + A4 + A5 \n C =~ C1 + ...")

# Reproducible orientation of the rotated loadings (sign + column order),
# or alignment to another fit of the same model (orientation transfer):
# align_loadings(results$fit_besem, target = "canonical")
# align_loadings(results$fit_besem, target = other_besem_fit)


# -- 8. Save tables + plots to disk -------------------------------------------
save_results(results, indices, output_folder = output_folder)


# =============================================================================
#  OPTIONAL 1 -- ESEM-within-CFA (EWC, Marsh et al. 2014)
#
#  Converts the rotated ESEM into an equivalent CFA model.  This is required
#  for most downstream extensions (latent regression, mediation, multi-group
#  with free factor variances, etc.) that standard ESEM rotations cannot fit.
# =============================================================================

# Auto-pick referent items (highest loading on each factor)
referents <- find_ewc_referents(results$fit_esem, spec)
print(referents)

# Generate the lavaan syntax (inspect / hand-edit if needed)
cat(ewc_syntax(results$fit_esem, spec))

# Fit the ESEM-within-CFA model
ewc_fit <- fit_ewc(results$fit_esem, spec)
print(ewc_fit)                            # compact CFI/TLI/RMSEA/SRMR summary

# Standardised loadings with primary (target) items highlighted
parameters(ewc_fit)                       # all loadings
# parameters(ewc_fit, suppress = 1.96)    # hide |z| < 1.96

# Full lavaan summary + the coloured loadings table
# summary(ewc_fit)


# =============================================================================
#  OPTIONAL 2 -- B-ESEM-within-CFA
#
#  Bifactor version of EWC.  Selects one referent per specific factor plus one
#  G-factor referent, fixes cross-loadings to their rotated values, and fits
#  the resulting CFA with orthogonal G + specific factors (Mplus BIFACTOR-
#  ESEM-WITHIN-CFA convention).
# =============================================================================

besem_referents <- find_ewc_referents(results$fit_besem, spec)
print(besem_referents)

# cat(ewc_syntax(results$fit_besem, spec))           # inspect generated syntax

besem_ewc_fit <- fit_ewc(results$fit_besem, spec)
print(besem_ewc_fit)
parameters(besem_ewc_fit)                          # coloured P/cross table
# summary(besem_ewc_fit)                           # full lavaan output + params


# =============================================================================
#  OPTIONAL 3 -- Factor scores
#
#  Export latent factor scores for each respondent (regression / Bartlett).
#  Works on ESEM, CFA, B-ESEM, and esem_invariance objects. Multi-group /
#  invariance results return a stacked data frame with a `group` column.
# =============================================================================

# Regression (BLUP) scores -- one row per respondent, one column per factor
scores_esem <- factor_scores(results$fit_esem, method = "regression")
head(scores_esem)

# Merge back to the original data
# bfi_with_scores <- cbind(bfi, scores_esem)

# Bartlett scores (unbiased; correct factor variance)
# scores_bart <- factor_scores(results$fit_esem, method = "bartlett")


# =============================================================================
#  OPTIONAL 4 -- Measurement invariance across groups
#
#  Requires `group` set in specify_model().  Fits configural / weak / strong /
#  strict levels and reports Delta CFI per step (Cheung & Rensvold 2002).
#  through = "varcov" or "means" adds levels 5 (latent variances/covariances)
#  and 6 (latent means) of Morin's sequence.  cores = 2 (default) fits the two
#  starting-value runs of each level on two worker processes; cores = 1 runs
#  them one after the other with identical results.
# =============================================================================

spec_mg <- do.call(specify_model, c(
  factors_list,
  list(data    = bfi,
       ordered = use_ordered,
       missing = missing_data,
       group   = "gender",
       label   = "BFI 5-factor (by gender)")
))

# Bifactor ESEM invariance (G + specific factors) -- the main attraction
# of this package. Swap to model = "esem" for the standard ESEM path.
# Validated against Mplus on BFI for G = 2 to 5 at all six levels (df
# identical, max |dCFI| = 0.001, standardized loadings within 0.001 once both
# solutions are expressed in one rotation orientation).
inv <- esem_invariance(spec_mg, model = "besem")
# inv <- esem_invariance(spec_mg, model = "esem")
# inv <- esem_invariance(spec_mg, model = "besem", through = "means")  # levels 5 + 6
# inv <- esem_invariance(spec_mg, model = "besem", cores = 1)          # sequential
print(inv)

# Partial invariance: release intercepts/thresholds (level = "strong") or
# residual variances (level = "strict") one at a time by score test, keeping
# two invariant anchors per factor; the printout lists every release.
# pinv <- partial_invariance(inv, level = "strong")

# Side-by-side check against Mplus (needs Mplus + MplusAutomation):
# run_mplus_besem_invariance(inv, "mplus_inv",
#   mplus_command = "C:/Program Files/Mplus/Mplus.exe",
#   group_labels = c("1" = "MALE", "2" = "FEMALE"))

# -- 4a. Access parameters at a specific invariance level --------------------
# inv$models has slots: configural, weak, strong, strict (plus varcov and
# means when fitted with through = "means"). Each is an esem_fit object
# that parameters() / summary() / fitMeasures() understand.

parameters(inv$models$strong)                   # STDYX loadings (strong level)
parameters(inv$models$strict)                   # STDYX loadings (strict level)
# summary(inv$models$strong, fit.measures = TRUE, standardized = TRUE)

# Save a level's loadings table to disk
ps_strict <- parameters(inv$models$strict, digits = 3)
# write.csv(ps_strict, "loadings_strict.csv", row.names = FALSE)


# -- 4b. Factor scores at a specific level -----------------------------------
# Auto-selected passing level: the most constrained fitted level (up to
# "means" when through = "means") with dCFI >= -0.010 at its own transition
# and an admissible solution (no negative variance, see inv$notes). Override
# with level = "configural" / "weak" / "strong" / "strict" / "varcov" /
# "means" to force a specific level.
scores_inv <- factor_scores(inv)                # auto level
# scores_strict <- factor_scores(inv, level = "strict")
# scores_bart   <- factor_scores(inv, method = "bartlett", level = "strong")

# Multi-group output prepends a `group` column. Save:
# write.csv(scores_inv, "scores_invariance.csv", row.names = FALSE)
# saveRDS(scores_inv, "scores_invariance.rds")


# -- End of template ----------------------------------------------------------
