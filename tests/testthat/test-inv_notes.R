# lavaan's raw warnings ("lavaan->lav_test_sb(): could not invert information
# matrix ...") must not reach the user during esem_invariance(): they are caught
# per level (on the master and on the workers), translated into plain notes,
# stored in inv$notes and printed after the table. Package warnings pass through.
data(HolzingerSwineford1939, package = "lavaan")
hs <- HolzingerSwineford1939
f0 <- lavaan::cfa("visual =~ x1 + x2 + x3\n speed =~ x7 + x8 + x9", data = hs)
lav_warn <- function() warning("lavaan->lav_test_sb():  \n   could not invert information matrix needed for robust test statistic", call. = FALSE)

test_that(".with_notes captures lavaan warnings and lets other warnings through", {
  r <- expect_no_warning(.with_notes({ lav_warn(); 1 }))
  expect_equal(r$value, 1)
  expect_match(r$notes, "robust test statistic could not be computed", fixed = TRUE)
  expect_warning(.with_notes(warning("mine")), "mine")
  expect_match(.translate_lavaan_warning("lavaan->lav_object_post_check():  \n   covariance matrix of latent variables is not positive definite;"),
               "not positive definite", fixed = TRUE)
})

test_that(".fit_with_retry returns the kept start's notes instead of raising, sequential and parallel", {
  skip_on_cran()
  fit_fn <- function(ctrl = NULL, opts = NULL) {
    if (is.null(opts$start)) lav_warn()          # only the default start warns
    list(lavaan_fit = f0, tag = if (is.null(opts$start)) "default" else "alt")
  }
  seq_fit <- expect_no_warning(.fit_with_retry(fit_fn, verbose = FALSE))
  expect_identical(seq_fit$tag, "default")        # equal fx: default kept
  expect_match(attr(seq_fit, "notes"), "robust test statistic", fixed = TRUE)
  cl <- .make_cluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  par_fit <- expect_no_warning(.fit_with_retry(fit_fn, verbose = FALSE, cl = cl))
  expect_identical(attr(par_fit, "notes"), attr(seq_fit, "notes"))
})

test_that("print.esem_invariance shows the notes and the conclusion line", {
  fits <- list(configural = list(lavaan_fit = f0))
  tbl <- .build_invariance_table(fits, list(), "configural", c(configural = "1. Configural"), FALSE)
  inv <- structure(list(table = tbl, models = fits,
                        lrt = list(), spec = list(label = "t", group = "g", group_levels = 1:2, ordered = NULL),
                        model = "esem", notes = list(strict = "robust test statistic could not be computed (singular information matrix)")),
                   class = "esem_invariance")
  out <- paste(capture.output(print(inv)), collapse = "\n")
  expect_match(out, "Estimation notes", fixed = TRUE)
  expect_match(out, "strict: robust test statistic", fixed = TRUE)
  expect_match(out, "concluded", fixed = TRUE)
  expect_match(out, "No level beyond configural is supported", fixed = TRUE)
  expect_match(out, "inv$models$configural", fixed = TRUE)
  concl <- paste(.inv_conclusion(list(models = list(configural = 1, weak = 1),
                                      table = data.frame(Model = c("1. Configural", "2. Weak (metric)"), dCFI = c(NA, 0)),
                                      notes = list(weak = "a negative variance estimate (Heywood case)"))), collapse = " ")
  expect_match(concl, "weak is supported by dCFI but its solution is inadmissible", fixed = TRUE)
  expect_match(concl, "No admissible level beyond configural", fixed = TRUE)
  concl2 <- paste(.inv_conclusion(list(models = list(configural = 1, weak = 1),
                                       table = data.frame(Model = c("1. Configural", "2. Weak (metric)"), dCFI = c(NA, 0)),
                                       notes = list(weak = "robust test statistic could not be computed (singular information matrix)"))), collapse = " ")
  expect_match(concl2, "most constrained level supported by dCFI >= -0.010 is: weak", fixed = TRUE)
  expect_match(concl2, "parameters(inv$models$weak)", fixed = TRUE)
  expect_match(concl2, "Note at that level: robust test statistic", fixed = TRUE)
})

test_that(".inv_admissibility names the group and parameter of every negative variance and non-PD psi", {
  pt <- data.frame(lhs = c("G", "ODI3", "G", "ODI3", "G", "G"), op = c("~~", "~~", "~~", "~~", "=~", "~~"),
                   rhs = c("G", "ODI3", "G", "ODI3", "ODI3", "F1"), group = c(1, 1, 2, 2, 2, 2),
                   est = c(1, 1, -0.12, -0.03, 0.5, 0.2), stringsAsFactors = FALSE)
  psi <- list(diag(2), matrix(c(1, 1.2, 1.2, 1), 2))    # group 2 not positive definite
  out <- .inv_admissibility_pt(pt, c("BE", "NL"), psi)
  expect_match(out[1], "negative variance estimate (Heywood case): NL (G latent variance -0.120; ODI3 residual variance -0.030)", fixed = TRUE)
  expect_match(out[2], "latent covariance matrix (psi) is not positive definite: NL", fixed = TRUE)
  expect_length(.inv_admissibility_pt(pt[pt$group == 1, ], "BE", psi[1]), 0)
  expect_length(.inv_inadmissible(out), 2)
  data(HolzingerSwineford1939, package = "lavaan")
  expect_length(.inv_admissibility(lavaan::cfa("visual =~ x1 + x2 + x3", data = HolzingerSwineford1939, group = "school")), 0)
})

test_that(".inv_supported_level picks the most constrained level whose own dCFI passes", {
  tbl <- data.frame(Model = c("1. Configural", "2. Weak (metric)", "3. Strong (scalar)", "4. Strict"),
                    dCFI = c(NA, -0.002, -0.020, -0.001), stringsAsFactors = FALSE)
  models <- list(configural = 1, weak = 1, strong = 1, strict = 1)
  expect_identical(.inv_supported_level(tbl, models), "strict")
  expect_identical(.inv_supported_level(tbl, models[c("configural", "weak", "strong")]), "weak")
  # every fitted level counts, including varcov and means (through = "means")
  tbl6 <- rbind(tbl, data.frame(Model = c("5. Latent var/cov", "6. Latent means"), dCFI = c(-0.001, -0.047)))
  models6 <- c(models, list(varcov = 1, means = 1))
  expect_identical(.inv_supported_level(tbl6, models6), "varcov")
  tbl6$dCFI[6] <- 0
  expect_identical(.inv_supported_level(tbl6, models6), "means")
  concl6 <- paste(.inv_conclusion(list(models = models6, table = tbl6, notes = list())), collapse = " ")
  expect_match(concl6, "Levels considered: configural, weak, strong, strict, varcov, means", fixed = TRUE)
  # an inadmissible solution (negative variance, non-PD matrix) is never recommended
  heywood <- list(strict = "a negative variance estimate (Heywood case)")
  expect_identical(.inv_supported_level(tbl6, models6[1:4], notes = heywood), "weak")
  expect_identical(.inv_supported_level(tbl, models, notes = heywood), "weak")
  npd <- list(strict = "latent covariance matrix (psi) is not positive definite in at least one group",
              weak = "robust test statistic could not be computed (singular information matrix)")
  expect_identical(.inv_supported_level(tbl, models, notes = npd), "weak")   # test-statistic note is not inadmissibility
  tbl$dCFI <- c(NA, -0.02, -0.02, -0.02)
  expect_identical(.inv_supported_level(tbl, models), "configural")
  expect_identical(.inv_supported_level(tbl, models, cutoff = -0.03), "strict")
  concl <- paste(.inv_conclusion(list(models = models, table = tbl, notes = heywood), cutoff = -0.03), collapse = " ")
  expect_match(concl, "strict is supported by dCFI but its solution is inadmissible (a negative variance estimate (Heywood case))", fixed = TRUE)
  expect_match(concl, "supported by dCFI >= -0.030 is: strong", fixed = TRUE)
  # strict inadmissible and levels 5/6 not fitted: point to through = "varcov" / "means"
  expect_match(concl, 'Tip: fit through = "varcov" or through = "means"', fixed = TRUE)
  concl_v <- paste(.inv_conclusion(list(models = c(models, list(varcov = 1)), table = tbl, notes = heywood), cutoff = -0.03), collapse = " ")
  expect_no_match(concl_v, "Tip:", fixed = TRUE)
})

test_that("factor_scores(level = 'auto') skips an inadmissible level", {
  skip_on_cran()
  data(HolzingerSwineford1939, package = "lavaan")
  spec <- specify_model(Vis = c("x1", "x2", "x3"), Tex = c("x4", "x5", "x6"), Spd = c("x7", "x8", "x9"),
                        data = HolzingerSwineford1939, group = "school", label = "HS")
  inv <- suppressWarnings(esem_invariance(spec, model = "esem", through = "strict", verbose = FALSE, cores = 1L))
  auto <- suppressMessages(.fs_select_invariance_level(inv, "auto"))$level
  inv$notes[[auto]] <- "a negative variance estimate (Heywood case)"
  auto2 <- suppressMessages(.fs_select_invariance_level(inv, "auto"))$level
  expect_false(identical(auto2, auto))
  expect_identical(auto2, .inv_supported_level(inv$table, inv$models, notes = inv$notes))
})
