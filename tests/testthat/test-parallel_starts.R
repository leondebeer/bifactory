# The two starts of every invariance level (lavaan default and start = "simple") are
# independent fits, so they run on two worker processes when a cluster is supplied;
# the result must be the fit the sequential rule picks.
data(HolzingerSwineford1939, package = "lavaan")
hs <- HolzingerSwineford1939
misfit <- "visual =~ x1 + x2 + x3\n speed =~ x7 + x8 + x9"

test_that(".make_cluster gives workers the package and the master's library path", {
  skip_on_cran()
  cl <- .make_cluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  expect_s3_class(cl, "cluster")
  # a closure over a package internal must run on the workers (installed or load_all'ed)
  f <- lavaan::cfa(misfit, data = hs)
  out <- parallel::parLapply(cl, list(list(lavaan_fit = f)), function(x) .inv_fx(x))
  expect_equal(out[[1]], .inv_fx(list(lavaan_fit = f)))
  lp <- parallel::clusterCall(cl, function() .libPaths()[1])
  expect_identical(lp[[1]], .libPaths()[1])
})

test_that(".fit_with_retry runs both starts on the cluster and applies the same rule", {
  skip_on_cran()
  cl <- .make_cluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  f <- lavaan::cfa(misfit, data = hs)
  # the alternative start is made 1 % deeper (a real well), so it must be chosen
  fit_fn <- function(ctrl = NULL, opts = NULL) {
    g <- f
    if (!is.null(opts$start)) g@optim$fx <- f@optim$fx * (1 - 1e-2)
    list(lavaan_fit = g, tag = if (is.null(opts$start)) "default" else "alt",
         pid = Sys.getpid())
  }
  par_fit <- .fit_with_retry(fit_fn, verbose = FALSE, cl = cl)
  seq_fit <- .fit_with_retry(fit_fn, verbose = FALSE)
  expect_identical(par_fit$tag, "alt")
  expect_identical(seq_fit$tag, "alt")
  expect_false(par_fit$pid == Sys.getpid())      # it was fitted on a worker
  expect_equal(.inv_fx(par_fit), .inv_fx(seq_fit))
})

test_that("esem_invariance(cores = 2) reproduces the sequential fit statistics", {
  skip_on_cran()
  hs_ord <- hs
  for (it in paste0("x", 1:9))
    hs_ord[[it]] <- cut(hs[[it]], quantile(hs[[it]], c(0, .25, .5, .75, 1)), include.lowest = TRUE, labels = FALSE)
  spec_o <- specify_model(Vis = c("x1", "x2", "x3"), Tex = c("x4", "x5", "x6"), Spd = c("x7", "x8", "x9"),
                          data = hs_ord, group = "school", ordered = TRUE, label = "HS ord")
  inv1 <- suppressWarnings(esem_invariance(spec_o, model = "besem", through = "strict", verbose = FALSE, cores = 1L))
  inv2 <- suppressWarnings(esem_invariance(spec_o, model = "besem", through = "strict", verbose = FALSE, cores = 2L))
  expect_equal(inv2$table$chisq, inv1$table$chisq, tolerance = 1e-6)
  expect_equal(inv2$table$df, inv1$table$df)
  for (lv in names(inv1$models))
    expect_equal(inv2$models[[lv]]$lavaan_fit@optim$fx, inv1$models[[lv]]$lavaan_fit@optim$fx, tolerance = 1e-8)
})
