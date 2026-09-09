# The "Fitting <level> ..." line of esem_invariance() is printed in bold blue
# where ANSI colour is available (same detection as parameters()), so the
# start of each level stands out in a long log; plain text otherwise.
test_that(".ansi_colour follows the bifactory.ansi_colour option", {
  withr::with_options(list(bifactory.ansi_colour = TRUE), expect_true(.ansi_colour()))
  withr::with_options(list(bifactory.ansi_colour = FALSE),
    withr::with_envvar(c(RSTUDIO = "", TERM = "dumb"), expect_false(.ansi_colour())))
})

test_that("a verbose esem_invariance() run ends with the print(inv) reminder", {
  skip_on_cran()
  data(HolzingerSwineford1939, package = "lavaan")
  spec <- specify_model(Vis = c("x1", "x2", "x3"), Tex = c("x4", "x5", "x6"),
                        data = HolzingerSwineford1939, group = "school", label = "HS")
  msgs <- character()
  withCallingHandlers(
    suppressWarnings(esem_invariance(spec, model = "esem", cores = 1L)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_identical(trimws(msgs[length(msgs)]), "Remember: print(inv)")
})

test_that(".fit_header colours the Fitting line only when colour is on", {
  withr::with_options(list(bifactory.ansi_colour = TRUE), {
    h <- .fit_header("2. Weak (metric)")
    expect_match(h, "\033[1;34m", fixed = TRUE)
    expect_match(h, "Fitting 2. Weak (metric) ...", fixed = TRUE)
    expect_match(h, "\033[0m", fixed = TRUE)
  })
  withr::with_options(list(bifactory.ansi_colour = FALSE),
    withr::with_envvar(c(RSTUDIO = "", TERM = "dumb"),
      expect_identical(.fit_header("4. Strict"), sprintf("  Fitting %-22s", "4. Strict ..."))))
})
