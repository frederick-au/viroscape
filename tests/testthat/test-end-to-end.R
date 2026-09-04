test_that("the bundled example data is installed and readable", {
  s <- example_sequences()
  env <- example_environment()

  expect_s3_class(s, "vs_sequences")
  expect_gt(length(s$seqs), 100)
  expect_true(all(c("country", "year", "agricultural_land_pct",
                    "livestock_index") %in% names(env)))
  expect_false(anyNA(env$agricultural_land_pct))
  expect_gte(length(unique(env$country)), 10)
})

test_that("a subsampled end-to-end run produces a validated model", {
  s <- small_example(60)
  aln <- align_sequences(s, method = "none")
  msel <- select_substitution_model(aln, models = c("JTT", "WAG"), criterion = "BIC")
  d <- distance_from_consensus(aln, metric = "ml_aa", model = msel)
  ds <- join_environment(d, example_environment())
  sel <- suppressWarnings(forward_select(ds, base_terms = c("country", "year"),
                                         calibrate = FALSE, verbose = FALSE))

  expect_s3_class(sel, "vs_selection")
  expect_true(nrow(sel$steps) > 0)
  expect_true(all(c("round", "predictor", "accepted", "reason") %in% names(sel$steps)))

  fit <- model_fit_summary(sel$model)
  expect_true(fit$adj_r_squared > 0)
})

test_that("the full example analysis recovers the predictors planted in the data", {
  skip_on_cran()
  a <- example_analysis(verbose = FALSE)

  expect_s3_class(a, "vs_analysis")
  expect_true("agricultural_land_pct" %in% a$selection$selected)
  expect_true("livestock_index" %in% a$selection$selected)

  co <- a$coefficients
  expect_gt(co$estimate[co$term == "agricultural_land_pct"], 0)
  expect_lt(co$estimate[co$term == "livestock_index"], 0)
  expect_true(all(a$lrt$below_alpha))
  expect_gt(a$fit$adj_r_squared, 0.5)
})

test_that("a report can be generated from an analysis", {
  a <- small_analysis()
  f <- tempfile(fileext = ".txt")
  report_analysis(a, file = f)
  txt <- readLines(f)

  expect_true(file.exists(f))
  expect_true(any(grepl("Forward selection steps", txt)))
  expect_true(any(grepl("Adjusted R-squared", txt)))
})

test_that("a report built on the simulated example says so, and says it loudly", {
  a <- small_analysis()
  txt <- report_analysis(a)

  expect_true(any(grepl("SIMULATED DATA", txt)))
  expect_true(any(grepl("planted effects", txt)))
  expect_true(any(grepl("was NOT calibrated", txt)))
  expect_true(any(grepl("per standard deviation", txt)))
})
