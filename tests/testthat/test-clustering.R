test_that("clustering is detected and quantified", {
  ds <- noise_dataset(seed = 11)
  cs <- cluster_summary(ds)

  expect_equal(cs$n_clusters, 50)
  expect_equal(cs$mean_size, 20)
  expect_gt(cs$icc, 0.8)
  expect_gt(cs$design_effect, 10)
})

test_that("aggregating collapses to one row per country-year and keeps weights", {
  ds <- noise_dataset(seed = 12)
  ag <- aggregate_to_cluster(ds)

  expect_equal(nrow(ag$data), 50)
  expect_equal(sum(ag$data$.n), 1000)
  expect_equal(ag$level, "cluster")
  expect_equal(ag$weights, ".n")
  ## the predictors are constant within a cluster, so averaging must not move them
  expect_equal(sort(unique(round(ag$data$x1, 6))),
               sort(unique(round(ds$data$x1, 6))))
})

test_that("sequence-level selection on clustered data is loudly flagged", {
  ds <- noise_dataset(seed = 13)
  expect_warning(
    forward_select(ds, base_terms = c("country", "year"), calibrate = FALSE,
                   verbose = FALSE),
    "design effect"
  )
})

## The regression test for the defect that mattered most: with no effect in the
## data at all, the procedure must not be allowed to report a discovery.
test_that("pure noise does not survive calibration", {
  ds <- noise_dataset(seed = 14)
  ag <- aggregate_to_cluster(ds)
  cal <- calibrate_selection(ag, base_terms = c("country", "year"),
                             B = 99, seed = 99)

  expect_s3_class(cal, "vs_calibration")
  ## a nominal -2 is reached by chance far more than 5% of the time here, which
  ## is exactly why the delta threshold cannot be read as a test
  expect_lt(cal$threshold_05, -2)
  expect_gt(cal$p_value, 0.05)
})

test_that("a real planted effect does survive calibration", {
  skip_on_cran()
  ds <- aggregate_to_cluster(suppressWarnings(example_analysis(level = "sequence"))$dataset)
  cal <- calibrate_selection(ds, base_terms = c("country", "year"),
                             B = 99, seed = 7)
  expect_lte(cal$p_value, 0.05)
  expect_lt(cal$observed, cal$threshold_05)
})

test_that("cluster-robust standard errors are larger than model-based ones", {
  skip_if_not_installed("sandwich")
  ds <- noise_dataset(seed = 15)
  m <- stats::lm(distance ~ country + year + x1, data = ds$data)

  classical <- tidy_coefficients(m, vcov = "classical", dataset = ds)
  robust <- tidy_coefficients(m, vcov = "cluster", dataset = ds)

  expect_equal(classical$estimate, robust$estimate)
  expect_gt(mean(robust$std_error / classical$std_error), 1.5)
  expect_match(attr(robust, "vcov_type"), "cluster-robust")
})

test_that("standardised coefficients report their unit and can be unscaled", {
  ds <- small_dataset()
  m <- stats::lm(distance ~ country + year + agricultural_land_pct, data = ds$data)
  co <- tidy_coefficients(m, dataset = ds)
  i <- which(co$term == "agricultural_land_pct")

  ## the predictor is standardised, so the coefficient is per SD and the table
  ## has to say so - this is the units ambiguity the analysis kept tripping on
  expect_equal(co$unit[i], "standard deviation")
  expect_false(is.na(co$sd[i]))
  expect_equal(co$unit[co$term == "year"], "natural unit")

  nat <- unscale(co, ds)
  expect_equal(nat$estimate[i], co$estimate[i] / co$sd[i])
  expect_equal(nat$conf.low[i], co$conf.low[i] / co$sd[i])
  expect_true(all(nat$unit == "natural unit"))
})

test_that("calibration runs inside forward_select and tightens the threshold", {
  ds <- aggregate_to_cluster(noise_dataset(seed = 21))
  sel <- forward_select(ds, base_terms = c("country", "year"),
                        calibrate = 99, verbose = FALSE)

  expect_s3_class(sel$calibration, "vs_calibration")
  ## the honest threshold on data shaped like this is far stricter than -2, and
  ## it has to be applied during the search, not offered afterwards
  expect_lt(sel$thresholds$delta, -2)
  expect_equal(sel$thresholds$delta_nominal, -2)
  expect_equal(sel$thresholds$delta, sel$calibration$threshold_05)
})

test_that("calibrate = FALSE leaves the nominal threshold alone", {
  ds <- aggregate_to_cluster(noise_dataset(seed = 22))
  sel <- forward_select(ds, base_terms = c("country", "year"),
                        calibrate = FALSE, verbose = FALSE)
  expect_null(sel$calibration)
  expect_equal(sel$thresholds$delta, -2)
})

## Finding 2 of the 0.3.0 review: the stability guard vetoed a real effect
## because adding it moved the year coefficient. Base terms are controls, not
## estimands, so movement in them cannot be evidence against a candidate.
test_that("a real effect is not vetoed by movement in a forced base term", {
  set.seed(4)
  cl <- expand.grid(country = LETTERS[1:5], year = 2003:2012)
  P <- matrix(stats::rnorm(nrow(cl) * 5), nrow(cl))
  colnames(P) <- paste0("pred", 1:5)
  cl <- cbind(cl, P)
  cl$u <- stats::rnorm(nrow(cl), 0, 0.01)
  d <- cl[rep(seq_len(nrow(cl)), each = 20), ]
  ## a genuine effect, and deliberately no year effect at all, so that year's
  ## near-zero coefficient moves a long way proportionally when pred1 enters
  d$distance <- 0.02 + 0.02 * d$pred1 + d$u + stats::rnorm(nrow(d), 0, 0.003)
  d$country <- factor(d$country); d$u <- NULL; rownames(d) <- NULL
  d$.cluster <- interaction(d[c("country", "year")], drop = TRUE, sep = " / ")
  ds <- structure(list(data = tibble::as_tibble(d), predictors = paste0("pred", 1:5),
                       response = "distance", keys = c("country", "year"),
                       scaled = FALSE, cluster = c("country", "year"),
                       weights = NULL, level = "sequence"), class = "vs_dataset")
  cs <- cluster_summary(ds); ds$icc <- cs$icc; ds$design_effect <- cs$design_effect
  ag <- aggregate_to_cluster(ds)

  sel <- forward_select(ag, base_terms = c("country", "year"), calibrate = 99,
                        verbose = FALSE)
  expect_true("pred1" %in% sel$selected)
  ## the movement is still reported, it just does not veto
  expect_true("base_max_shift" %in% names(sel$steps))
})

test_that("a selection that contradicts its own calibration says so", {
  set.seed(4)
  cl <- expand.grid(country = LETTERS[1:5], year = 2003:2012)
  P <- matrix(stats::rnorm(nrow(cl) * 5), nrow(cl))
  colnames(P) <- paste0("pred", 1:5)
  cl <- cbind(cl, P); cl$u <- stats::rnorm(nrow(cl), 0, 0.01)
  d <- cl[rep(seq_len(nrow(cl)), each = 20), ]
  d$distance <- 0.02 + 0.02 * d$pred1 + d$u + stats::rnorm(nrow(d), 0, 0.003)
  d$country <- factor(d$country); d$u <- NULL; rownames(d) <- NULL
  d$.cluster <- interaction(d[c("country", "year")], drop = TRUE, sep = " / ")
  ds <- structure(list(data = tibble::as_tibble(d), predictors = paste0("pred", 1:5),
                       response = "distance", keys = c("country", "year"),
                       scaled = FALSE, cluster = c("country", "year"),
                       weights = NULL, level = "sequence"), class = "vs_dataset")
  cs <- cluster_summary(ds); ds$icc <- cs$icc; ds$design_effect <- cs$design_effect

  ## the strict scope reinstates the veto; the object must not then report
  ## "nothing selected" while its calibration says the winner beat the null
  sel <- forward_select(aggregate_to_cluster(ds), base_terms = c("country", "year"),
                        calibrate = 99, stability_scope = "all", verbose = FALSE)
  if (!length(sel$selected) && isTRUE(sel$calibration$p_value < 0.05)) {
    expect_false(is.null(sel$conflict))
    expect_output(print(sel), "disagrees with its own calibration")
  } else {
    expect_true("pred1" %in% sel$selected)
  }
})

test_that("allow_unstable_base is deprecated in favour of an accurate name", {
  ds <- aggregate_to_cluster(noise_dataset(seed = 31))
  expect_warning(
    forward_select(ds, base_terms = c("country", "year"), calibrate = FALSE,
                   allow_unstable_base = TRUE, verbose = FALSE),
    "deprecated"
  )
})

test_that("weighted cluster refitting reproduces the sequence-level estimate", {
  set.seed(6)
  cl <- expand.grid(country = LETTERS[1:6], year = 2003:2012)
  cl$x <- stats::rnorm(nrow(cl))
  cl$u <- stats::rnorm(nrow(cl), 0, 0.01)
  n_i <- sample(5:30, nrow(cl), TRUE)          # uneven, as real sampling is
  d <- cl[rep(seq_len(nrow(cl)), times = n_i), ]
  d$distance <- 0.02 + 0.02 * d$x + d$u + stats::rnorm(nrow(d), 0, 0.003)
  d$country <- factor(d$country); d$u <- NULL; rownames(d) <- NULL
  d$.cluster <- interaction(d[c("country", "year")], drop = TRUE, sep = " / ")
  ds <- structure(list(data = tibble::as_tibble(d), predictors = "x",
                       response = "distance", keys = c("country", "year"),
                       scaled = FALSE, cluster = c("country", "year"),
                       weights = NULL, level = "sequence"), class = "vs_dataset")

  seq_fit <- stats::lm(distance ~ country + x, data = ds$data)
  wt <- aggregate_to_cluster(ds, weight = TRUE)
  unwt <- aggregate_to_cluster(ds, weight = FALSE)
  w_fit <- stats::lm(distance ~ country + x, data = wt$data, weights = wt$data$.n)
  u_fit <- stats::lm(distance ~ country + x, data = unwt$data)

  ## size-weighted cluster means reproduce the sequence-level coefficient
  ## exactly, so weighting is what decides whether the number moves at all
  expect_equal(unname(stats::coef(w_fit)["x"]),
               unname(stats::coef(seq_fit)["x"]), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(unname(stats::coef(u_fit)["x"]),
                                unname(stats::coef(seq_fit)["x"]),
                                tolerance = 1e-6)))
})

test_that("holdout_estimate selects and estimates on different clusters", {
  set.seed(8)
  cl <- expand.grid(country = LETTERS[1:8], year = 2003:2012)
  P <- matrix(stats::rnorm(nrow(cl) * 4), nrow(cl))
  colnames(P) <- paste0("pred", 1:4)
  cl <- cbind(cl, P); cl$u <- stats::rnorm(nrow(cl), 0, 0.008)
  d <- cl[rep(seq_len(nrow(cl)), each = 10), ]
  d$distance <- 0.02 + 0.015 * d$pred1 + d$u + stats::rnorm(nrow(d), 0, 0.003)
  d$country <- factor(d$country); d$u <- NULL; rownames(d) <- NULL
  d$.cluster <- interaction(d[c("country", "year")], drop = TRUE, sep = " / ")
  ds <- structure(list(data = tibble::as_tibble(d), predictors = paste0("pred", 1:4),
                       response = "distance", keys = c("country", "year"),
                       scaled = FALSE, cluster = c("country", "year"),
                       weights = NULL, level = "sequence"), class = "vs_dataset")

  ho <- holdout_estimate(aggregate_to_cluster(ds), base_terms = c("country", "year"),
                         repeats = 8, seed = 2)
  expect_true("pred1" %in% ho$term)
  row <- ho[ho$term == "pred1", ]
  expect_true(is.finite(row$estimate_heldout))
  ## the held-out estimate is a different number from the one the search saw
  expect_false(isTRUE(all.equal(row$estimate_heldout, row$estimate_selection_half)))
  ## and it carries its units, as the other coefficient tables do
  expect_true(all(c("unit", "sd") %in% names(ho)))
  expect_equal(row$unit, "natural unit")
  expect_true(all(c("times_selected", "n_splits") %in% names(ho)))
  expect_equal(row$n_splits, 8L)
  expect_lte(row$times_selected, row$n_splits)
})
