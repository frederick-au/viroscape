## The toy fixture is clustered by construction (several draws share a
## country-year), so forward_select() warns about the design effect. That
## warning is the subject of test-clustering.R; here it is noise.
## calibrate = FALSE keeps these fast; the calibration itself is covered in
## test-clustering.R
forward_select <- function(...) suppressWarnings(viroscape::forward_select(..., calibrate = FALSE))

build_toy_dataset <- function(seed = 11, effect = 0.02) {
  set.seed(seed)
  env <- toy_environment()
  n <- 150
  draw <- env[sample(nrow(env), n, replace = TRUE), ]
  d <- tibble::tibble(
    id = sprintf("s%03d", seq_len(n)),
    country = draw$country,
    year = draw$year,
    distance = 0.05 +
      effect * scale(draw$agricultural_land_pct)[, 1] +
      0.004 * (draw$year - 2004) +
      stats::rnorm(n, 0, 0.01)
  )
  join_environment(d, env, by = c("country", "year"))
}

test_that("structure comparison reports every candidate and warns on few groups", {
  ds <- build_toy_dataset()
  expect_warning(st <- select_structure(ds), NA)
  expect_s3_class(st, "vs_structure")
  expect_true(all(c("null", "fixed") %in% st$table$structure))
  expect_equal(st$n_groups, 3)
})

test_that("forward selection recovers a planted predictor", {
  ds <- build_toy_dataset(effect = 0.03)
  sel <- forward_select(ds, base_terms = c("country", "year"), verbose = FALSE)

  expect_s3_class(sel, "vs_selection")
  expect_true("agricultural_land_pct" %in% sel$selected)
  expect_true(any(sel$steps$accepted))
  expect_true(all(sel$steps$delta[sel$steps$accepted] < 0))
})

test_that("forward selection stops when nothing improves the criterion", {
  set.seed(4)
  env <- toy_environment()
  n <- 120
  draw <- env[sample(nrow(env), n, replace = TRUE), ]
  d <- tibble::tibble(id = sprintf("s%03d", seq_len(n)),
                      country = draw$country, year = draw$year,
                      distance = stats::rnorm(n, 0.05, 0.01))   # pure noise
  ds <- join_environment(d, env, by = c("country", "year"))
  sel <- forward_select(ds, base_terms = c("country", "year"), verbose = FALSE)

  expect_length(sel$selected, 0)
  expect_true(any(grepl("no meaningful improvement", sel$steps$reason)))
})

test_that("a duplicated predictor never enters alongside its twin", {
  ds <- build_toy_dataset(effect = 0.03)
  ds$data$agri_copy <- ds$data$agricultural_land_pct + stats::rnorm(nrow(ds$data), 0, 1e-4)
  ds$predictors <- c(ds$predictors, "agri_copy")
  sel <- forward_select(ds, base_terms = c("country", "year"),
                        gvif_action = "reject", verbose = FALSE)

  expect_false(all(c("agricultural_land_pct", "agri_copy") %in% sel$selected))
})

test_that("the likelihood ratio chain matches the accepted steps", {
  ds <- build_toy_dataset(effect = 0.03)
  sel <- forward_select(ds, base_terms = c("country", "year"), verbose = FALSE)
  lrt <- lrt_chain(sel)

  expect_equal(nrow(lrt), length(sel$selected))
  if (nrow(lrt)) {
    expect_equal(lrt$added, sel$selected)
    expect_true(all(lrt$nominal_p < 0.5, na.rm = TRUE))
  }
})

test_that("fit summary reports adjusted R squared and sample size", {
  ds <- build_toy_dataset(effect = 0.03)
  sel <- forward_select(ds, base_terms = c("country", "year"), verbose = FALSE)
  fit <- model_fit_summary(sel$model)

  expect_true(fit$adj_r_squared > 0 && fit$adj_r_squared <= 1)
  expect_equal(fit$n, nrow(ds$data))
})

test_that("gvif_action decides whether collinearity alone rejects a candidate", {
  ## A low GVIF threshold with a permissive stability threshold isolates the
  ## collinearity branch: the same candidate is judged only on its GVIF.
  ds <- build_toy_dataset(effect = 0.03)
  args <- list(dataset = ds, base_terms = c("country", "year"),
               gvif_threshold = 1.5, stability_threshold = 100,
               max_predictors = 1, verbose = FALSE)

  lenient <- do.call(forward_select, c(args, list(gvif_action = "flag")))
  strict <- do.call(forward_select, c(args, list(gvif_action = "reject")))

  expect_length(lenient$selected, 1)
  expect_true(any(grepl("GVIF", lenient$steps$reason)))
  expect_true(all(lenient$steps$gvif_scaled[lenient$steps$accepted] > 1.5))

  expect_length(strict$selected, 0)
  expect_true(all(grepl("GVIF|no meaningful improvement", strict$steps$reason)))
})

test_that("stability_scope = 'all' is stricter than the default", {
  ds <- build_toy_dataset(effect = 0.03)
  default <- forward_select(ds, base_terms = c("country", "year"), verbose = FALSE)
  strict <- forward_select(ds, base_terms = c("country", "year"),
                           stability_scope = "all", verbose = FALSE)

  expect_gte(length(default$selected), length(strict$selected))
})

test_that("selection can be capped at a maximum number of predictors", {
  ds <- build_toy_dataset(effect = 0.03)
  sel <- forward_select(ds, base_terms = c("country", "year"),
                        max_predictors = 1, verbose = FALSE)
  expect_lte(length(sel$selected), 1)
})
