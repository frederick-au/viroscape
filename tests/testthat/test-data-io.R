test_that("the cache stores and returns a value without re-evaluating", {
  dir <- tempfile(); dir.create(dir)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0
  key <- vs_cache_key("test", 1:3)

  a <- with_cache(key, { counter$n <- counter$n + 1; 42 }, cache_dir = dir, quiet = TRUE)
  b <- with_cache(key, { counter$n <- counter$n + 1; 42 }, cache_dir = dir, quiet = TRUE)

  expect_equal(a, 42)
  expect_equal(b, 42)
  expect_equal(counter$n, 1)

  cc <- with_cache(key, { counter$n <- counter$n + 1; 99 }, cache_dir = dir,
                   refresh = TRUE, quiet = TRUE)
  expect_equal(cc, 99)
  expect_equal(counter$n, 2)
})

test_that("cache keys are deterministic and input sensitive", {
  expect_equal(vs_cache_key("a", 1), vs_cache_key("a", 1))
  expect_false(identical(vs_cache_key("a", 1), vs_cache_key("a", 2)))
})

test_that("environmental data reads and joins on country and year", {
  env <- toy_environment()
  f <- tempfile(fileext = ".csv")
  utils::write.csv(env, f, row.names = FALSE)
  back <- read_environment(f)

  expect_equal(nrow(back), nrow(env))
  expect_true(all(c("country", "year") %in% names(back)))
  expect_type(back$year, "integer")
})

test_that("joining drops rows whose predictors are missing", {
  env <- toy_environment()
  env$agricultural_land_pct[1] <- NA
  draw <- env[rep(1:3, each = 4), ]
  d <- tibble::tibble(id = sprintf("s%02d", seq_len(nrow(draw))),
                      country = draw$country, year = draw$year,
                      distance = stats::runif(nrow(draw)))

  expect_warning(ds <- join_environment(d, env, by = c("country", "year")))
  expect_true(nrow(ds$data) < nrow(d))
  expect_false(anyNA(ds$data$agricultural_land_pct))
})

test_that("standardising leaves the response and year alone", {
  env <- toy_environment()
  env$distance <- stats::runif(nrow(env))
  out <- standardise(env)
  expect_equal(out$year, env$year)
  expect_equal(out$distance, env$distance)
  expect_equal(mean(out$agricultural_land_pct), 0, tolerance = 1e-8)
})

test_that("sampling summary counts sequences per country", {
  d <- data.frame(country = rep(c("A", "B"), c(5, 3)),
                  year = c(2004:2008, 2010:2012))
  s <- sampling_summary(d)
  expect_equal(s$by_country$sequences, c(5, 3))
  expect_equal(s$by_country$years_sampled, c(5, 3))
})
