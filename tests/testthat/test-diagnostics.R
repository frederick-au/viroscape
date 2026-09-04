test_that("GVIF is scaled by degrees of freedom and matches car where available", {
  set.seed(5)
  d <- data.frame(
    y = stats::rnorm(80),
    g = factor(sample(letters[1:4], 80, TRUE)),
    x1 = stats::rnorm(80)
  )
  d$x2 <- d$x1 * 0.9 + stats::rnorm(80, 0, 0.1)
  m <- stats::lm(y ~ g + x1 + x2, data = d)

  gv <- gvif_table(m)
  expect_true(all(c("term", "gvif", "df", "gvif_scaled") %in% names(gv)))
  expect_setequal(gv$term, c("g", "x1", "x2"))
  expect_equal(gv$df[gv$term == "g"], 3)
  expect_gt(gv$gvif_scaled[gv$term == "x1"], 2)   # x1 and x2 are near-duplicates

  manual <- viroscape:::gvif_manual(m)
  expect_equal(manual$gvif_scaled[order(manual$term)],
               gv$gvif_scaled[order(gv$term)], tolerance = 1e-6)
})

test_that("a sign flip on a near-zero coefficient is not treated as instability", {
  set.seed(6)
  n <- 200
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$noise <- stats::rnorm(n, 0, 0.001)            # essentially no effect
  d$y <- 2 * d$x1 + stats::rnorm(n, 0, 0.5)

  m0 <- stats::lm(y ~ x1 + noise, data = d)
  m1 <- stats::lm(y ~ x1 + noise + x2, data = d)
  st <- coef_stability(m0, m1)

  expect_true(all(c("sign_flip", "lost_significance", "destabilising") %in% names(st)))
  expect_false(isTRUE(attr(st, "any_destabilising")))
})

test_that("grouping contrasts are exempt from the proportional shift rule", {
  set.seed(8)
  n <- 150
  g <- factor(sample(letters[1:5], n, TRUE))
  level_mean <- stats::setNames(stats::rnorm(5, 0, 2), letters[1:5])
  z <- level_mean[as.character(g)] + stats::rnorm(n, 0, 0.05)
  d <- data.frame(g = g, z = as.numeric(z))
  d$y <- 0.5 * d$z + stats::rnorm(n, 0, 0.4)

  m0 <- stats::lm(y ~ g, data = d)
  m1 <- stats::lm(y ~ g + z, data = d)

  strict <- coef_stability(m0, m1)
  exempt <- coef_stability(m0, m1, group_terms = "g")
  expect_true(all(exempt$is_group))
  expect_false(isTRUE(attr(exempt, "any_destabilising")))
  expect_true(sum(strict$substantial) >= sum(exempt$substantial & !exempt$is_group))
})
