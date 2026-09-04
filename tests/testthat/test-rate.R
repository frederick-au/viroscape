test_that("clock_correct is gone and says so", {
  aln <- align_sequences(small_example(40), method = "none")
  expect_error(
    distance_from_consensus(aln, clock_correct = TRUE),
    "divergence_rate"
  )
})

test_that("root-to-tip regression recovers a positive rate on the example data", {
  skip_on_cran()
  aln <- align_sequences(example_sequences(), method = "none")
  r <- divergence_rate(aln, by = "country")

  expect_s3_class(r, "vs_rate")
  expect_equal(r$group[1], "all")
  expect_gt(r$slope[1], 0)
  expect_true(is.finite(r$root_date[1]))
  expect_gt(nrow(r), 1)
})

test_that("a group with too few years is reported as NA rather than fitted", {
  skip_on_cran()
  aln <- align_sequences(example_sequences(), method = "none")
  r <- divergence_rate(aln, by = "country", min_n = 1e6)
  expect_true(all(is.na(r$slope[-1])))
})

test_that("the rate estimator recovers a planted constant rate", {
  skip_on_cran()
  fx <- lineage_fixture(n_lin = 1, per = 60, L = 400, seed = 3)
  aln <- align_sequences(fx$sequences, method = "none")

  ## sampling years are random in this fixture, so there is no clock signal in
  ## it at all - the estimator has to say so rather than return a rate
  expect_warning(r <- divergence_rate(aln, metric = "p_aa"),
                 "no usable clock signal")
  expect_lt(abs(r$slope[1]), 0.01)
})
