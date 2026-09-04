test_that("plot helpers return ggplot objects", {
  skip_if_not_installed("ggplot2")
  a <- example_analysis(verbose = FALSE)

  expect_s3_class(plot_sampling(a), "ggplot")
  expect_s3_class(plot_model_selection(a), "ggplot")
  expect_s3_class(plot_distances(a), "ggplot")
  expect_s3_class(plot_coefficients(a), "ggplot")
  expect_s3_class(plot_gvif(a), "ggplot")
  expect_s3_class(plot_selection_path(a$selection), "ggplot")
  expect_s3_class(plot_predictor(a, a$dataset$predictors[1]), "ggplot")
})
