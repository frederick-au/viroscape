test_that("the Shiny app object and its UI construct without error", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  app_dir <- system.file("app", package = "viroscape")
  expect_true(nzchar(app_dir))
  expect_true(file.exists(file.path(app_dir, "app.R")))

  env <- new.env(parent = globalenv())
  obj <- source(file.path(app_dir, "app.R"), local = env)$value
  expect_s3_class(obj, "shiny.appobj")

  ## rendering the UI catches malformed tags and bad bslib arguments
  html <- as.character(htmltools::renderTags(env$ui)$html)
  expect_true(nchar(html) > 1000)
  expect_true(grepl("Forward selection|Selection", html))
})
