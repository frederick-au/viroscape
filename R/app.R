#' Launch the viroscape explorer
#'
#' An interactive front-end over the whole pipeline: choose a data source,
#' filter countries, years and hosts, swap the divergence metric and
#' substitution model, and rerun forward selection with different guards to see
#' how sensitive the conclusions are to each choice.
#'
#' @param analysis Optional `vs_analysis` to open with. Defaults to the bundled
#'   example data.
#' @param launch.browser Open a browser window.
#' @param ... Passed to [shiny::runApp()].
#' @return Invisibly, the result of [shiny::runApp()].
#' @export
launch_app <- function(analysis = NULL, launch.browser = interactive(), ...) {
  need_pkg("shiny", "launch the explorer")
  need_pkg("bslib", "launch the explorer")
  need_pkg("DT", "launch the explorer")
  app_dir <- system.file("app", package = "viroscape")
  if (!nzchar(app_dir)) vs_abort("The Shiny app could not be found in the installed package.")
  if (!is.null(analysis)) {
    options(viroscape.initial_analysis = analysis)
    on.exit(options(viroscape.initial_analysis = NULL), add = TRUE)
  }
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}
