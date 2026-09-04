vs_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.1),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      axis.title = ggplot2::element_text(colour = "grey25"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' Sampling effort by country and year
#'
#' Reproduces the bubble plot used to show how unevenly sequences are deposited
#' across the region, which is the context needed to read country coefficients
#' honestly.
#'
#' @param x A `vs_dataset`, `vs_analysis` or data frame with `country`/`year`.
#' @return A ggplot object.
#' @export
plot_sampling <- function(x) {
  d <- extract_data(x)
  counts <- dplyr::count(d, country, year, name = "n")
  ggplot2::ggplot(counts, ggplot2::aes(x = year, y = country, size = n)) +
    ggplot2::geom_point(colour = "#2c6e8f", alpha = 0.75) +
    ggplot2::scale_size_area(max_size = 12, name = "Sequences") +
    ggplot2::labs(
      title = "Sequence sampling across country and year",
      subtitle = "Bubble area is the number of sequences deposited",
      x = NULL, y = NULL
    ) +
    vs_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "grey92"))
}

#' Information criterion path across forward selection
#'
#' @param selection A `vs_selection` or `vs_analysis`.
#' @return A ggplot object.
#' @export
plot_selection_path <- function(selection) {
  if (inherits(selection, "vs_analysis")) selection <- selection$selection
  steps <- selection$steps
  acc <- steps[which(steps$accepted), , drop = FALSE]
  crit <- tolower(selection$criterion)
  base_val <- stats::AIC(selection$chain[[1]])
  if (selection$criterion == "BIC") base_val <- stats::BIC(selection$chain[[1]])
  path <- tibble::tibble(
    step = seq_len(nrow(acc) + 1L) - 1L,
    predictor = c("base model", acc$predictor),
    value = c(base_val, acc[[crit]])
  )
  path$label <- c("", sprintf("%+.1f", diff(path$value)))
  ggplot2::ggplot(path, ggplot2::aes(x = stats::reorder(predictor, step), y = value, group = 1)) +
    ggplot2::geom_line(colour = "#2c6e8f", linewidth = 0.8) +
    ggplot2::geom_point(size = 3, colour = "#2c6e8f") +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -1, size = 3.4, colour = "grey30") +
    ggplot2::labs(
      title = sprintf("%s at each accepted step", selection$criterion),
      subtitle = "Labels show the change relative to the previous step",
      x = NULL, y = selection$criterion
    ) +
    ggplot2::expand_limits(y = max(path$value) + 0.05 * diff(range(path$value))) +
    vs_theme()
}

#' Collinearity diagnostics
#'
#' @param x A fitted model, `vs_selection` or `vs_analysis`.
#' @param threshold Threshold line to draw.
#' @return A ggplot object.
#' @export
plot_gvif <- function(x, threshold = 10) {
  gv <- if (inherits(x, "vs_analysis")) x$gvif
        else if (inherits(x, "vs_selection")) gvif_table(x$model)
        else gvif_table(x)
  if (!nrow(gv)) return(empty_plot("No collinearity diagnostics available"))
  gv$flag <- gv$gvif_scaled > threshold
  ggplot2::ggplot(gv, ggplot2::aes(x = stats::reorder(term, gvif_scaled),
                                   y = gvif_scaled, fill = flag)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "#2c6e8f", `TRUE` = "#c1553b"),
                               guide = "none") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Generalised variance inflation",
      subtitle = sprintf("Scaled GVIF; dashed line at the %g threshold", threshold),
      x = NULL, y = expression(GVIF^{1/(2 * df)})
    ) +
    vs_theme()
}

#' Coefficient estimates with confidence intervals
#'
#' @param x A fitted model, `vs_selection` or `vs_analysis`.
#' @param drop_intercept Hide the intercept, which is usually on a different scale.
#' @param level Confidence level.
#' @return A ggplot object.
#' @export
plot_coefficients <- function(x, drop_intercept = TRUE, level = 0.95) {
  co <- if (inherits(x, "vs_analysis")) x$coefficients else tidy_coefficients(x, level)
  if (drop_intercept) co <- co[co$term != "(Intercept)", , drop = FALSE]
  if (!nrow(co)) return(empty_plot("No coefficients to display"))
  co$sig <- !is.na(co$p_value) & co$p_value < 0.05
  ggplot2::ggplot(co, ggplot2::aes(x = estimate, y = stats::reorder(term, estimate))) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = conf.low, xmax = conf.high),
                            colour = "grey35", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = sig), size = 2.8) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#2c6e8f"),
                                 guide = "none") +
    ggplot2::labs(
      title = "Final model coefficients",
      subtitle = sprintf("%.0f%% confidence intervals; positive values indicate greater divergence from consensus", level * 100),
      x = "Estimate", y = NULL
    ) +
    vs_theme()
}

#' Distribution of divergence from consensus
#'
#' @param x A `vs_dataset`, `vs_distances` or `vs_analysis`.
#' @param by Optional grouping column, e.g. `"country"` or `"host_class"`.
#' @return A ggplot object.
#' @export
plot_distances <- function(x, by = "country") {
  d <- extract_data(x)
  if (!is.null(by) && by %in% names(d)) {
    ggplot2::ggplot(d, ggplot2::aes(x = distance, y = .data[[by]])) +
      ggplot2::geom_boxplot(outlier.alpha = 0.35, fill = "#2c6e8f", alpha = 0.25,
                            colour = "#20506a", width = 0.55) +
      ggplot2::labs(title = "Divergence from regional consensus",
                    x = "Genetic distance", y = NULL) +
      vs_theme()
  } else {
    ggplot2::ggplot(d, ggplot2::aes(x = distance)) +
      ggplot2::geom_histogram(bins = 30, fill = "#2c6e8f", colour = "white") +
      ggplot2::labs(title = "Divergence from regional consensus",
                    x = "Genetic distance", y = "Sequences") +
      vs_theme()
  }
}

#' Relationship between a predictor and divergence
#'
#' @param x A `vs_dataset` or `vs_analysis`.
#' @param predictor Column name of the predictor.
#' @param colour_by Optional grouping column for point colour.
#' @return A ggplot object.
#' @export
plot_predictor <- function(x, predictor, colour_by = "country") {
  d <- extract_data(x)
  if (!predictor %in% names(d)) vs_abort("Predictor {.field {predictor}} not found.")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[predictor]], y = distance))
  if (!is.null(colour_by) && colour_by %in% names(d)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data[[colour_by]]),
                                 alpha = 0.7, size = 1.9)
  } else {
    p <- p + ggplot2::geom_point(alpha = 0.7, size = 1.9, colour = "#2c6e8f")
  }
  p +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         colour = "grey20", linewidth = 0.7) +
    ggplot2::labs(title = sprintf("Divergence against %s", predictor),
                  x = predictor, y = "Genetic distance", colour = NULL) +
    vs_theme()
}

#' Substitution model ranking
#'
#' @param x A `vs_model_selection` or `vs_analysis`.
#' @return A ggplot object.
#' @export
plot_model_selection <- function(x) {
  if (inherits(x, "vs_analysis")) x <- x$model_selection
  tab <- x$table
  crit <- x$criterion
  tab$best <- tab$model == x$best
  ggplot2::ggplot(tab, ggplot2::aes(x = stats::reorder(model, -.data[[crit]]),
                                    y = .data[[crit]], fill = best)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::coord_flip(ylim = c(min(tab[[crit]]) * 0.999, max(tab[[crit]]) * 1.001)) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#2c6e8f"),
                               guide = "none") +
    ggplot2::labs(title = sprintf("Substitution model ranking by %s", crit),
                  subtitle = sprintf("Selected: %s", x$best),
                  x = NULL, y = crit) +
    vs_theme()
}

extract_data <- function(x) {
  if (inherits(x, "vs_analysis")) return(x$dataset$data)
  if (inherits(x, "vs_dataset")) return(x$data)
  tibble::as_tibble(x)
}

empty_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg, colour = "grey45") +
    ggplot2::theme_void()
}
