#' Join divergence estimates to environmental predictors
#'
#' @param distances A `vs_distances` tibble from [distance_from_consensus()].
#' @param environment A country-year (or country-admin1-year) tibble.
#' @param by Join keys.
#' @param scale Standardise numeric predictors so coefficients are comparable.
#'   Coefficients are then per standard deviation, not per natural unit; the
#'   centre and scale are kept and reported by [tidy_coefficients()], and
#'   [unscale()] converts back.
#' @param drop_incomplete Drop rows with missing predictor values, which model
#'   comparison by AIC requires: models fitted to different numbers of rows are
#'   not comparable.
#' @param cluster Columns defining the level at which the predictors actually
#'   vary. Defaults to the join keys, which is almost always right: an
#'   environmental indicator is a country-year quantity, and every sequence from
#'   that country-year carries the same value. Sequence-level inference ignores
#'   this and is badly anti-conservative; see [cluster_summary()] and
#'   [aggregate_to_cluster()].
#' @return A `vs_dataset` object.
#' @export
join_environment <- function(distances, environment,
                             by = c("country", "year"),
                             scale = TRUE, drop_incomplete = TRUE,
                             cluster = by) {
  missing_keys <- setdiff(by, intersect(names(distances), names(environment)))
  if (length(missing_keys)) {
    vs_abort("Join key{?s} {.field {missing_keys}} missing from one of the inputs.")
  }
  distances <- tibble::as_tibble(distances)
  ## the environmental table is authoritative for any column it shares with the
  ## divergence table, so drop the duplicates rather than creating .x/.y pairs
  overlap <- setdiff(intersect(names(distances), names(environment)), by)
  if (length(overlap)) {
    vs_inform("Overriding {length(overlap)} column{?s} already present in the divergence table: {.field {overlap}}.")
    distances <- distances[setdiff(names(distances), overlap)]
  }
  d <- dplyr::left_join(distances, environment, by = by)

  predictors <- setdiff(names(environment), c(by, "iso3", "country_name"))
  predictors <- predictors[vapply(d[predictors], is.numeric, logical(1))]

  n_before <- nrow(d)
  if (drop_incomplete) {
    d <- d[stats::complete.cases(d[c("distance", predictors, by)]), , drop = FALSE]
  }
  dropped <- n_before - nrow(d)
  if (dropped > 0) {
    vs_warn("Dropped {dropped} row{?s} with missing predictor or key values.")
  }
  if (!nrow(d)) vs_abort("No rows remain after joining; check that country names and years match.")

  d$country <- factor(d$country)
  if (scale) d <- standardise(d, cols = predictors)

  cluster <- intersect(cluster %||% character(), names(d))
  if (length(cluster)) {
    d$.cluster <- interaction(d[cluster], drop = TRUE, sep = " / ")
  }

  out <- structure(
    list(data = d, predictors = predictors, response = "distance",
         keys = by, scaled = scale, cluster = cluster, weights = NULL,
         level = "sequence",
         metric = attr(distances, "metric"),
         model = attr(distances, "model")),
    class = "vs_dataset"
  )
  cs <- cluster_summary(out)
  out$icc <- cs$icc
  out$design_effect <- cs$design_effect
  out
}

#' How clustered is a dataset?
#'
#' Environmental predictors vary at the country-year level, while the response
#' is measured per sequence, and sequences from the same country-year are near
#' identical by descent. The intraclass correlation says how much of the
#' response variance is between clusters rather than within them; the design
#' effect says how many effectively independent observations the dataset
#' contains relative to its row count. A design effect of 5 means 500 sequences
#' carry about as much information as 100 independent ones - and every
#' criterion computed per sequence (AIC, GVIF, likelihood ratio, `summary()`
#' standard errors) is wrong by roughly that factor.
#'
#' @param dataset A `vs_dataset`.
#' @return A list with `n`, `n_clusters`, `mean_size`, `icc` and
#'   `design_effect`.
#' @export
cluster_summary <- function(dataset) {
  stopifnot(inherits(dataset, "vs_dataset"))
  d <- dataset$data
  empty <- list(n = nrow(d), n_clusters = NA_integer_, mean_size = NA_real_,
                icc = NA_real_, design_effect = NA_real_)
  if (!length(dataset$cluster) || !".cluster" %in% names(d)) return(empty)
  g <- droplevels(as.factor(d$.cluster))
  y <- d[[dataset$response]]
  k <- nlevels(g)
  n <- length(y)
  if (k < 2L || n <= k) return(utils::modifyList(empty, list(n_clusters = k)))
  ni <- as.numeric(table(g))
  mu <- tapply(y, g, mean)
  msb <- sum(ni * (mu - mean(y))^2) / (k - 1)
  msw <- sum((y - stats::ave(y, g))^2) / (n - k)
  n0 <- (n - sum(ni^2) / n) / (k - 1)
  s2b <- max((msb - msw) / n0, 0)
  icc <- if ((s2b + msw) > 0) s2b / (s2b + msw) else NA_real_
  list(n = n, n_clusters = k, mean_size = n / k, icc = icc,
       design_effect = 1 + (n / k - 1) * icc)
}

#' Aggregate a dataset to the level at which the predictors vary
#'
#' Collapses to one row per cluster, averaging the response. This is the level
#' the environmental predictors live at, so it is the level at which AIC,
#' standard errors and likelihood ratios mean what they say. The number of
#' sequences behind each row is retained as a regression weight.
#'
#' This is a necessary correction, not a sufficient one. Forward selection over
#' several candidates is still a search, and still needs
#' [calibrate_selection()] before any delta is interpreted.
#'
#' @param dataset A `vs_dataset`.
#' @param weight Weight the aggregated rows by the number of sequences behind
#'   them.
#' @return A `vs_dataset` with one row per cluster.
#' @export
aggregate_to_cluster <- function(dataset, weight = TRUE) {
  stopifnot(inherits(dataset, "vs_dataset"))
  if (!length(dataset$cluster)) vs_abort("This dataset has no cluster definition.")
  d <- dataset$data
  keys <- dataset$cluster
  num <- unique(c(dataset$response, dataset$predictors,
                  intersect("year", names(d))))
  num <- num[vapply(d[num], is.numeric, logical(1))]

  split_key <- interaction(d[keys], drop = TRUE, sep = "\r")
  agg <- lapply(split(d, split_key), function(part) {
    row <- part[1, keys, drop = FALSE]
    for (nm in num) row[[nm]] <- mean(part[[nm]], na.rm = TRUE)
    row$.n <- nrow(part)
    row
  })
  agg <- dplyr::bind_rows(agg)
  if ("country" %in% names(agg)) agg$country <- factor(as.character(agg$country))
  agg$.cluster <- interaction(agg[keys], drop = TRUE, sep = " / ")
  rownames(agg) <- NULL

  out <- dataset
  agg <- tibble::as_tibble(agg)
  ## the standardisation parameters describe the predictors, which aggregation
  ## does not change, so they must survive it or unscale() has nothing to work with
  attr(agg, "vs_scaling") <- attr(dataset$data, "vs_scaling")
  out$data <- agg
  out$weights <- if (weight) ".n" else NULL
  out$level <- "cluster"
  cs <- cluster_summary(out)
  out$icc <- cs$icc
  out$design_effect <- cs$design_effect
  out
}

#' @export
print.vs_dataset <- function(x, ...) {
  cli::cli_h3("vs_dataset")
  cli::cli_text("{nrow(x$data)} observation{?s}, {length(x$predictors)} candidate predictor{?s}")
  cli::cli_text("Response: {.field {x$response}} (metric {x$metric})")
  n_ctry <- nlevels(x$data$country)
  cli::cli_text("Groups: {n_ctry} countr{?y/ies}, years {min(x$data$year)}-{max(x$data$year)}")
  if (n_ctry < 10) {
    cli::cli_alert_warning(
      "Only {n_ctry} group{?s}: variance components for a country random effect will be unreliable (rule of thumb: 10+)."
    )
  }
  cli::cli_text("Predictors: {paste(x$predictors, collapse = ', ')}")
  if (isTRUE(x$scaled)) {
    cli::cli_text("Predictors are standardised: coefficients are per standard deviation, not per natural unit.")
  }
  if (identical(x$level, "sequence") && is.finite(x$design_effect %||% NA_real_) &&
      x$design_effect > 1.5) {
    cli::cli_alert_warning(c(
      "Clustered by {.field {paste(x$cluster, collapse = ' x ')}}: ICC {round(x$icc, 3)}, design effect {round(x$design_effect, 1)}."
    ))
    cli::cli_text("{.emph Sequence-level AIC, standard errors and likelihood ratios are anti-conservative by roughly that factor. See {.fn aggregate_to_cluster} and {.fn calibrate_selection}.}")
  }
  invisible(x)
}

#' Summarise sampling effort by country and year
#'
#' Uneven sampling is the single biggest threat to interpreting country
#' coefficients, so this is worth looking at before any modelling.
#'
#' @param x A `vs_dataset`, `vs_distances` or data frame with `country`/`year`.
#' @return A tibble of counts by country and year, with per-country totals and
#'   the number of distinct years sampled.
#' @export
sampling_summary <- function(x) {
  d <- if (inherits(x, "vs_dataset")) x$data else tibble::as_tibble(x)
  by_cy <- dplyr::count(d, country, year, name = "n")
  totals <- dplyr::summarise(
    dplyr::group_by(by_cy, country),
    sequences = sum(n),
    years_sampled = dplyr::n_distinct(year),
    first_year = min(year), last_year = max(year),
    .groups = "drop"
  )
  attr(by_cy, "totals") <- totals
  list(by_country_year = by_cy, by_country = totals)
}
