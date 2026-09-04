#' Calibrate a forward selection against its own null
#'
#' Forward selection is a search, and a search finds something even when there
#' is nothing to find. On data with five countries, ten years, twenty sequences
#' per country-year and seven pure-noise predictors, [forward_select()] retains
#' at least one predictor in 100% of runs and a mean of four and a half of the
#' seven. Aggregating to the cluster level with [aggregate_to_cluster()] brings
#' that down to roughly a third of runs - better, and still nowhere near the 5%
#' the delta threshold implies, because seven candidates tested at
#' `delta_threshold = -2` is about seven tests at alpha 0.046.
#'
#' There is no analytic fix for this. The empirical one is to permute the
#' predictors against the response and see how large a first-round improvement
#' the procedure produces when, by construction, no predictor has any effect.
#' Permutation is at the *cluster* level - whole country-year predictor vectors
#' are reassigned to other country-years - so the permuted data keep the
#' clustering and the collinearity of the real data and destroy only the
#' association with the response.
#'
#' @param dataset A `vs_dataset`.
#' @param base_terms Terms always retained, or a `vs_structure`.
#' @param candidates Candidate predictors. Defaults to all.
#' @param random_term Optional random-effect term.
#' @param criterion `"AIC"` or `"BIC"`.
#' @param strata Columns within which clusters are permuted. Defaults to any
#'   grouping term that is already in `base_terms` - usually `country`. Holding
#'   the base terms fixed is what makes the null match the model actually being
#'   fitted: country is a fixed effect, so the question is whether the
#'   *within-country, year-to-year* variation in a predictor tracks divergence,
#'   not whether countries differ. Permuting across countries would answer an
#'   easier question and give an anti-conservative null. Set `character()` to
#'   permute freely.
#' @param B Number of permutations.
#' @param seed Optional seed, for a reproducible null.
#' @return A `vs_calibration` object: the observed best first-round delta, the
#'   permutation null, the empirical p-value, and the delta threshold that would
#'   actually give a 5% false-positive rate on data shaped like this.
#' @examples
#' \donttest{
#' ds <- example_analysis()$dataset
#' calibrate_selection(ds, B = 99)
#' }
#' @export
calibrate_selection <- function(dataset, base_terms = NULL, candidates = NULL,
                                random_term = NULL, criterion = c("AIC", "BIC"),
                                strata = NULL, B = 500, seed = NULL) {
  stopifnot(inherits(dataset, "vs_dataset"))
  criterion <- match.arg(criterion)
  if (inherits(base_terms, "vs_structure")) {
    random_term <- random_term %||% base_terms$random_term
    base_terms <- base_terms$base_terms
  }
  base_terms <- base_terms %||% intersect(c("country", "year"), names(dataset$data))
  candidates <- setdiff(candidates %||% dataset$predictors, base_terms)
  if (!length(candidates)) vs_abort("No candidate predictors to calibrate.")
  if (!is.null(seed)) set.seed(seed)

  if (is.null(strata)) {
    strata <- Filter(function(nm) {
      nm %in% names(dataset$data) &&
        (is.factor(dataset$data[[nm]]) || is.character(dataset$data[[nm]]))
    }, intersect(base_terms, dataset$cluster))
  }
  strata <- intersect(strata, names(dataset$data))

  observed <- best_delta(dataset, base_terms, candidates, random_term, criterion)
  null <- vapply(seq_len(B), function(i) {
    best_delta(permute_clusters(dataset, candidates, strata), base_terms,
               candidates, random_term, criterion)
  }, numeric(1))
  null <- null[is.finite(null)]
  if (!length(null)) vs_abort("Every permutation failed to fit.")

  p <- (1 + sum(null <= observed)) / (1 + length(null))
  structure(
    list(observed = observed, null = null, p_value = p, B = length(null),
         criterion = criterion, candidates = candidates, strata = strata,
         base_terms = base_terms, level = dataset$level %||% "sequence",
         threshold_05 = unname(stats::quantile(null, 0.05, names = FALSE)),
         nominal_threshold = -2),
    class = "vs_calibration"
  )
}

## Best (most negative) first-round delta over all candidates.
best_delta <- function(dataset, base_terms, candidates, random_term, criterion) {
  d <- dataset$data
  resp <- dataset$response
  crit_fun <- if (criterion == "AIC") stats::AIC else stats::BIC
  base <- try(fit_model(d, resp, base_terms, random_term,
                        weights = model_weights(dataset)), silent = TRUE)
  if (inherits(base, "try-error")) return(NA_real_)
  c0 <- crit_fun(base)
  deltas <- vapply(candidates, function(p) {
    fit <- try(fit_model(d, resp, c(base_terms, p), random_term,
                         weights = model_weights(dataset)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    crit_fun(fit) - c0
  }, numeric(1))
  if (all(is.na(deltas))) return(NA_real_)
  min(deltas, na.rm = TRUE)
}

## Reassign whole cluster-level predictor vectors to other clusters, keeping the
## response and the cluster sizes exactly as they are. Donors are drawn within
## `strata`, so anything already in the model (country, typically) is held fixed.
permute_clusters <- function(dataset, candidates, strata = character()) {
  d <- dataset$data
  if (!".cluster" %in% names(d)) {
    d[candidates] <- d[sample(nrow(d)), candidates, drop = FALSE]
    dataset$data <- d
    return(dataset)
  }
  cl <- as.character(d$.cluster)
  lv <- unique(cl)
  first <- match(lv, cl)

  donor <- seq_along(lv)
  if (length(strata)) {
    block <- as.character(interaction(d[strata], drop = TRUE, sep = "\r"))[first]
    for (b in unique(block)) {
      i <- which(block == b)
      if (length(i) > 1L) donor[i] <- i[sample(length(i))]
    }
  } else {
    donor <- sample(donor)
  }
  idx <- first[donor][match(cl, lv)]
  d[candidates] <- d[idx, candidates, drop = FALSE]
  dataset$data <- d
  dataset
}

model_weights <- function(dataset) {
  if (is.null(dataset$weights)) return(NULL)
  dataset$data[[dataset$weights]]
}

#' @export
print.vs_calibration <- function(x, ...) {
  cli::cli_h3("vs_calibration")
  cli::cli_text("{x$B} cluster-level permutations at the {.strong {x$level}} level, criterion {x$criterion}.")
  if (length(x$strata)) {
    cli::cli_text("Permuted within {.field {paste(x$strata, collapse = ' x ')}}, so the base terms are held fixed.")
  }
  cli::cli_text("Observed best first-round delta: {.strong {round(x$observed, 2)}}")
  cli::cli_text("Permutation null: median {round(stats::median(x$null), 2)}, 5th percentile {round(x$threshold_05, 2)}")
  cli::cli_text("Empirical p-value: {.strong {signif(x$p_value, 3)}}")
  if (x$threshold_05 < x$nominal_threshold) {
    cli::cli_alert_warning(
      "A delta of {x$nominal_threshold} is reached by chance in more than 5% of permutations; the honest threshold here is {round(x$threshold_05, 1)}."
    )
  }
  if (x$p_value > 0.05) {
    cli::cli_alert_warning("The observed improvement is within what the search produces on permuted data. It is not evidence of an effect.")
  }
  invisible(x)
}
