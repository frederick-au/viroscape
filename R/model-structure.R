#' Compare structural specifications for the grouping variable
#'
#' Before any environmental predictor is considered, the country term has to be
#' settled: absent, a fixed effect, or a random intercept. A random effect is
#' the more natural specification but its variance component is only reliable
#' with roughly ten or more groups, so this function reports the group count
#' alongside the comparison and warns when the random-effect option is being
#' fitted on too few.
#'
#' All models are fitted by maximum likelihood, not REML, because REML
#' likelihoods are not comparable across models with different fixed effects.
#'
#' @param dataset A `vs_dataset` from [join_environment()].
#' @param group Grouping variable, usually `"country"`.
#' @param time Time covariate, usually `"year"`. Set `NULL` to skip.
#' @param criterion `"AIC"` or `"BIC"`.
#' @param force_time Retain the time covariate whether or not it improves fit.
#'   This is the default, and it is deliberate: divergence from a consensus
#'   increases mechanically with sampling date, so year is a confounder of every
#'   environmental predictor that trends over time, not a candidate competing
#'   for a place in the model. A confounder earns its place by being a
#'   confounder, not by improving AIC. Set `FALSE` to restore the 0.1.0
#'   behaviour, which dropped year whenever it failed to improve AIC by 2.
#' @return A `vs_structure` object with the comparison `table`, the chosen
#'   `base_terms`, and whether a random effect was selected.
#' @export
select_structure <- function(dataset, group = "country", time = "year",
                             criterion = c("AIC", "BIC"), force_time = TRUE) {
  stopifnot(inherits(dataset, "vs_dataset"))
  criterion <- match.arg(criterion)
  d <- dataset$data
  resp <- dataset$response
  n_groups <- length(unique(d[[group]]))

  wt <- model_weights(dataset)
  fits <- list(null = fit_model(d, resp, character(), weights = wt))
  fits$fixed <- fit_model(d, resp, group, weights = wt)
  if (is_installed("lme4") && n_groups >= 3) {
    fits$random <- try(fit_model(d, resp, character(), sprintf("(1|%s)", group),
                                 weights = wt), silent = TRUE)
    if (inherits(fits$random, "try-error")) fits$random <- NULL
  }

  tab <- dplyr::bind_rows(lapply(names(fits), function(nm) {
    tibble::tibble(structure = nm, logLik = as.numeric(stats::logLik(fits[[nm]])),
                   df = stats::AIC(fits[[nm]]) |> attr("df") %||% NA_real_,
                   AIC = stats::AIC(fits[[nm]]), BIC = stats::BIC(fits[[nm]]))
  }))
  tab$df <- vapply(fits[tab$structure], function(m) attr(stats::logLik(m), "df"), numeric(1))
  tab <- tab[order(tab[[criterion]]), ]
  tab$delta <- tab[[criterion]] - min(tab[[criterion]])

  chosen <- tab$structure[1]
  if (chosen == "random" && n_groups < 10) {
    vs_warn(c(
      "A random intercept for {.field {group}} was selected on {criterion} but there are only {n_groups} group{?s}.",
      "i" = "Variance components are generally unreliable below about 10 groups; consider {.val fixed}."
    ))
  }

  base_terms <- switch(chosen,
    null = character(),
    fixed = group,
    random = character()
  )
  random_term <- if (chosen == "random") sprintf("(1|%s)", group) else NULL

  time_test <- NULL
  if (!is.null(time) && time %in% names(d)) {
    m0 <- fit_model(d, resp, base_terms, random_term, weights = wt)
    m1 <- fit_model(d, resp, c(base_terms, time), random_term, weights = wt)
    time_test <- tibble::tibble(
      added = time,
      AIC_before = stats::AIC(m0), AIC_after = stats::AIC(m1),
      delta_AIC = stats::AIC(m1) - stats::AIC(m0)
    )
    lrt <- suppressWarnings(try(stats::anova(m0, m1), silent = TRUE))
    if (!inherits(lrt, "try-error")) {
      p <- utils::tail(lrt[[ncol(lrt)]], 1)
      time_test$p_value <- as.numeric(p)
    }
    if (force_time) {
      base_terms <- c(base_terms, time)
      if (!isTRUE(time_test$delta_AIC < -2)) {
        vs_inform("{.field {time}} did not improve fit (dAIC = {round(time_test$delta_AIC, 1)}) but is retained as a confounder.")
      }
    } else if (isTRUE(time_test$delta_AIC < -2)) {
      base_terms <- c(base_terms, time)
    } else {
      vs_warn(c(
        "{.field {time}} was dropped because it did not improve fit (dAIC = {round(time_test$delta_AIC, 1)}).",
        "i" = "Divergence accumulates with sampling date, so {.field {time}} confounds any predictor with a trend. Prefer {.code force_time = TRUE}."
      ))
    }
  }

  structure(
    list(table = tab, chosen = chosen, base_terms = base_terms,
         random_term = random_term, n_groups = n_groups,
         time_test = time_test, group = group, criterion = criterion),
    class = "vs_structure"
  )
}

fit_model <- function(data, response, terms, random_term = NULL, reml = FALSE,
                      weights = NULL) {
  rhs <- c(terms, random_term)
  if (!length(rhs)) rhs <- "1"
  f <- stats::as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
  if (!is.null(random_term)) {
    need_pkg("lme4", "fit random-effect models")
    return(lme4::lmer(f, data = data, REML = reml, weights = weights))
  }
  stats::lm(f, data = data, weights = weights)
}

#' @export
print.vs_structure <- function(x, ...) {
  cli::cli_h3("vs_structure")
  cli::cli_text("Grouping variable {.field {x$group}} with {x$n_groups} level{?s}")
  print(as.data.frame(x$table), row.names = FALSE, digits = 6)
  cli::cli_text("Chosen: {.strong {x$chosen}}")
  cli::cli_text("Base terms: {.val {if (length(x$base_terms)) x$base_terms else 'intercept only'}}")
  if (!is.null(x$time_test)) {
    cli::cli_text("Time covariate: dAIC = {round(x$time_test$delta_AIC, 1)}")
  }
  invisible(x)
}
