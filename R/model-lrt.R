#' Sequential likelihood ratio tests across the nested model chain
#'
#' The models along a forward-selection path are nested, so the likelihood ratio
#' statistic can be computed. Nesting is necessary for a valid test but it is
#' not sufficient, and this is the trap: the terms were *chosen on the same
#' data*, by searching over candidates and keeping the winner. A statistic
#' computed on the selected path does not have the chi-squared null distribution
#' the p-value assumes, and the p-values it produces are anti-conservative -
#' often severely so when the data are clustered.
#'
#' They are reported as `nominal_p` for exactly that reason. Treat them as a
#' description of the fitted path, not as evidence against a null hypothesis,
#' unless the model was pre-specified or the search was calibrated with
#' [calibrate_selection()].
#'
#' @param selection A `vs_selection` object from [forward_select()].
#' @param alpha Threshold used for the `below_alpha` column. It is deliberately
#'   not called `significant`.
#' @return A tibble with one row per addition: `added`, `statistic`, `df`,
#'   `nominal_p` and `below_alpha`.
#' @export
lrt_chain <- function(selection, alpha = 0.05) {
  stopifnot(inherits(selection, "vs_selection"))
  chain <- selection$chain
  if (length(chain) < 2) {
    return(tibble::tibble(added = character(), statistic = numeric(),
                          df = numeric(), nominal_p = numeric(),
                          below_alpha = logical()))
  }
  added <- setdiff(selection$terms, selection$base_terms)
  rows <- lapply(seq_len(length(chain) - 1), function(i) {
    a <- suppressWarnings(try(stats::anova(chain[[i]], chain[[i + 1]]), silent = TRUE))
    if (inherits(a, "try-error")) {
      return(tibble::tibble(added = added[i], statistic = NA_real_,
                            df = NA_real_, nominal_p = NA_real_))
    }
    extract_anova_row(a, added[i])
  })
  out <- dplyr::bind_rows(rows)
  out$below_alpha <- !is.na(out$nominal_p) & out$nominal_p < alpha
  attr(out, "post_selection") <- TRUE
  out
}

extract_anova_row <- function(a, label) {
  a <- as.data.frame(a)
  last <- nrow(a)
  pcol <- grep("^Pr\\(", names(a), value = TRUE)
  scol <- intersect(c("F", "Chisq", "F value", "LRT"), names(a))
  dcol <- intersect(c("Df", "Chi Df", "Res.Df"), names(a))
  tibble::tibble(
    added = label,
    statistic = if (length(scol)) as.numeric(a[[scol[1]]][last]) else NA_real_,
    df = if (length(dcol)) as.numeric(a[[dcol[1]]][last]) else NA_real_,
    nominal_p = if (length(pcol)) as.numeric(a[[pcol[1]]][last]) else NA_real_
  )
}

#' Tidy coefficient table with confidence intervals
#'
#' Standard errors are cluster-robust by default. The environmental predictors
#' vary only at the country-year level while the response is per sequence, so
#' the model-based standard errors are too small by roughly the square root of
#' the design effect (see [cluster_summary()]). `sandwich::vcovCL` fixes the
#' variance; it does not fix the model selection that chose the terms.
#'
#' With few clusters - and five countries times twenty years is few - even the
#' cluster-robust estimator is optimistic. `clubSandwich::vcovCR` with `type =
#' "CR2"` is the better choice there, and can be supplied through `vcov`.
#'
#' @param model A fitted model, or a `vs_selection` (which carries its dataset,
#'   and therefore the cluster definition, with it).
#' @param level Confidence level.
#' @param vcov `"cluster"` for `sandwich::vcovCL`, `"classical"` for the
#'   model-based matrix, or a variance-covariance matrix / function supplied
#'   directly.
#' @param dataset Optional `vs_dataset` supplying the cluster definition when
#'   `model` is a bare fitted model.
#' @return A tibble with `term`, `estimate`, `std_error`, `statistic`,
#'   `p_value`, `conf.low`, `conf.high`, and `unit`/`sd` describing what one
#'   unit of each predictor means.
#' @export
tidy_coefficients <- function(model, level = 0.95,
                              vcov = c("cluster", "classical"),
                              dataset = NULL) {
  if (inherits(model, "vs_selection")) {
    dataset <- dataset %||% model$dataset
    model <- model$model
  }
  if (is.character(vcov)) vcov <- match.arg(vcov)

  est <- safe_coef(model)
  s <- stats::coef(summary(model))
  se <- s[, 2]
  vcov_type <- "classical"

  V <- NULL
  if (is.matrix(vcov)) {
    V <- vcov; vcov_type <- "supplied"
  } else if (is.function(vcov)) {
    V <- vcov(model); vcov_type <- "supplied"
  } else if (identical(vcov, "cluster")) {
    cl <- cluster_vector(model, dataset)
    if (is.null(cl)) {
      vs_warn("No cluster definition available; falling back to model-based standard errors.")
    } else if (!is_installed("sandwich")) {
      vs_warn(c("{.pkg sandwich} is not installed; falling back to model-based standard errors.",
                "i" = "These are anti-conservative on clustered data. {.code install.packages(\"sandwich\")}"))
    } else if (inherits(model, "merMod")) {
      vs_warn("Cluster-robust variances are not applied to mixed models here; using the model-based matrix.")
    } else {
      V <- try(sandwich::vcovCL(model, cluster = cl, type = "HC1"), silent = TRUE)
      if (inherits(V, "try-error")) {
        vs_warn("{.fn sandwich::vcovCL} failed; falling back to model-based standard errors.")
        V <- NULL
      } else {
        vcov_type <- sprintf("cluster-robust (%d clusters)", nlevels(droplevels(as.factor(cl))))
      }
    }
  }

  if (!is.null(V)) {
    keep <- intersect(names(est), rownames(V))
    se <- sqrt(diag(V))[match(names(est), rownames(V))]
    names(se) <- names(est)
  }
  stat <- est / se
  dfres <- if (inherits(model, "merMod")) Inf else stats::df.residual(model)
  if (!is.null(V) && !is.null(dataset) && length(dataset$cluster)) {
    cl <- cluster_vector(model, dataset)
    if (!is.null(cl)) dfres <- max(nlevels(droplevels(as.factor(cl))) - 1L, 1L)
  }
  pv <- 2 * stats::pt(-abs(stat), df = dfres)
  crit <- stats::qt(1 - (1 - level) / 2, df = dfres)

  sds <- scaling_sd(dataset, names(est))
  out <- tibble::tibble(
    term = names(est),
    estimate = unname(est),
    std_error = unname(se),
    statistic = unname(stat),
    p_value = unname(pv),
    conf.low = unname(est - crit * se),
    conf.high = unname(est + crit * se),
    unit = unname(ifelse(is.na(sds), "natural unit", "standard deviation")),
    sd = unname(sds)
  )
  attr(out, "vcov_type") <- vcov_type
  attr(out, "df") <- dfres
  out
}

cluster_vector <- function(model, dataset) {
  if (is.null(dataset) || !length(dataset$cluster)) return(NULL)
  d <- dataset$data
  if (!".cluster" %in% names(d)) return(NULL)
  mf <- try(stats::model.frame(model), silent = TRUE)
  cl <- d$.cluster
  if (!inherits(mf, "try-error") && nrow(mf) != length(cl)) {
    idx <- as.integer(rownames(mf))
    if (!anyNA(idx) && max(idx) <= length(cl)) cl <- cl[idx] else return(NULL)
  }
  droplevels(as.factor(cl))
}

scaling_sd <- function(dataset, terms) {
  out <- stats::setNames(rep(NA_real_, length(terms)), terms)
  if (is.null(dataset) || !isTRUE(dataset$scaled)) return(out)
  pars <- attr(dataset$data, "vs_scaling")
  if (!length(pars)) return(out)
  for (nm in names(pars)) {
    if (nm %in% terms) out[[nm]] <- unname(pars[[nm]]["scale"])
  }
  out
}

#' Convert standardised coefficients back to natural units
#'
#' [join_environment()] standardises predictors by default, so a coefficient is
#' the change in divergence per standard deviation of the predictor. That is the
#' right default for comparing predictors with each other and the wrong one for
#' reporting an effect size, and confusing the two is a documented hazard of
#' this analysis. This converts.
#'
#' @param coefficients A tibble from [tidy_coefficients()].
#' @param dataset The `vs_dataset` the model was fitted to. Optional if
#'   `coefficients` already carries an `sd` column.
#' @return The same tibble with `estimate`, `std_error` and the interval bounds
#'   expressed per natural unit of each predictor.
#' @export
unscale <- function(coefficients, dataset = NULL) {
  sds <- coefficients$sd %||% scaling_sd(dataset, coefficients$term)
  if (all(is.na(sds))) {
    vs_warn("Nothing to unscale: no standardisation parameters were recorded.")
    return(coefficients)
  }
  f <- ifelse(is.na(sds), 1, 1 / sds)
  for (col in c("estimate", "std_error", "conf.low", "conf.high")) {
    if (col %in% names(coefficients)) coefficients[[col]] <- coefficients[[col]] * f
  }
  coefficients$unit <- rep("natural unit", length(sds))
  coefficients$sd <- sds
  attr(coefficients, "unscaled") <- TRUE
  coefficients
}

#' Analysis of variance table for the final model
#'
#' @param selection A `vs_selection` object.
#' @return A tibble of sequential (Type I) tests.
#' @export
anova_table <- function(selection) {
  m <- if (inherits(selection, "vs_selection")) selection$model else selection
  a <- suppressWarnings(try(stats::anova(m), silent = TRUE))
  if (inherits(a, "try-error")) return(tibble::tibble())
  df <- as.data.frame(a)
  df$term <- rownames(df)
  tibble::as_tibble(df[c("term", setdiff(names(df), "term"))])
}
