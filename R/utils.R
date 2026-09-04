#' @keywords internal
#' @importFrom rlang .data
#' @importFrom stats formula setNames
#' @importFrom utils head tail
"_PACKAGE"

## quiet R CMD check for tidy-eval column names
utils::globalVariables(c(
  "accession", "country", "year", "host", "strain", "subtype", "segment",
  "distance", "predictor", "aic", "delta_aic", "gvif", "accepted", "reason",
  "term", "estimate", "conf.low", "conf.high", "n", "model", "logLik", "df",
  "bic", "delta_bic", "round", "value", "indicator", "statistic", "p_value",
  "id", "step", "criterion_value", "delta", "sign_flip", "max_prop_shift", ".",
  "sig", "flag", "best", "label", "gvif_scaled", "iso3", "prop_shift",
  "significant_before", "is_group", "se_shift", "substantial", "destabilising",
  "lost_significance", "estimate_before", "estimate_after"
))

`%||%` <- function(x, y) if (is.null(x)) y else x

## cli captures the calling environment for glue interpolation, so the wrappers
## must forward it explicitly or messages lose access to local variables.
vs_abort <- function(msg, ..., .envir = parent.frame()) {
  cli::cli_abort(msg, ..., .envir = .envir, call = NULL)
}
vs_warn <- function(msg, ..., .envir = parent.frame()) {
  cli::cli_warn(msg, ..., .envir = .envir)
}
vs_inform <- function(msg, ..., .envir = parent.frame()) {
  cli::cli_inform(msg, ..., .envir = .envir)
}

is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

need_pkg <- function(pkg, why) {
  if (!is_installed(pkg)) {
    vs_abort(c(
      "The {.pkg {pkg}} package is required to {why}.",
      "i" = "Install it with {.run install.packages(\"{pkg}\")}."
    ))
  }
  invisible(TRUE)
}

#' Standard ISO3 codes for the countries used in the source study
#'
#' @return A named character vector mapping country name to ISO3 code.
#' @export
sea_countries <- function() {
  c(Cambodia = "KHM", Indonesia = "IDN", Malaysia = "MYS",
    Philippines = "PHL", Thailand = "THA")
}

#' Country name to ISO3 lookup used across the package
#'
#' Covers the countries most represented in avian influenza surveillance data.
#' Unmatched names are returned as `NA` with a warning.
#'
#' @param x Character vector of country names.
#' @return Character vector of ISO3 codes.
#' @export
country_to_iso3 <- function(x) {
  map <- c(
    Cambodia = "KHM", Indonesia = "IDN", Malaysia = "MYS", Philippines = "PHL",
    Thailand = "THA", Vietnam = "VNM", "Viet Nam" = "VNM", Laos = "LAO",
    "Lao PDR" = "LAO", Myanmar = "MMR", Burma = "MMR", Singapore = "SGP",
    Brunei = "BRN", "Timor-Leste" = "TLS", China = "CHN", India = "IND",
    Bangladesh = "BGD", Nepal = "NPL", Pakistan = "PAK", "Sri Lanka" = "LKA",
    Japan = "JPN", "South Korea" = "KOR", "Korea" = "KOR", Taiwan = "TWN",
    "Hong Kong" = "HKG", Mongolia = "MNG", Egypt = "EGY", Nigeria = "NGA",
    "South Africa" = "ZAF", Turkey = "TUR", Russia = "RUS", Germany = "DEU",
    France = "FRA", Netherlands = "NLD", "United Kingdom" = "GBR",
    "United States" = "USA", USA = "USA", Canada = "CAN", Mexico = "MEX",
    Chile = "CHL", Peru = "PER", Brazil = "BRA", Australia = "AUS"
  )
  out <- unname(map[as.character(x)])
  if (anyNA(out) && any(!is.na(x))) {
    bad <- unique(x[is.na(out) & !is.na(x)])
    if (length(bad)) {
      vs_warn("No ISO3 code for {.val {bad}}; supply {.arg iso3} manually.")
    }
  }
  out
}

#' Z-standardise numeric columns
#'
#' Standardisation makes regression coefficients comparable across predictors
#' measured on different scales, which is what makes the effect sizes reported
#' by [forward_select()] interpretable side by side.
#'
#' @param data A data frame.
#' @param cols Character vector of columns to scale. Defaults to all numeric
#'   columns except those in `exclude`.
#' @param exclude Columns never to scale.
#' @return `data` with the selected columns standardised. The centre and scale
#'   used are kept in the `"vs_scaling"` attribute so predictions can be mapped
#'   back to the original units.
#' @export
standardise <- function(data, cols = NULL, exclude = c("year", "distance")) {
  if (is.null(cols)) {
    cols <- names(data)[vapply(data, is.numeric, logical(1))]
    cols <- setdiff(cols, exclude)
  }
  cols <- intersect(cols, names(data))
  pars <- list()
  for (nm in cols) {
    mu <- mean(data[[nm]], na.rm = TRUE)
    sdv <- stats::sd(data[[nm]], na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) {
      vs_warn("Column {.field {nm}} has zero variance and was left unscaled.")
      next
    }
    data[[nm]] <- (data[[nm]] - mu) / sdv
    pars[[nm]] <- c(center = mu, scale = sdv)
  }
  attr(data, "vs_scaling") <- pars
  data
}
