#' Generalised variance inflation factors
#'
#' Standard VIF is not applicable to a predictor occupying several degrees of
#' freedom, such as a multi-level country factor. The generalised VIF treats
#' each predictor as a block:
#' \deqn{GVIF_j = \frac{\det(R_{jj})\det(R_{(-j)(-j)})}{\det(R)}}
#' and \eqn{GVIF^{1/(2 df_j)}} rescales it so that values are comparable across
#' predictors with differing degrees of freedom (it is the analogue of
#' \eqn{\sqrt{VIF}}, so a threshold of 10 corresponds to VIF = 100).
#'
#' Uses \pkg{car} when available and falls back to an internal implementation
#' otherwise.
#'
#' @param model A fitted model.
#' @return A tibble with `term`, `gvif`, `df` and `gvif_scaled`.
#' @export
gvif_table <- function(model) {
  if (is_installed("car")) {
    v <- try(car::vif(model), silent = TRUE)
    if (!inherits(v, "try-error")) return(tidy_car_vif(v))
  }
  gvif_manual(model)
}

tidy_car_vif <- function(v) {
  if (is.null(dim(v))) {
    return(tibble::tibble(term = names(v), gvif = as.numeric(v),
                          df = 1, gvif_scaled = sqrt(as.numeric(v))))
  }
  tibble::tibble(
    term = rownames(v),
    gvif = as.numeric(v[, 1]),
    df = as.numeric(v[, 2]),
    gvif_scaled = as.numeric(v[, 3])
  )
}

gvif_manual <- function(model) {
  V <- try(stats::vcov(model), silent = TRUE)
  if (inherits(V, "try-error")) return(empty_gvif())
  V <- as.matrix(V)
  mm <- stats::model.matrix(model)
  asgn <- attr(mm, "assign")
  nms <- colnames(mm)
  if (length(nms) && nms[1] == "(Intercept)") {
    keep <- -1
    V <- V[keep, keep, drop = FALSE]
    asgn <- asgn[keep]
  }
  if (!length(asgn) || ncol(V) < 2) return(empty_gvif())
  R <- stats::cov2cor(V)
  detR <- det(R)
  labs <- attr(stats::terms(model), "term.labels")
  ids <- sort(unique(asgn))
  rows <- lapply(ids, function(j) {
    idx <- which(asgn == j)
    if (length(idx) == ncol(R)) return(NULL)
    g <- det(as.matrix(R[idx, idx, drop = FALSE])) *
      det(as.matrix(R[-idx, -idx, drop = FALSE])) / detR
    tibble::tibble(term = labs[j], gvif = g, df = length(idx),
                   gvif_scaled = g^(1 / (2 * length(idx))))
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) empty_gvif() else out
}

empty_gvif <- function() {
  tibble::tibble(term = character(), gvif = numeric(),
                 df = numeric(), gvif_scaled = numeric())
}

#' Compare coefficients before and after adding a predictor
#'
#' Elevated GVIF alone does not establish that a model is unusable. What
#' matters is whether adding a predictor destabilises the estimates already in
#' the model: a coefficient changing sign, or shifting by a large proportion,
#' is direct evidence that the two predictors are not separately identified.
#'
#' Two refinements keep the rule from firing on changes that carry no
#' information.
#'
#' First, a sign flip only counts when the coefficient also loses significance.
#' A term estimated at effectively zero will change sign on any perturbation,
#' and rejecting a predictor for that is noise, not diagnosis.
#'
#' Second, contrasts belonging to a grouping factor named in `group_terms` are
#' exempt. Adding a predictor that varies at the group level necessarily
#' re-partitions variance between that predictor and the group contrasts, and
#' when the predictor genuinely explains the group differences those contrasts
#' *should* collapse towards zero. Treating that as instability penalises the
#' predictor precisely when it works. What matters is whether the terms you
#' intend to interpret stay put.
#'
#' @param before,after Fitted models, `after` nesting `before`.
#' @param threshold Proportional shift treated as substantial (0.5 = 50\%).
#' @param alpha Significance level used to decide whether a coefficient was
#'   meaningfully non-zero before the addition.
#' @param group_terms Names of terms whose contrasts should be exempt from the
#'   proportional-shift rule, typically the grouping factor.
#' @param ignore_intercept Exclude the intercept from the comparison.
#' @return A tibble of shared terms with `estimate_before`, `estimate_after`,
#'   `prop_shift`, `sign_flip`, `lost_significance` and `destabilising`; plus
#'   attributes `any_destabilising`, `any_sign_flip` and `max_prop_shift`.
#' @export
coef_stability <- function(before, after, threshold = 0.5, alpha = 0.05,
                           group_terms = character(),
                           ignore_intercept = TRUE) {
  cb <- safe_coef(before)
  ca <- safe_coef(after)
  shared <- intersect(names(cb), names(ca))
  if (ignore_intercept) shared <- setdiff(shared, "(Intercept)")
  if (!length(shared)) {
    out <- tibble::tibble(term = character(), estimate_before = numeric(),
                          estimate_after = numeric(), prop_shift = numeric(),
                          sign_flip = logical(), lost_significance = logical(),
                          is_group = logical(), destabilising = logical())
    attr(out, "any_destabilising") <- FALSE
    attr(out, "any_sign_flip") <- FALSE
    attr(out, "max_prop_shift") <- 0
    return(out)
  }
  pb <- safe_pvals(before)[shared]
  pa <- safe_pvals(after)[shared]
  seb <- safe_se(before)[shared]
  out <- tibble::tibble(
    term = shared,
    estimate_before = unname(cb[shared]),
    estimate_after = unname(ca[shared]),
    prop_shift = abs(unname(ca[shared]) - unname(cb[shared])) /
      pmax(abs(unname(cb[shared])), .Machine$double.eps),
    sign_flip = sign(unname(cb[shared])) != sign(unname(ca[shared])) &
      unname(cb[shared]) != 0 & unname(ca[shared]) != 0,
    lost_significance = !is.na(pb) & !is.na(pa) & pb < 0.05 & pa >= 0.05
  )
  ## A proportional shift is unstable in the literal sense when the baseline
  ## coefficient is near zero, so it is paired with the shift measured in
  ## standard errors: a move has to be both proportionally large and larger
  ## than the uncertainty it is being judged against.
  out$se_shift <- abs(out$estimate_after - out$estimate_before) /
    pmax(unname(seb), .Machine$double.eps)
  out$significant_before <- !is.na(pb) & pb < alpha
  out$is_group <- if (length(group_terms)) {
    grepl(paste0("^(", paste(group_terms, collapse = "|"), ")"), out$term)
  } else rep(FALSE, nrow(out))
  out$substantial <- out$prop_shift > threshold & out$se_shift > 1
  ## a sign flip counts only when significance is lost with it; a large shift
  ## counts only for a term that is not a grouping contrast
  out$destabilising <- !out$is_group &
    ((out$sign_flip & out$lost_significance) |
       (out$substantial & out$significant_before))
  attr(out, "any_destabilising") <- any(out$destabilising, na.rm = TRUE)
  attr(out, "any_sign_flip") <- any(out$sign_flip & out$significant_before &
                                      !out$is_group, na.rm = TRUE)
  attr(out, "max_prop_shift") <- zero_if_infinite(
    suppressWarnings(max(out$prop_shift[out$significant_before & !out$is_group],
                         na.rm = TRUE)))
  attr(out, "any_lost_significance") <- any(out$lost_significance, na.rm = TRUE)
  out
}

safe_se <- function(m) {
  s <- try(stats::coef(summary(m)), silent = TRUE)
  if (inherits(s, "try-error")) {
    return(stats::setNames(rep(NA_real_, length(safe_coef(m))), names(safe_coef(m))))
  }
  stats::setNames(s[, 2], rownames(s))
}

zero_if_infinite <- function(x) if (!is.finite(x)) 0 else x

safe_coef <- function(m) {
  if (inherits(m, "merMod")) return(lme4::fixef(m))
  stats::coef(m)
}

safe_pvals <- function(m) {
  s <- try(stats::coef(summary(m)), silent = TRUE)
  if (inherits(s, "try-error") || !"Pr(>|t|)" %in% colnames(s)) {
    return(stats::setNames(rep(NA_real_, length(safe_coef(m))), names(safe_coef(m))))
  }
  stats::setNames(s[, "Pr(>|t|)"], rownames(s))
}
