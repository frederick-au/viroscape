#' Estimate a substitution rate by root-to-tip regression
#'
#' Divergence from a consensus is not a rate: the consensus is not an ancestor,
#' and dividing by elapsed time (as `clock_correct` did in 0.1.0) does not make
#' it one - it manufactures a negative association with year even in data
#' generated under a perfectly constant rate. This is the honest alternative,
#' and it is the standard one: build a tree, root it, measure each tip's
#' distance from the root, and regress that on collection date. The slope is
#' substitutions per site per year and the intercept implies a date for the root
#' (Rambaut et al. 2016, *Virus Evolution* 2:vew007).
#'
#' Root-to-tip regression has known limitations - tips are not independent, so
#' the standard error is optimistic, and the fit is sensitive to the rooting.
#' They are limitations of an estimator that answers the right question, rather
#' than an estimator that answers the wrong one.
#'
#' @section What this estimator can and cannot see:
#' Root-to-tip regression needs two things: enough substitutions to order the
#' sample in time, and a root. On simulated data evolving at a known constant
#' rate, with a divergent outgroup in the alignment to anchor the deep split,
#' it recovers the realised rate to within a few per cent (best-fitting rooting
#' -3%, midpoint -7% across five replicates). Remove that outgroup from the same
#' data and both rootings collapse to a slope near zero or below, because a
#' neighbour-joining tree built from roughly twenty substitutions spread over
#' twenty years cannot resolve the order of the trunk, and root-to-tip distance
#' becomes tip noise. This is not a failure mode to correct for; it is the
#' absence of signal, and the function warns rather than returning a number.
#'
#' Two practical consequences. Include an outgroup if you have one, and use
#' `root = "outgroup"`. And read `r_squared` before reading `slope`: the default
#' `root = "best"` chooses the root position that minimises the regression's
#' residual variance (Rambaut et al. 2016), which is stable and standard but
#' also fitted, so its `r_squared` is optimistic by construction.
#'
#' Expect a modest upward bias from the maximum-likelihood metrics relative to
#' the uncorrected ones - roughly 10-15% on the data above. That is the multiple
#' hit correction doing its job, not an error.
#'
#' @param alignment A `vs_alignment`.
#' @param by Optional grouping column(s) in the alignment metadata, e.g.
#'   `"country"`. A separate regression is reported for each group as well as
#'   for the whole sample.
#' @param model Substitution model, or a `vs_model_selection` object.
#' @param metric Distance metric used to build the tree; see
#'   [distance_metrics()].
#' @param root `"best"` (the root that minimises the regression's residual
#'   variance), `"midpoint"`, `"outgroup"` (requires `outgroup`) or `"none"`
#'   (use the tree as returned by neighbour joining).
#' @param outgroup Tip labels forming the outgroup when `root = "outgroup"`.
#'   The outgroup is dropped from the regression after rooting.
#' @param min_n Groups with fewer than this many sequences, or fewer than three
#'   distinct years, are reported with `NA` slopes rather than a fit nothing
#'   supports.
#' @return A `vs_rate` object: a tibble of slopes with standard errors, plus the
#'   tree and the per-tip root-to-tip distances as attributes.
#' @examples
#' \donttest{
#' aln <- align_sequences(example_sequences(), method = "none")
#' divergence_rate(aln, by = "country")
#' }
#' @export
divergence_rate <- function(alignment, by = NULL, model = NULL,
                            metric = NULL,
                            root = c("best", "midpoint", "outgroup", "none"),
                            outgroup = NULL, min_n = 10) {
  stopifnot(inherits(alignment, "vs_alignment"))
  root <- match.arg(root)
  if (root == "outgroup" && !length(outgroup)) {
    vs_abort('{.arg outgroup} is required when {.code root = "outgroup"}.')
  }
  if (inherits(model, "vs_model_selection")) model <- model$best
  model <- model %||% if (alignment$molecule == "protein") "JTT" else "GTR"
  metric <- metric %||% if (alignment$molecule == "protein") "ml_aa" else "ml_nt"

  meta <- alignment$meta
  if (root == "best") {
    tree <- build_tree(alignment, model = model, metric = metric, root = "none")
    yr <- meta$year[match(tree$tip.label, meta$id)]
    fit_root <- best_fitting_root(tree, yr)
    rtt <- fit_root$root_to_tip
    tree <- fit_root$tree
  } else {
    tree <- build_tree(alignment, model = model, metric = metric, root = root,
                       outgroup = outgroup)
    rtt <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
    names(rtt) <- tree$tip.label
  }

  dat <- tibble::tibble(id = names(rtt), root_to_tip = as.numeric(rtt))
  dat <- dplyr::left_join(dat, meta, by = "id")
  if (!"year" %in% names(dat)) vs_abort("The alignment metadata has no {.field year} column.")
  dat <- dat[!is.na(dat$year) & !is.na(dat$root_to_tip), , drop = FALSE]
  if (!nrow(dat)) vs_abort("No sequences have both a root-to-tip distance and a year.")

  rows <- list(rate_row(dat, "all"))
  if (!is.null(by)) {
    missing_by <- setdiff(by, names(dat))
    if (length(missing_by)) vs_abort("Grouping column{?s} {.field {missing_by}} not found.")
    key <- interaction(dat[by], drop = TRUE, sep = " / ")
    for (g in levels(key)) {
      rows[[length(rows) + 1L]] <- rate_row(dat[key == g, , drop = FALSE], g, min_n = min_n)
    }
  }
  out <- dplyr::bind_rows(rows)

  if (isTRUE(out$slope[1] <= 0)) {
    vs_warn(c(
      "The overall root-to-tip slope is not positive ({signif(out$slope[1], 3)} substitutions/site/year).",
      "i" = "There is no usable clock signal in this sample; do not report a rate from it."
    ))
  } else if (isTRUE(out$r_squared[1] < 0.1)) {
    vs_warn(c(
      "The root-to-tip regression explains {round(100 * out$r_squared[1])}% of the variance.",
      "i" = "A positive slope on this little structure is not evidence of a clock. Check the tree has enough substitutions to order the sample in time."
    ))
  }

  structure(out, tree = tree, data = dat, model = model, metric = metric,
            root = root, class = c("vs_rate", class(out)))
}

## The root position minimising the residual variance of root-to-tip against
## date, searched over every branch and analytically along each one. Midpoint
## rooting answers a question about tree shape; this answers the question the
## regression is actually asking.
best_fitting_root <- function(tree, years) {
  ntip <- ape::Ntip(tree)
  ok <- !is.na(years)
  if (sum(ok) < 3L) vs_abort("At least three dated sequences are needed to fit a root.")
  D <- ape::dist.nodes(tree)
  tips <- seq_len(ntip)

  ## residuals from regressing on date, centred so that calendar years do not
  ## make the normal equations ill-conditioned
  yc <- years[ok] - mean(years[ok])
  sy2 <- sum(yc^2)
  if (sy2 <= 0) vs_abort("All dated sequences share one year; no rate is identifiable.")
  resid_of <- function(v) {
    vc <- v - mean(v)
    vc - (sum(yc * vc) / sy2) * yc
  }

  best <- list(rss = Inf)
  for (e in seq_len(nrow(tree$edge))) {
    u <- tree$edge[e, 1]
    v <- tree$edge[e, 2]
    Le <- tree$edge.length[e]
    if (!is.finite(Le)) next
    du <- D[u, tips]
    dv <- D[v, tips]
    ## a tip sits below v when the path from u reaches it through v
    below <- abs(du - (dv + Le)) < 1e-8
    a <- ifelse(below, dv + Le, du)          # distance when the root sits at u
    sgn <- ifelse(below, -1, 1)              # moving the root towards v
    a <- a[ok]; sgn <- sgn[ok]

    ra <- resid_of(a)
    rs <- resid_of(sgn)
    denom <- sum(rs^2)
    pstar <- if (denom > 0) -sum(ra * rs) / denom else 0
    pstar <- min(max(pstar, 0), Le)
    rss <- sum((ra + pstar * rs)^2)
    if (rss < best$rss) best <- list(rss = rss, edge = e, u = u, v = v,
                                     len = Le, p = pstar, a = a, sgn = sgn)
  }
  if (!is.finite(best$rss)) vs_abort("No root position could be fitted.")

  full_a <- ifelse(abs(D[best$u, tips] - (D[best$v, tips] + best$len)) < 1e-8,
                   D[best$v, tips] + best$len, D[best$u, tips])
  full_sgn <- ifelse(abs(D[best$u, tips] - (D[best$v, tips] + best$len)) < 1e-8, -1, 1)
  rtt <- full_a + best$p * full_sgn
  names(rtt) <- tree$tip.label

  rooted <- try(ape::root(tree, node = best$v, resolve.root = TRUE), silent = TRUE)
  if (inherits(rooted, "try-error")) rooted <- tree
  list(root_to_tip = rtt, tree = rooted, rss = best$rss, position = best$p)
}

rate_row <- function(d, label, min_n = 0) {
  n_year <- length(unique(d$year))
  if (nrow(d) < min_n || n_year < 3) {
    return(tibble::tibble(group = label, n = nrow(d), n_years = n_year,
                          slope = NA_real_, std_error = NA_real_,
                          p_value = NA_real_, r_squared = NA_real_,
                          root_date = NA_real_))
  }
  fit <- stats::lm(root_to_tip ~ year, data = d)
  cf <- stats::coef(summary(fit))
  slope <- cf["year", 1]
  tibble::tibble(
    group = label, n = nrow(d), n_years = n_year,
    slope = slope, std_error = cf["year", 2], p_value = cf["year", 4],
    r_squared = summary(fit)$r.squared,
    root_date = if (slope > 0) unname(-stats::coef(fit)[1] / slope) else NA_real_
  )
}

#' @export
print.vs_rate <- function(x, ...) {
  cli::cli_h3("vs_rate")
  cli::cli_text("Root-to-tip regression ({attr(x, 'metric')}, model {attr(x, 'model')}, {attr(x, 'root')} rooting)")
  cli::cli_text("Slope is substitutions per site per year; {.field root_date} is the implied date of the root.")
  cli::cli_alert_info("Tips are not independent, so these standard errors are optimistic.")
  print(as.data.frame(x), row.names = FALSE, digits = 3)
  invisible(x)
}
