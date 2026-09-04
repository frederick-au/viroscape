#' Build the tree the distance metrics imply
#'
#' Neighbour joining on the same distance matrix the rest of the package uses,
#' so the tree, the divergence values and the rate estimate are all derived from
#' one set of assumptions rather than three.
#'
#' @inheritParams divergence_rate
#' @return A `phylo` object whose tip labels are sequence ids.
#' @export
build_tree <- function(alignment, model = NULL, metric = NULL,
                       root = c("midpoint", "outgroup", "none"),
                       outgroup = NULL) {
  stopifnot(inherits(alignment, "vs_alignment"))
  root <- match.arg(root)
  if (inherits(model, "vs_model_selection")) model <- model$best
  model <- model %||% if (alignment$molecule == "protein") "JTT" else "GTR"
  metric <- metric %||% if (alignment$molecule == "protein") "ml_aa" else "ml_nt"

  d <- compute_distance(alignment$matrix, metric, model, alignment$molecule)
  tree <- ape::nj(d)
  if (root == "midpoint") {
    need_pkg("phangorn", "midpoint-root the tree")
    tree <- phangorn::midpoint(tree)
  } else if (root == "outgroup") {
    if (!length(outgroup)) vs_abort('{.arg outgroup} is required when {.code root = "outgroup"}.')
    missing_og <- setdiff(outgroup, tree$tip.label)
    if (length(missing_og)) {
      vs_abort("Outgroup tip{?s} {.val {missing_og}} {?is/are} not in the alignment.")
    }
    tree <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)
    tree <- ape::drop.tip(tree, outgroup)
  }
  attr(tree, "vs_model") <- model
  attr(tree, "vs_metric") <- metric
  tree
}

#' Assign sequences to clades
#'
#' Groups sequences by patristic distance on the tree. This exists because
#' lineage is usually the largest single source of structure in divergence data:
#' sequences in the same clade share their history, so they are not independent
#' observations of anything, and an environmental coefficient that survives
#' country and year may still be a lineage effect wearing a country's clothes.
#' Having a clade factor lets you find out - put it in `base_terms`, use it as a
#' `random_term`, or make it the cluster.
#'
#' Groups are cut from the tree itself, at a depth from the root, so every group
#' is monophyletic by construction. `k` is reached by searching for the depth
#' that produces that many groups; the number of groups only ever increases with
#' depth, so the search is well behaved, but an exact `k` is not always
#' attainable and the nearest is returned. If you have a curated clade
#' assignment - H5 clade designations, for instance - use that instead: pass it
#' straight to [add_clades()].
#'
#' A tree with no internal structure - as the bundled example has, since its
#' sequences are each drawn independently from one ancestor - will not yield
#' useful clades, and it should not: there are none to find.
#'
#' @param alignment A `vs_alignment`. Optional if `tree` is supplied.
#' @param k Target number of clades. Ignored if `h` is given.
#' @param h Depth from the root at which to cut.
#' @param tree Optional pre-built tree from [build_tree()].
#' @param ... Passed to [build_tree()].
#' @return A named factor of clade labels, one per sequence id, with the cut
#'   depth in a `height` attribute.
#' @export
assign_clades <- function(alignment = NULL, k = 6, h = NULL, tree = NULL, ...) {
  if (is.null(tree) && is.null(alignment)) {
    vs_abort("Supply either {.arg alignment} or {.arg tree}.")
  }
  tree <- tree %||% build_tree(alignment, ...)
  ntip <- ape::Ntip(tree)
  if (!is.null(k) && k >= ntip) vs_abort("{.arg k} must be smaller than the number of sequences.")
  depth <- ape::node.depth.edgelength(tree)

  groups_at <- function(hh) {
    pd <- depth[tree$edge[, 1]]
    cd <- depth[tree$edge[, 2]]
    starts <- tree$edge[, 2][pd < hh & cd >= hh]
    tips <- lapply(starts, function(nd) {
      if (nd <= ntip) nd else phangorn::Descendants(tree, nd, type = "tips")[[1]]
    })
    covered <- unlist(tips)
    ## a tip whose whole root path is shallower than the cut is its own group
    loners <- setdiff(seq_len(ntip), covered)
    c(tips, as.list(loners))
  }

  if (is.null(h)) {
    need_pkg("phangorn", "cut a tree into clades")
    lo <- 0
    hi <- max(depth)
    for (i in seq_len(60)) {
      mid <- (lo + hi) / 2
      n <- length(groups_at(mid))
      if (n < k) lo <- mid else hi <- mid
      if (n == k) break
    }
    h <- (lo + hi) / 2
  }
  need_pkg("phangorn", "cut a tree into clades")
  gl <- groups_at(h)

  idx <- integer(ntip)
  for (i in seq_along(gl)) idx[gl[[i]]] <- i
  ## label clades by size, largest first, so clade01 is the dominant lineage
  ord <- order(-tabulate(idx, nbins = length(gl)))
  relabel <- match(idx, ord)
  out <- factor(sprintf("clade%02d", relabel),
                levels = sprintf("clade%02d", sort(unique(relabel))))
  names(out) <- tree$tip.label
  attr(out, "height") <- h
  out
}

#' Add a clade factor to a dataset
#'
#' @param dataset A sequence-level `vs_dataset`.
#' @param clades A named factor or character vector of clade labels, named by
#'   sequence id - from [assign_clades()] or from your own curated assignment.
#' @param cluster Treat clade as the clustering variable as well, so
#'   [cluster_summary()] reports the lineage ICC and [aggregate_to_cluster()]
#'   collapses by clade.
#' @return The dataset with a `clade` column.
#' @export
add_clades <- function(dataset, clades, cluster = FALSE) {
  stopifnot(inherits(dataset, "vs_dataset"))
  d <- dataset$data
  if (!"id" %in% names(d)) vs_abort("The dataset has no {.field id} column to match clades on.")
  if (is.null(names(clades))) vs_abort("{.arg clades} must be named by sequence id.")
  hit <- match(d$id, names(clades))
  if (anyNA(hit)) {
    vs_warn("{sum(is.na(hit))} row{?s} had no clade assignment and will be {.val NA}.")
  }
  d$clade <- factor(as.character(clades)[hit])
  dataset$data <- d
  if (cluster) {
    dataset$cluster <- "clade"
    dataset$data$.cluster <- factor(as.character(d$clade))
    cs <- cluster_summary(dataset)
    dataset$icc <- cs$icc
    dataset$design_effect <- cs$design_effect
  }
  dataset
}

#' Generalised least squares with a phylogenetic correlation structure
#'
#' Clustering by country-year handles replication. It does not handle the
#' tree. Sequences in the same lineage are correlated by descent whatever
#' country-year they were sampled in, and divergence from a consensus is
#' precisely the kind of response that inherits that correlation. This fits the
#' same regression under a Brownian (or Pagel) covariance derived from the tree,
#' and reports the coefficients beside their ordinary least squares
#' counterparts, because what matters is not the phylogenetic fit on its own but
#' how much of an environmental effect survives it.
#'
#' The comparison is the output. If a coefficient goes to zero once the tree is
#' in the model, the association was lineage structure, and no amount of
#' cluster-robust standard error would have told you.
#'
#' @section Read the clade column first:
#' A third fit is reported alongside: the same regression with an explicit clade
#' factor. On simulated data whose divergence is driven entirely by six discrete
#' clades, with no environmental effect at all, the clade factor returns the
#' correct answer (p = 0.19) while PGLS removes 97.6% of the spurious
#' coefficient but leaves a remainder that still clears 0.05. That is not a bug
#' in the fit: Brownian and Pagel structures spread covariance smoothly along
#' branch lengths, so they approximate block structure rather than absorbing it,
#' and discrete clades are block structure. Where the two disagree, believe the
#' clade factor and treat PGLS as the secondary check. `shrinkage_pgls` is the
#' number worth reporting from it - a coefficient that loses 98% of its
#' magnitude has been explained by descent whatever its p-value says.
#'
#' @param dataset A **sequence-level** `vs_dataset` - one row per tip. Aggregated
#'   datasets have no tips to map onto and are rejected.
#' @param tree A `phylo` from [build_tree()], or a `vs_alignment` to build one
#'   from.
#' @param terms Right-hand-side terms. Defaults to the dataset's predictors plus
#'   `year`.
#' @param correlation `"brownian"` fixes the phylogenetic signal; `"pagel"`
#'   estimates Pagel's lambda, which is the more honest default when you do not
#'   know how much signal there is, at the cost of a slower and less stable fit.
#' @param clades Clade assignment for the third, categorical fit - a named
#'   factor as returned by [assign_clades()]. `NULL` derives one from the tree
#'   with `clade_k` groups; `NA` skips the clade fit.
#' @param clade_k Number of clades to derive when `clades` is `NULL`.
#' @param ... Passed to [nlme::gls()].
#' @return A `vs_phylo_fit`: the fitted `gls`, the equivalent OLS fit, the clade
#'   fit, and a `comparison` tibble of coefficients under all three.
#' @examples
#' \donttest{
#' aln <- align_sequences(example_sequences(), method = "none")
#' ds <- join_environment(distance_from_consensus(aln), example_environment())
#' phylo_gls(ds, aln, terms = c("year", "agricultural_land_pct"))
#' }
#' @export
phylo_gls <- function(dataset, tree, terms = NULL,
                      correlation = c("brownian", "pagel"),
                      clades = NULL, clade_k = 6, ...) {
  stopifnot(inherits(dataset, "vs_dataset"))
  correlation <- match.arg(correlation)
  need_pkg("nlme", "fit a phylogenetic GLS")
  if (identical(dataset$level, "cluster")) {
    vs_abort(c(
      "{.fn phylo_gls} needs one row per sequence, and this dataset is aggregated to clusters.",
      "i" = "Pass the sequence-level dataset; the tree has no tips to map onto cluster means."
    ))
  }
  if (inherits(tree, "vs_alignment")) tree <- build_tree(tree)
  if (!inherits(tree, "phylo")) vs_abort("{.arg tree} must be a {.cls phylo} or a {.cls vs_alignment}.")

  d <- as.data.frame(dataset$data)
  if (!"id" %in% names(d)) vs_abort("The dataset has no {.field id} column to match tips on.")
  shared <- intersect(d$id, tree$tip.label)
  if (length(shared) < 3L) vs_abort("Fewer than three sequences are shared between the dataset and the tree.")
  if (anyDuplicated(d$id)) vs_abort("Sequence ids must be unique to map onto tips.")
  dropped <- setdiff(tree$tip.label, shared)
  if (length(dropped)) tree <- ape::drop.tip(tree, dropped)
  d <- d[match(tree$tip.label, d$id), , drop = FALSE]
  rownames(d) <- d$id

  terms <- terms %||% unique(c(intersect("year", names(d)), dataset$predictors))
  terms <- intersect(terms, names(d))
  if (!length(terms)) vs_abort("No usable terms.")
  f <- stats::reformulate(terms, dataset$response)

  cs <- if (correlation == "pagel") {
    ape::corPagel(1, phy = tree, form = ~id)
  } else {
    ape::corBrownian(1, phy = tree, form = ~id)
  }
  fit <- nlme::gls(f, data = d, correlation = cs, method = "ML", ...)
  ols <- stats::lm(f, data = d)

  cf <- summary(fit)$tTable
  co <- stats::coef(summary(ols))
  keep <- intersect(rownames(cf), rownames(co))
  comparison <- tibble::tibble(
    term = keep,
    estimate_ols = unname(co[keep, 1]),
    p_ols = unname(co[keep, 4]),
    estimate_pgls = unname(cf[keep, 1]),
    std_error_pgls = unname(cf[keep, 2]),
    p_pgls = unname(cf[keep, 4])
  )
  comparison$survives_pgls <- comparison$p_pgls < 0.05
  comparison$shrinkage_pgls <- 1 - abs(comparison$estimate_pgls) / abs(comparison$estimate_ols)

  ## the categorical counterpart: discrete clades absorbed exactly rather than
  ## approximated by a smooth covariance
  clade_fit <- NULL
  if (!identical(clades, NA)) {
    cl <- clades %||% try(assign_clades(NULL, k = clade_k, tree = tree), silent = TRUE)
    if (!inherits(cl, "try-error") && !is.null(cl)) {
      d$clade <- factor(as.character(cl)[match(d$id, names(cl))])
      if (nlevels(d$clade) > 1L) {
        clade_fit <- stats::lm(stats::reformulate(c(terms, "clade"), dataset$response),
                               data = d)
        cc <- stats::coef(summary(clade_fit))
        k2 <- match(comparison$term, rownames(cc))
        comparison$estimate_clade <- unname(cc[k2, 1])
        comparison$p_clade <- unname(cc[k2, 4])
        comparison$shrinkage_clade <-
          1 - abs(comparison$estimate_clade) / abs(comparison$estimate_ols)
      }
    }
  }

  structure(
    list(model = fit, ols = ols, clade_model = clade_fit,
         comparison = comparison, tree = tree,
         correlation = correlation, terms = terms,
         lambda = if (correlation == "pagel") {
           tryCatch(unname(stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE)),
                    error = function(e) NA_real_)
         } else NA_real_),
    class = "vs_phylo_fit"
  )
}

#' @export
print.vs_phylo_fit <- function(x, ...) {
  cli::cli_h3("vs_phylo_fit")
  cli::cli_text("{ape::Ntip(x$tree)} tips, {x$correlation} correlation{if (!is.na(x$lambda)) paste0(' (lambda = ', signif(x$lambda, 3), ')') else ''}")
  cmp <- x$comparison
  real <- cmp$term != "(Intercept)"
  worst <- suppressWarnings(max(cmp$shrinkage_pgls[real], na.rm = TRUE))
  if (is.finite(worst)) {
    ## one decimal: 99.6% and 100% are different claims, and rounding to the
    ## latter reads as complete elimination
    cli::cli_text("Largest coefficient shrinkage once the tree is in the model: {.strong {sprintf('%.1f%%', 100 * worst)}}")
  }
  cli::cli_text("Coefficients before and after:")
  print(as.data.frame(cmp), row.names = FALSE, digits = 3)

  lost <- cmp$term[!cmp$survives_pgls & cmp$p_ols < 0.05 & real]
  if (length(lost)) {
    cli::cli_alert_warning("{.field {lost}} {?is/are} significant under OLS and not under the phylogenetic model: {?that association is/those associations are} lineage structure.")
  }
  if ("p_clade" %in% names(cmp)) {
    disagree <- cmp$term[real & !is.na(cmp$p_clade) &
                           cmp$survives_pgls & cmp$p_clade >= 0.05]
    if (length(disagree)) {
      cli::cli_alert_warning("{.field {disagree}} survive{?s/} the phylogenetic model but not an explicit clade factor.")
      cli::cli_text("{.emph A smooth covariance approximates block structure; a clade factor absorbs it. Believe the clade column.}")
    }
  }
  invisible(x)
}
