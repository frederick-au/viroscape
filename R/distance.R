#' Available divergence metrics
#'
#' Amino acid metrics reproduce the published analysis. Nucleotide metrics are
#' the extension that recovers synonymous substitutions, which amino acid
#' distances cannot see.
#'
#' @return A tibble describing each metric.
#' @export
distance_metrics <- function() {
  tibble::tibble(
    metric = c("ml_aa", "p_aa", "ml_nt", "tn93", "k80", "raw"),
    molecule = c("protein", "protein", "nucleotide", "nucleotide",
                 "nucleotide", "both"),
    description = c(
      "Maximum likelihood distance under the selected AA substitution matrix",
      "Uncorrected proportion of differing amino acid sites",
      "Maximum likelihood distance under the selected nucleotide model",
      "Tamura-Nei 93 distance (transition/transversion and base composition bias)",
      "Kimura 2-parameter distance",
      "Count of differing sites, uncorrected"
    )
  )
}

#' Divergence of each sequence from a consensus
#'
#' Appends the consensus to the alignment as an extra taxon, computes the
#' distance matrix under the chosen metric, and returns the consensus column.
#' No molecular clock correction is applied: the aim is regression on
#' environmental predictors rather than phylogenetic inference, and avoiding a
#' strict clock keeps the metric robust to rate variation across lineages.
#'
#' @section Why there is no clock correction here:
#' 0.1.0 offered `clock_correct = TRUE`, which divided the distance by
#' `pmax(year - min(year), 1)`. That is not a rate. The consensus is not an
#' ancestor, `min(year)` is not a root date, and the flooring at 1 gives the
#' two earliest cohorts the same denominator. On data simulated under a
#' perfectly constant rate the option produces a strong, overwhelmingly
#' significant *negative* association with year - it manufactures the artefact
#' its name promises to remove. Use [divergence_rate()], which roots a tree and
#' regresses root-to-tip distance on collection date.
#'
#' @param alignment A `vs_alignment` object.
#' @param consensus A `vs_consensus` object, or `NULL` to build one internally.
#' @param metric One of [distance_metrics()].
#' @param model Substitution model for the ML metrics. Accepts a
#'   `vs_model_selection` object, in which case its chosen model is used.
#' @param clock_correct Removed in 0.2.0. Setting it to `TRUE` is an error;
#'   use [divergence_rate()] instead. See the note below.
#' @param ... Passed to the underlying distance function.
#' @return A tibble with one row per sequence: `id`, `distance` and the metadata
#'   columns carried through from the alignment.
#' @export
distance_from_consensus <- function(alignment, consensus = NULL,
                                    metric = "ml_aa", model = NULL,
                                    clock_correct = FALSE, ...) {
  stopifnot(inherits(alignment, "vs_alignment"))
  metric <- match.arg(metric, distance_metrics()$metric)
  if (isTRUE(clock_correct)) {
    vs_abort(c(
      "{.arg clock_correct} was removed in 0.2.0 because it did not estimate a rate.",
      "x" = "Dividing divergence by elapsed years induces a negative year effect even under a constant rate.",
      "i" = "Use {.fn divergence_rate} for a root-to-tip rate estimate."
    ))
  }
  consensus <- consensus %||% build_consensus(alignment)
  if (inherits(model, "vs_model_selection")) model <- model$best
  model <- model %||% if (alignment$molecule == "protein") "JTT" else "GTR"

  check_metric_molecule(metric, alignment$molecule)

  mat <- rbind(alignment$matrix, `__consensus__` = consensus$sequence)
  d <- compute_distance(mat, metric, model, alignment$molecule, ...)
  dm <- as.matrix(d)
  dist_to_consensus <- dm["__consensus__", setdiff(rownames(dm), "__consensus__")]

  out <- tibble::tibble(
    id = names(dist_to_consensus),
    distance = as.numeric(dist_to_consensus)
  )
  out <- dplyr::left_join(out, alignment$meta, by = "id")

  attr(out, "metric") <- metric
  attr(out, "model") <- model
  class(out) <- c("vs_distances", class(out))
  out
}

check_metric_molecule <- function(metric, molecule) {
  spec <- distance_metrics()
  need <- spec$molecule[spec$metric == metric]
  if (need != "both" && need != molecule) {
    vs_abort(c(
      "Metric {.val {metric}} requires {need} sequences but the alignment is {molecule}.",
      "i" = "Retrieve nucleotide sequences with {.code molecule = \"nucleotide\"}, or pick an amino acid metric."
    ))
  }
}

compute_distance <- function(mat, metric, model, molecule, ...) {
  if (metric %in% c("ml_aa", "ml_nt")) {
    type <- if (molecule == "nucleotide") "DNA" else "AA"
    phy <- phangorn::phyDat(mat, type = type)
    if (type == "DNA") return(nt_ml_distance(phy, mat, model, ...))
    return(phangorn::dist.ml(phy, model = model, ...))
  }
  if (metric == "p_aa") {
    phy <- phangorn::phyDat(mat, type = "AA")
    return(phangorn::dist.hamming(phy, ratio = TRUE))
  }
  if (metric == "raw") {
    type <- if (molecule == "nucleotide") "DNA" else "AA"
    phy <- phangorn::phyDat(mat, type = type)
    return(phangorn::dist.hamming(phy, ratio = FALSE))
  }
  bin <- ape::as.DNAbin(tolower(mat))
  ape_model <- switch(metric, tn93 = "TN93", k80 = "K80")
  ape::dist.dna(bin, model = ape_model, pairwise.deletion = TRUE, ...)
}

#' Which distance a nucleotide substitution model maps onto
#'
#' Model *selection* and model-based *distances* do not support the same set of
#' models, and the gap is not obvious. `phangorn::optim.pml()` will happily fit
#' GTR, HKY, SYM and the rest, so they are legitimate entries in the ranking
#' table - but `phangorn::dist.ml()` implements only JC69 and F81 for DNA, and
#' its `bf` and `Q` arguments are ignored for nucleotides. There is no
#' closed-form GTR distance to fall back on either. So a selected model has to
#' be mapped onto the nearest model that a distance function actually
#' implements, and 0.4.0 and earlier simply passed the name straight through,
#' which failed with `'arg' should be one of "JC69", "F81", ...` for every
#' nucleotide model except those two.
#'
#' @section What "closest" means here:
#' Two axes could define it and they disagree, so the rule is stated rather than
#' left implicit: **the base-frequency assumption is preserved first, and rate
#' structure matched as closely as possible within that.**
#'
#' A model's frequency assumption changes the estimator systematically - an
#' equal-frequency model and an empirical-frequency one answer slightly
#' different questions - whereas extra rate classes mostly buy efficiency. So
#' `SYM` (six rates, equal base frequencies) maps to `K81` (three rates, equal
#' base frequencies) rather than to the richer `TN93`, which would silently
#' swap equal frequencies for empirical ones. `GTR` (six rates, empirical
#' frequencies) maps to `TN93` for the same reason in the other direction.
#'
#' @param model A model name, e.g. `"GTR"`.
#' @return A list with `engine` (`"phangorn"` or `"ape"`), `model` (the name
#'   that engine will accept), `exact` (whether the requested model is the one
#'   being used), `family` (`"equal"` or `"empirical"` base frequencies) and
#'   `requested`.
#' @examples
#' nt_distance_model("K81")   # exact: ape implements it
#' nt_distance_model("SYM")   # inexact, but keeps equal base frequencies
#' @export
nt_distance_model <- function(model) {
  m <- toupper(as.character(model)[1])
  ## engine, model, exact, base-frequency family
  map <- list(
    ## equal base frequencies
    JC     = list("phangorn", "JC69", TRUE,  "equal"),
    JC69   = list("phangorn", "JC69", TRUE,  "equal"),
    K80    = list("ape",      "K80",  TRUE,  "equal"),
    K2P    = list("ape",      "K80",  TRUE,  "equal"),
    K81    = list("ape",      "K81",  TRUE,  "equal"),
    K3P    = list("ape",      "K81",  TRUE,  "equal"),
    K3ST   = list("ape",      "K81",  TRUE,  "equal"),
    TPM1   = list("ape",      "K81",  TRUE,  "equal"),
    TPM2   = list("ape",      "K81",  FALSE, "equal"),
    TPM3   = list("ape",      "K81",  FALSE, "equal"),
    TIM1E  = list("ape",      "K81",  FALSE, "equal"),
    TIM2E  = list("ape",      "K81",  FALSE, "equal"),
    TIM3E  = list("ape",      "K81",  FALSE, "equal"),
    TVME   = list("ape",      "K81",  FALSE, "equal"),
    SYM    = list("ape",      "K81",  FALSE, "equal"),
    ## empirical base frequencies
    F81    = list("phangorn", "F81",  TRUE,  "empirical"),
    T92    = list("ape",      "T92",  TRUE,  "empirical"),
    F84    = list("ape",      "F84",  TRUE,  "empirical"),
    HKY    = list("ape",      "F84",  FALSE, "empirical"),
    TRN    = list("ape",      "TN93", TRUE,  "empirical"),
    TN93   = list("ape",      "TN93", TRUE,  "empirical"),
    TPM1U  = list("ape",      "TN93", FALSE, "empirical"),
    TPM2U  = list("ape",      "TN93", FALSE, "empirical"),
    TPM3U  = list("ape",      "TN93", FALSE, "empirical"),
    TIM1   = list("ape",      "TN93", FALSE, "empirical"),
    TIM2   = list("ape",      "TN93", FALSE, "empirical"),
    TIM3   = list("ape",      "TN93", FALSE, "empirical"),
    TIM    = list("ape",      "TN93", FALSE, "empirical"),
    TVM    = list("ape",      "TN93", FALSE, "empirical"),
    GTR    = list("ape",      "TN93", FALSE, "empirical")
  )
  hit <- map[[m]] %||% list("ape", "TN93", FALSE, "empirical")
  list(engine = hit[[1]], model = hit[[2]], exact = hit[[3]],
       family = hit[[4]], requested = m)
}

nt_ml_distance <- function(phy, mat, model, ...) {
  spec <- nt_distance_model(model)
  if (!spec$exact) {
    freq <- if (spec$family == "equal") "equal base frequencies" else "empirical base frequencies"
    vs_inform(c(
      "No distance is implemented for {.strong {spec$requested}}; using {.strong {spec$model}}, the closest one that keeps its {freq}.",
      "i" = "{.strong {spec$requested}} still governs the model ranking and any likelihood-based step; only the pairwise distance is substituted."
    ))
  }
  if (spec$engine == "phangorn") {
    return(phangorn::dist.ml(phy, model = spec$model, ...))
  }
  bin <- ape::as.DNAbin(tolower(mat))
  ape::dist.dna(bin, model = spec$model, pairwise.deletion = TRUE)
}

#' Pairwise distance matrix for an alignment
#'
#' @inheritParams distance_from_consensus
#' @return A `dist` object.
#' @export
pairwise_distances <- function(alignment, metric = "ml_aa", model = NULL, ...) {
  stopifnot(inherits(alignment, "vs_alignment"))
  if (inherits(model, "vs_model_selection")) model <- model$best
  model <- model %||% if (alignment$molecule == "protein") "JTT" else "GTR"
  check_metric_molecule(metric, alignment$molecule)
  compute_distance(alignment$matrix, metric, model, alignment$molecule, ...)
}

#' @export
print.vs_distances <- function(x, ...) {
  cli::cli_h3("vs_distances")
  cli::cli_text("{nrow(x)} sequence{?s}; metric {.strong {attr(x, 'metric')}} (model {attr(x, 'model')})")
  cli::cli_text("Distance range: {round(min(x$distance, na.rm = TRUE), 4)} - {round(max(x$distance, na.rm = TRUE), 4)}")
  NextMethod()
}
