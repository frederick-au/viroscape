#' Candidate substitution models
#'
#' @param molecule `"protein"` or `"nucleotide"`.
#' @return Character vector of model names understood by \pkg{phangorn}.
#' @export
candidate_models <- function(molecule = c("protein", "nucleotide")) {
  molecule <- match.arg(molecule)
  if (molecule == "protein") {
    c("JTT", "WAG", "LG", "Dayhoff", "Blosum62", "VT")
  } else {
    c("JC", "F81", "K80", "HKY", "SYM", "GTR")
  }
}

#' Select a substitution model by information criterion
#'
#' Fits each candidate model to the alignment on a fixed neighbour-joining tree
#' and ranks them. Holding the tree fixed across candidates is what makes the
#' comparison fair: only the substitution matrix varies. BIC is the default
#' because its penalty scales with alignment length, so it is more conservative
#' than AIC about selecting an unnecessarily parameterised matrix.
#'
#' @section Why the criterion rarely matters here:
#' The protein candidates - JTT, WAG, LG, Dayhoff, Blosum62, VT - are all fixed
#' empirical matrices. With the tree and branch lengths held fixed, the only
#' free parameter is the gamma shape, so every candidate has the same degrees of
#' freedom, the penalty term is a constant, and AIC, BIC and the raw
#' log-likelihood all produce the identical ranking. Switching `criterion`
#' cannot change which model is chosen unless `optimise_edges = TRUE` or the
#' candidate set mixes models with different parameter counts (as the nucleotide
#' set does). It is worth knowing before spending any effort defending the
#' choice.
#'
#' @param alignment A `vs_alignment` object.
#' @param models Candidate models. Defaults to [candidate_models()].
#' @param criterion `"BIC"` or `"AIC"`.
#' @param gamma Model rate heterogeneity across sites.
#' @param k Number of discrete gamma rate categories.
#' @param base_model Model used to build the fixed starting tree.
#' @param optimise_edges Re-optimise branch lengths for each candidate. `FALSE`
#'   (the default) holds them at the neighbour-joining estimates, which is much
#'   cheaper and keeps the comparison on identical footing.
#' @param empirical_freq Use amino acid / base frequencies observed in the data
#'   rather than the model's own assumed frequencies.
#' @return A `vs_model_selection` object with a ranked `table`, the chosen
#'   `best` model, and the fixed `tree`.
#' @export
select_substitution_model <- function(alignment, models = NULL,
                                      criterion = c("BIC", "AIC"),
                                      gamma = TRUE, k = 4,
                                      base_model = NULL,
                                      optimise_edges = FALSE,
                                      empirical_freq = TRUE) {
  stopifnot(inherits(alignment, "vs_alignment"))
  criterion <- match.arg(criterion)
  models <- models %||% candidate_models(alignment$molecule)
  base_model <- base_model %||% if (alignment$molecule == "protein") "LG" else "JC"

  phy <- alignment$phy
  vs_inform("Building fixed neighbour-joining tree under {.strong {base_model}}.")
  ## dist.ml names DNA models differently from optim.pml ("JC69", not "JC"), and
  ## implements only two of them; map before calling it
  d0 <- if (alignment$molecule == "nucleotide") {
    nt_ml_distance(phy, alignment$matrix, base_model)
  } else {
    phangorn::dist.ml(phy, model = base_model)
  }
  tree <- ape::nj(d0)
  tree <- ape::multi2di(tree)
  tree$edge.length <- pmax(tree$edge.length, 1e-8)

  n_sites <- attr(phy, "nr")
  rows <- lapply(models, function(m) {
    fit <- try(fit_one_model(phy, tree, m, gamma, k, optimise_edges,
                             empirical_freq, alignment$molecule), silent = TRUE)
    if (inherits(fit, "try-error")) {
      vs_warn("Model {.strong {m}} failed to fit and was skipped.")
      return(NULL)
    }
    ll <- as.numeric(stats::logLik(fit))
    npar <- attr(stats::logLik(fit), "df")
    tibble::tibble(
      model = m, logLik = ll, df = npar,
      AIC = -2 * ll + 2 * npar,
      BIC = -2 * ll + log(n_sites) * npar
    )
  })
  tab <- dplyr::bind_rows(rows)
  if (!nrow(tab)) vs_abort("No candidate substitution model could be fitted.")

  tab <- tab[order(tab[[criterion]]), ]
  tab$delta <- tab[[criterion]] - min(tab[[criterion]])
  tab$rank <- seq_len(nrow(tab))

  structure(
    list(table = tab, best = tab$model[1], criterion = criterion,
         tree = tree, gamma = gamma, k = k, n_sites = n_sites,
         molecule = alignment$molecule),
    class = "vs_model_selection"
  )
}

fit_one_model <- function(phy, tree, model, gamma, k, optimise_edges,
                          empirical_freq, molecule) {
  args <- list(tree = tree, data = phy)
  if (gamma) args$k <- k
  fit <- do.call(phangorn::pml, args)
  opt_args <- list(
    object = fit, model = model,
    optNni = FALSE, optEdge = optimise_edges,
    optGamma = gamma, optInv = FALSE,
    optBf = if (molecule == "nucleotide") empirical_freq else FALSE,
    control = phangorn::pml.control(trace = 0)
  )
  if (molecule == "protein" && empirical_freq) opt_args$optBf <- FALSE
  suppressWarnings(do.call(phangorn::optim.pml, opt_args))
}

#' @export
print.vs_model_selection <- function(x, ...) {
  cli::cli_h3("vs_model_selection")
  cli::cli_text("Best model by {x$criterion}: {.strong {x$best}} ({x$n_sites} sites)")
  print(as.data.frame(x$table), row.names = FALSE, digits = 6)
  if (length(unique(x$table$df)) == 1L) {
    cli::cli_text("{.emph Every candidate has {x$table$df[1]} parameters, so AIC, BIC and log-likelihood rank them identically; the criterion is not doing any work here.}")
  }
  invisible(x)
}
