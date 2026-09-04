#' Bundled example data
#'
#' A small synthetic H5-like haemagglutinin dataset paired with real World Bank
#' indicators for the five Southeast Asian countries in the source study. It
#' exists so that the pipeline, the tests and the Shiny app all run with no
#' network access and no Bioconductor dependency. The sequences are simulated,
#' not real isolates: use them to learn the interface, not to draw biological
#' conclusions.
#'
#' @return `example_sequences()` returns a `vs_sequences` object;
#'   `example_environment()` returns a country-year tibble.
#' @export
example_sequences <- function() {
  path <- system.file("extdata", "example_sequences.fasta", package = "viroscape")
  if (!nzchar(path)) vs_abort("Example sequences are not installed.")
  read_sequences(path, molecule = "protein")
}

#' @rdname example_sequences
#' @export
example_environment <- function() {
  path <- system.file("extdata", "example_environment.csv", package = "viroscape")
  if (!nzchar(path)) vs_abort("Example environment data is not installed.")
  read_environment(path)
}

#' Run the pipeline on the bundled example data
#'
#' The sequences are simulated from the environmental predictors with planted
#' effects, so this recovers a construction rather than a finding. Reports built
#' from it carry a banner saying so.
#'
#' @param metric Divergence metric.
#' @param criterion Information criterion.
#' @param level `"cluster"` (default) or `"sequence"`; see [run_analysis()].
#' @param calibrate Calibrate the selection against a permutation null.
#' @param B Permutations used when `calibrate = TRUE`.
#' @param verbose Print progress.
#' @param ... Passed to [forward_select()].
#' @return A `vs_analysis` object.
#' @export
example_analysis <- function(metric = "ml_aa", criterion = "AIC",
                             level = c("cluster", "sequence"),
                             calibrate = FALSE, B = 199,
                             verbose = FALSE, ...) {
  level <- match.arg(level)
  seqs <- example_sequences()
  aln <- align_sequences(seqs, method = "none")
  msel <- select_substitution_model(aln, criterion = "BIC")
  cons <- build_consensus(aln)
  dists <- distance_from_consensus(aln, cons, metric = metric, model = msel)
  env <- example_environment()
  ds <- join_environment(dists, env)
  model_ds <- if (level == "cluster") aggregate_to_cluster(ds) else ds
  struct <- select_structure(model_ds, criterion = criterion)
  cal <- if (calibrate) {
    calibrate_selection(model_ds, base_terms = struct, criterion = criterion, B = B)
  } else NULL
  sel <- forward_select(model_ds, base_terms = struct, criterion = criterion,
                        calibrate = calibrate, calibration = cal,
                        verbose = verbose, ...)
  structure(
    list(sequences = seqs, alignment = aln, model_selection = msel,
         consensus = cons, distances = dists, environment = env,
         dataset = ds, model_dataset = model_ds, level = level,
         structure = struct, selection = sel, calibration = cal,
         lrt = lrt_chain(sel), coefficients = tidy_coefficients(sel),
         gvif = gvif_table(sel$model), fit = model_fit_summary(sel$model),
         sampling = sampling_summary(ds)),
    class = "vs_analysis"
  )
}
