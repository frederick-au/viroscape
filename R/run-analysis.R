#' Run the whole pipeline end to end
#'
#' A single call covering retrieval through to a validated final model. Every
#' stage is also available individually if you want to intervene between steps.
#'
#' @param sequences A `vs_sequences` object, or a path to a FASTA file, or
#'   `NULL` to fetch from NCBI using `countries` and `years`.
#' @param environment A country-year tibble, a path to a CSV, or `NULL` to fetch
#'   from the World Bank.
#' @param countries,years Used when retrieving live data.
#' @param molecule `"protein"` or `"nucleotide"`.
#' @param metric Divergence metric; see [distance_metrics()].
#' @param aligner Alignment backend; see [available_aligners()].
#' @param hosts Optional host classes to retain.
#' @param variables Environmental indicators to consider.
#' @param criterion Information criterion for model selection.
#' @param gvif_threshold,stability_threshold,delta_threshold Guards passed to
#'   [forward_select()].
#' @param cache_dir,refresh Caching controls.
#' @param subtype Subtype guard applied after retrieval; see
#'   [filter_sequences()].
#' @param level `"cluster"` runs the model selection at the country-year level,
#'   which is where the environmental predictors vary and therefore the only
#'   level at which the information criteria count independent observations.
#'   `"sequence"` restores the 0.1.0 behaviour, which treats every sequence as
#'   an independent draw and retains pure-noise predictors in essentially every
#'   run; see [aggregate_to_cluster()].
#' @param calibrate Calibrate the selection against a cluster-level permutation
#'   null; see [calibrate_selection()].
#' @param B Number of permutations when `calibrate = TRUE`.
#' @param verbose Print progress.
#' @return A `vs_analysis` object holding every intermediate result.
#' @export
run_analysis <- function(sequences = NULL, environment = NULL,
                         countries = names(sea_countries()),
                         years = c(2003, 2022),
                         molecule = c("protein", "nucleotide"),
                         metric = NULL, aligner = "auto", hosts = NULL,
                         variables = NULL, criterion = c("AIC", "BIC"),
                         gvif_threshold = 10, stability_threshold = 0.5,
                         delta_threshold = -2, subtype = NULL,
                         level = c("cluster", "sequence"),
                         calibrate = TRUE, B = 499,
                         cache_dir = vs_cache_dir(), refresh = FALSE,
                         verbose = TRUE) {
  molecule <- match.arg(molecule)
  criterion <- match.arg(criterion)
  level <- match.arg(level)
  metric <- metric %||% if (molecule == "protein") "ml_aa" else "ml_nt"

  seqs <- resolve_sequences(sequences, countries, years, molecule, cache_dir, refresh)
  if (!is.null(hosts)) seqs <- filter_sequences(seqs, hosts = hosts)
  seqs <- filter_sequences(seqs, countries = countries, years = years,
                           subtype = subtype)

  if (verbose) print(seqs)

  aln <- align_sequences(seqs, method = aligner)
  msel <- select_substitution_model(aln, criterion = "BIC")
  if (verbose) print(msel)

  cons <- build_consensus(aln)
  dists <- distance_from_consensus(aln, cons, metric = metric, model = msel)

  env <- resolve_environment(environment, countries, years, variables, cache_dir, refresh)
  ds <- join_environment(dists, env)
  if (verbose) print(ds)

  model_ds <- if (level == "cluster") aggregate_to_cluster(ds) else ds
  if (level == "cluster" && verbose) {
    vs_inform("Modelling at the country-year level: {nrow(model_ds$data)} row{?s} from {nrow(ds$data)} sequence{?s}.")
  }

  struct <- select_structure(model_ds, criterion = criterion)
  if (verbose) print(struct)

  cal <- if (calibrate) {
    calibrate_selection(model_ds, base_terms = struct, criterion = criterion, B = B)
  } else NULL

  sel <- forward_select(model_ds, base_terms = struct, criterion = criterion,
                        gvif_threshold = gvif_threshold,
                        stability_threshold = stability_threshold,
                        delta_threshold = delta_threshold,
                        calibrate = calibrate, calibration = cal,
                        verbose = verbose)

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

resolve_sequences <- function(sequences, countries, years, molecule,
                              cache_dir, refresh) {
  if (inherits(sequences, "vs_sequences")) return(sequences)
  if (is.character(sequences) && length(sequences) == 1L) {
    return(read_sequences(sequences, molecule = molecule))
  }
  if (is.null(sequences)) {
    return(fetch_sequences(countries, years, molecule = molecule,
                           cache_dir = cache_dir, refresh = refresh))
  }
  vs_abort("{.arg sequences} must be a vs_sequences object, a file path, or NULL.")
}

resolve_environment <- function(environment, countries, years, variables,
                                cache_dir, refresh) {
  if (is.data.frame(environment)) return(tibble::as_tibble(environment))
  if (is.character(environment) && length(environment) == 1L) {
    return(read_environment(environment))
  }
  if (is.null(environment)) {
    return(fetch_environment(countries, years, variables,
                             cache_dir = cache_dir, refresh = refresh))
  }
  vs_abort("{.arg environment} must be a data frame, a file path, or NULL.")
}

#' @export
print.vs_analysis <- function(x, ...) {
  cli::cli_h2("viroscape analysis")
  cli::cli_text("{nrow(x$dataset$data)} sequence{?s} across {nlevels(x$dataset$data$country)} countr{?y/ies}")
  cli::cli_text("Substitution model: {.strong {x$model_selection$best}} | metric: {.strong {attr(x$distances, 'metric')}}")
  cli::cli_text("Selected predictor{?s}: {.val {if (length(x$selection$selected)) x$selection$selected else 'none'}}")
  cli::cli_text("Adjusted R-squared: {round(x$fit$adj_r_squared, 3)}")
  cat("\n")
  cli::cli_h3("Coefficients")
  print(as.data.frame(x$coefficients), row.names = FALSE, digits = 3)
  invisible(x)
}

#' Write a plain-text report of an analysis
#'
#' @param x A `vs_analysis` object.
#' @param file Output path. `NULL` returns the lines instead of writing.
#' @return The report lines, invisibly when written to a file.
#' @export
report_analysis <- function(x, file = NULL) {
  stopifnot(inherits(x, "vs_analysis"))
  fmt <- function(d) paste(utils::capture.output(print(as.data.frame(d), row.names = FALSE, digits = 4)), collapse = "\n")
  banner <- if (is_example_analysis(x)) c(
    "*******************************************************************",
    "*  SIMULATED DATA. The bundled example sequences are generated    *",
    "*  from the environmental predictors with planted effects         *",
    "*  (agricultural land positive, livestock index negative). Any    *",
    "*  result below recovers that construction and says nothing about *",
    "*  influenza. See data-raw/make-example-data.R.                   *",
    "*******************************************************************",
    ""
  ) else character()

  cs <- cluster_summary(x$dataset)
  caveats <- c(
    "",
    "How to read this report",
    sprintf("  Modelling level: %s (%d rows from %d sequences)",
            x$level %||% "sequence",
            nrow((x$model_dataset %||% x$dataset)$data), nrow(x$dataset$data)),
    sprintf("  Clustering: %d country-years, ICC %.3f, design effect %.1f",
            cs$n_clusters %||% NA_integer_, cs$icc %||% NA_real_,
            cs$design_effect %||% NA_real_),
    if (isTRUE(x$dataset$scaled))
      "  Coefficients are per standard deviation of the predictor, not per natural unit (see unscale())."
    else "  Coefficients are per natural unit.",
    "  Likelihood ratio p-values are post-selection and anti-conservative; they are labelled nominal_p.",
    "  Coefficients below are post-selection and biased away from zero (the winner's curse);",
    "  use refit_at_cluster() for a magnitude that was not chosen for its size.",
    if (is.null(x$calibration))
      "  The selection was NOT calibrated. Do not interpret the deltas."
    else sprintf("  Permutation calibration: empirical p = %.3g over %d permutations.",
                 x$calibration$p_value, x$calibration$B),
    ""
  )

  lines <- c(
    banner,
    "viroscape analysis report",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    caveats,
    sprintf("Sequences: %d across %d countries (%d-%d)",
            nrow(x$dataset$data), nlevels(x$dataset$data$country),
            min(x$dataset$data$year), max(x$dataset$data$year)),
    sprintf("Alignment: %d sites, backend %s", x$alignment$n_sites, x$alignment$method),
    sprintf("Substitution model: %s (selected by %s)", x$model_selection$best, x$model_selection$criterion),
    sprintf("Divergence metric: %s", attr(x$distances, "metric")),
    "",
    "Substitution model ranking", fmt(x$model_selection$table), "",
    "Sampling by country", fmt(x$sampling$by_country), "",
    "Structural comparison", fmt(x$structure$table), "",
    "Forward selection steps", fmt(x$selection$steps), "",
    "Sequential likelihood ratio tests", fmt(x$lrt), "",
    "Collinearity (scaled GVIF)", fmt(x$gvif), "",
    "Final coefficients", fmt(x$coefficients), "",
    sprintf("Adjusted R-squared: %.3f", x$fit$adj_r_squared),
    sprintf("Observations: %d", x$fit$n),
    if (!is.null(x$calibration))
      paste0("", utils::capture.output(print(x$calibration))) else NULL
  )
  lines <- lines[!vapply(lines, is.null, logical(1))]
  if (is.null(file)) return(lines)
  writeLines(lines, file)
  invisible(lines)
}

## Provenance: was this built from the bundled, simulated example?
is_example_analysis <- function(x) {
  src <- x$sequences$source %||% ""
  isTRUE(grepl("example_sequences", src, fixed = TRUE))
}
