#' Forward stepwise selection with collinearity and stability guards
#'
#' At each round every remaining candidate is added to the current model and
#' the models are compared on AIC. AIC rather than a likelihood ratio test is
#' used here because candidates within a round are not nested in one another,
#' and the likelihood ratio statistic only has its chi-squared null
#' distribution under nesting.
#'
#' The best candidate is then subjected to two further checks before it is
#' accepted, which is what stops the procedure chasing an AIC improvement that
#' comes from a predictor that is simply a re-expression of one already in the
#' model:
#'
#' * **Collinearity**: the scaled GVIF of every term is checked against
#'   `gvif_threshold`. By default an exceedance is recorded but does not by
#'   itself reject the candidate, because a high GVIF arising from structural
#'   correlation with a grouping factor is not the same problem as redundancy
#'   between two predictors. Set `gvif_action = "reject"` for the strict rule.
#' * **Stability**: among the terms you intend to interpret, no coefficient may
#'   change sign and lose significance, and no significant coefficient may shift
#'   by more than `stability_threshold` proportionally *and* by more than one
#'   standard error. By default *every base term* is excluded from this check,
#'   because base terms are what you are controlling for rather than what you
#'   are interpreting. A predictor that genuinely explains group differences
#'   collapses the grouping contrasts, and that is success; a real predictor
#'   correlated with time moves the `year` coefficient, and that is arithmetic.
#'   Movement in base terms is still reported - `base_max_shift` in the step log
#'   - it just does not veto. Set `stability_scope = "all"` for the strict
#'   reading, in which any shared coefficient can reject a candidate.
#'
#' Together these mean a candidate is rejected when it is both collinear with
#' something already present and visibly damaging to the estimates that matter.
#' Elevated GVIF on its own, arising from structural correlation with a
#' grouping factor, is recorded and reported rather than acted on.
#'
#' A candidate failing either check is rejected with a recorded reason and the
#' next best is tried. Selection stops when no remaining candidate improves AIC
#' by more than `delta_threshold`.
#'
#' @param dataset A `vs_dataset` from [join_environment()].
#' @param base_terms Terms always retained, e.g. `c("country", "year")`. If a
#'   `vs_structure` object is supplied, its choices are used.
#' @param candidates Candidate predictors. Defaults to all in the dataset.
#' @param random_term Optional random-effect term, e.g. `"(1|country)"`.
#' @param criterion `"AIC"` or `"BIC"`.
#' @param delta_threshold Improvement required to accept a predictor. The
#'   conventional -2 is the default.
#' @param gvif_threshold Scaled GVIF treated as severe collinearity.
#' @param gvif_action `"flag"` records an exceedance and lets the stability
#'   checks decide; `"reject"` rejects the candidate outright.
#' @param stability_threshold Proportional coefficient shift treated as
#'   destabilising.
#' @param stability_scope `"predictors"` exempts every base term from the
#'   stability check, judging only the coefficients you intend to interpret;
#'   `"all"` includes them, which is the stricter reading used in the published
#'   analysis and which will reject real predictors that happen to move a
#'   nuisance covariate.
#' @param allow_collinear_base If `TRUE`, a base term exceeding the GVIF
#'   threshold does not block selection. This matches the published analysis,
#'   where agricultural land's GVIF was inflated by structural correlation with
#'   country rather than by redundancy with another predictor.
#' @param allow_unstable_base Deprecated alias for `allow_collinear_base`. The
#'   old name was wrong: it never had anything to do with the stability check.
#' @param max_predictors Stop after this many environmental predictors.
#' @param calibrate Run the permutation calibration inline and use it. `TRUE`
#'   uses 199 permutations; a number sets the count; `FALSE` skips it. When it
#'   runs, `delta_threshold` is tightened to the calibrated 5% threshold if that
#'   is stricter, so the honest threshold is applied *during* the search rather
#'   than being available afterwards to whoever reads the documentation. The
#'   calibrated threshold is derived from the first round, where the search is
#'   over the most candidates, so applying it to later rounds is conservative.
#' @param calibration An existing `vs_calibration` from [calibrate_selection()].
#'   Supplying one skips the inline run and uses it instead.
#' @param verbose Print progress.
#' @return A `vs_selection` object.
#'
#' @section The selected coefficient is not an effect size:
#' The coefficient of a predictor that was chosen *because* it fitted well is
#' biased away from zero - the winner's curse. It is estimated on the same data
#' that selected it, so the sampling fluctuation that won it the round is inside
#' the estimate.
#'
#' How much this matters depends entirely on how marginal the effect is, and the
#' dependence is steep. Conditioning on the predictor actually being selected,
#' over twenty simulated datasets at each effect size:
#'
#' \preformatted{
#'   planted   selected   full-data estimate   held-out estimate
#'   0.0015      6/20       0.00340 (+127%)     0.00250 (+67%)
#'   0.0025     11/20       0.00393  (+57%)     0.00303 (+21%)
#' }
#'
#' A predictor that only sometimes wins carries an estimate inflated by more
#' than the effect itself; one that wins every time carries almost no inflation,
#' because nothing had to get lucky. [holdout_estimate()] reduced that inflation
#' in every regime measured here, but it does not remove it - these runs are
#' still conditioned on the predictor having been selected at all - and it
#' answers a slightly different question. Its documentation sets out which.
#'
#' This is inherent to selection, not a defect. Two things help, and they help
#' with different problems. [refit_at_cluster()] puts the estimate on the unit
#' at which the predictors vary, which fixes the standard errors and the degrees
#' of freedom - but note that when cluster sizes are equal and the predictors
#' are constant within a cluster, the weighted cluster-level coefficient is
#' *arithmetically identical* to the sequence-level one, so it does nothing at
#' all for the bias. Only [holdout_estimate()] addresses the bias, by selecting
#' on one half of the clusters and estimating on the other.
#'
#' @section What the deltas do and do not mean:
#' Every criterion here is computed on the rows of `dataset`. When those rows
#' are sequences and the predictors are country-year quantities, the rows are
#' not independent and the criterion is anti-conservative - severely so. On
#' simulated data with no environmental effect at all, this function retains at
#' least one pure-noise predictor in 100% of runs at the sequence level and
#' about 36% at the cluster level. Run [aggregate_to_cluster()] first, and
#' [calibrate_selection()] always.
#' @export
forward_select <- function(dataset, base_terms = NULL, candidates = NULL,
                           random_term = NULL, criterion = c("AIC", "BIC"),
                           delta_threshold = -2, gvif_threshold = 10,
                           gvif_action = c("flag", "reject"),
                           stability_threshold = 0.5,
                           stability_scope = c("predictors", "all"),
                           allow_collinear_base = TRUE,
                           allow_unstable_base = NULL,
                           max_predictors = Inf, calibrate = TRUE,
                           calibration = NULL, verbose = TRUE) {
  stopifnot(inherits(dataset, "vs_dataset"))
  criterion <- match.arg(criterion)
  gvif_action <- match.arg(gvif_action)
  stability_scope <- match.arg(stability_scope)
  if (!is.null(allow_unstable_base)) {
    vs_warn(c("{.arg allow_unstable_base} is deprecated; use {.arg allow_collinear_base}.",
              "i" = "It only ever controlled the collinearity check, never the stability one."))
    allow_collinear_base <- isTRUE(allow_unstable_base)
  }

  if (inherits(base_terms, "vs_structure")) {
    st <- base_terms
    base_terms <- st$base_terms
    random_term <- random_term %||% st$random_term
  }
  base_terms <- base_terms %||% intersect(c("country", "year"), names(dataset$data))
  candidates <- candidates %||% dataset$predictors
  candidates <- setdiff(candidates, base_terms)

  d <- dataset$data
  resp <- dataset$response
  wt <- model_weights(dataset)
  crit_fun <- if (criterion == "AIC") stats::AIC else stats::BIC

  if (identical(dataset$level %||% "sequence", "sequence") &&
      is.finite(dataset$design_effect %||% NA_real_) &&
      dataset$design_effect > 1.5) {
    vs_warn(c(
      "Selecting at the sequence level on data clustered by {.field {paste(dataset$cluster, collapse = ' x ')}} (design effect {round(dataset$design_effect, 1)}).",
      "x" = "{criterion} differences are inflated by roughly that factor; predictors will be retained that explain nothing.",
      "i" = "Run {.fn aggregate_to_cluster} first, and calibrate with {.fn calibrate_selection}."
    ))
  }

  ## Base terms are exempt from the stability veto unless the caller asks for
  ## the stricter scope. They are the covariates you decided in advance to keep,
  ## so movement in them is a consequence of adding a real predictor, not
  ## evidence against one - and `year` in particular will move whenever a
  ## candidate has any correlation with time, which is most of them.
  group_terms <- if (stability_scope == "all") character() else base_terms

  ## Calibrate first, so the threshold the search actually uses is the honest
  ## one. A warning that arrives after the step log has already been printed is
  ## a warning most readers will act on second, if at all.
  B <- if (isTRUE(calibrate)) 199L else if (is.numeric(calibrate)) as.integer(calibrate) else 0L
  if (is.null(calibration) && B > 0L && length(candidates)) {
    calibration <- try(calibrate_selection(
      dataset, base_terms = base_terms, candidates = candidates,
      random_term = random_term, criterion = criterion, B = B), silent = TRUE)
    if (inherits(calibration, "try-error")) {
      vs_warn(c("The permutation calibration failed; the search will use the nominal threshold.",
                "i" = "{conditionMessage(attr(calibration, 'condition'))}"))
      calibration <- NULL
    }
  }
  nominal_threshold <- delta_threshold
  if (!is.null(calibration) && is.finite(calibration$threshold_05) &&
      calibration$threshold_05 < delta_threshold) {
    delta_threshold <- calibration$threshold_05
    if (verbose) {
      cli::cli_alert_info("Calibration tightened the acceptance threshold from {nominal_threshold} to {round(delta_threshold, 1)}.")
    }
  }

  current <- fit_model(d, resp, base_terms, random_term, weights = wt)
  current_terms <- base_terms
  current_crit <- crit_fun(current)
  chain <- list(current)
  steps <- list()
  vetoed <- list()
  round_i <- 0L

  if (verbose) {
    cli::cli_h3("Forward selection")
    cli::cli_text("Base model {.code {deparse_terms(resp, current_terms, random_term)}}: {criterion} = {round(current_crit, 1)}")
  }

  remaining <- candidates
  while (length(remaining) && length(current_terms) - length(base_terms) < max_predictors) {
    round_i <- round_i + 1L
    trials <- lapply(remaining, function(p) {
      fit <- try(fit_model(d, resp, c(current_terms, p), random_term, weights = wt), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      tibble::tibble(predictor = p, criterion_value = crit_fun(fit),
                     delta = crit_fun(fit) - current_crit)
    })
    trials <- dplyr::bind_rows(trials)
    if (!nrow(trials)) break
    trials <- trials[order(trials$criterion_value), ]

    accepted <- NULL
    for (i in seq_len(nrow(trials))) {
      cand <- trials$predictor[i]
      delta <- trials$delta[i]

      if (delta > delta_threshold) {
        steps[[length(steps) + 1L]] <- step_row(
          round_i, cand, trials$criterion_value[i], delta, NA_real_, NA, NA_real_,
          FALSE, sprintf("no meaningful improvement (delta %s = %.1f)", criterion, delta)
        )
        if (verbose) {
          cli::cli_alert_info("Round {round_i}: best candidate {.field {cand}} gives delta{criterion} = {round(delta, 1)}; stopping.")
        }
        accepted <- NA
        break
      }

      fit <- fit_model(d, resp, c(current_terms, cand), random_term, weights = wt)
      gv <- gvif_table(fit)
      gv_new <- gv$gvif_scaled[gv$term == cand]
      offenders <- gv$term[gv$gvif_scaled > gvif_threshold]
      if (allow_collinear_base) offenders <- setdiff(offenders, base_terms)
      stab <- coef_stability(current, fit, threshold = stability_threshold,
                             group_terms = group_terms)
      flip <- isTRUE(attr(stab, "any_destabilising"))
      shift <- attr(stab, "max_prop_shift") %||% NA_real_
      ## exempt terms no longer veto, but their movement is still recorded
      base_shift <- suppressWarnings(max(stab$prop_shift[stab$is_group], na.rm = TRUE))
      if (!is.finite(base_shift)) base_shift <- NA_real_

      reason <- NA_character_
      if (length(offenders) && gvif_action == "reject") {
        reason <- sprintf("collinearity: scaled GVIF > %g for %s",
                          gvif_threshold, paste(offenders, collapse = ", "))
      } else if (flip) {
        culprits <- stab$term[stab$destabilising]
        reason <- sprintf("destabilised existing coefficients (%s)",
                          paste(culprits, collapse = ", "))
      }

      ok <- is.na(reason)
      note <- if (!ok) reason
        else if (length(offenders)) {
          sprintf("accepted; scaled GVIF > %g for %s but coefficients stayed stable",
                  gvif_threshold, paste(offenders, collapse = ", "))
        } else "accepted"

      steps[[length(steps) + 1L]] <- step_row(
        round_i, cand, trials$criterion_value[i], delta,
        if (length(gv_new)) gv_new else NA_real_,
        flip, shift, ok, note, base_shift
      )
      if (!ok) {
        vetoed[[length(vetoed) + 1L]] <- tibble::tibble(
          round = round_i, predictor = cand, delta = delta, reason = reason
        )
      }

      if (ok) {
        if (verbose) {
          cli::cli_alert_success("Round {round_i}: added {.field {cand}} (delta{criterion} = {round(delta, 1)}, scaled GVIF = {round(gv_new, 2)}).")
        }
        current <- fit
        current_terms <- c(current_terms, cand)
        current_crit <- trials$criterion_value[i]
        chain[[length(chain) + 1L]] <- fit
        accepted <- cand
        break
      } else if (verbose) {
        cli::cli_alert_warning("Round {round_i}: rejected {.field {cand}} - {reason}")
      }
    }

    if (is.null(accepted) || is.na(accepted)) break
    remaining <- setdiff(remaining, accepted)
  }

  steps_tbl <- dplyr::bind_rows(steps)
  names(steps_tbl)[names(steps_tbl) == "criterion_value"] <- tolower(criterion)
  vetoed_tbl <- dplyr::bind_rows(vetoed)

  ## A selection that kept nothing while its own calibration says the best
  ## candidate beat the permutation null is a contradiction, and it is the kind
  ## a reader will not notice unless it is stated.
  conflict <- NULL
  if (nrow(vetoed_tbl) && !length(setdiff(current_terms, base_terms)) &&
      !is.null(calibration) && isTRUE(calibration$p_value < 0.05)) {
    best <- vetoed_tbl[which.min(vetoed_tbl$delta), ]
    conflict <- list(
      predictor = best$predictor, delta = best$delta, reason = best$reason,
      calibration_p = calibration$p_value
    )
  }

  structure(
    list(model = current, formula = stats::formula(current),
         terms = current_terms, base_terms = base_terms,
         selected = setdiff(current_terms, base_terms),
         random_term = random_term, steps = steps_tbl, chain = chain,
         vetoed = vetoed_tbl, conflict = conflict,
         criterion = criterion, dataset = dataset, calibration = calibration,
         thresholds = list(delta = delta_threshold, delta_nominal = nominal_threshold,
                           gvif = gvif_threshold, stability = stability_threshold)),
    class = "vs_selection"
  )
}

step_row <- function(round_i, predictor, crit, delta, gvif, flip, shift, accepted,
                     reason, base_shift = NA_real_) {
  tibble::tibble(
    round = round_i, predictor = predictor, criterion_value = crit,
    delta = delta, gvif_scaled = gvif, sign_flip = flip,
    max_prop_shift = shift, base_max_shift = base_shift,
    accepted = accepted, reason = reason
  )
}

deparse_terms <- function(resp, terms, random_term) {
  rhs <- c(terms, random_term)
  if (!length(rhs)) rhs <- "1"
  paste(resp, "~", paste(rhs, collapse = " + "))
}

#' @export
print.vs_selection <- function(x, ...) {
  cli::cli_h3("vs_selection")
  cli::cli_text("Final model: {.code {deparse_terms(x$dataset$response, x$terms, x$random_term)}}")
  cli::cli_text("Selected predictor{?s}: {.val {if (length(x$selected)) x$selected else 'none'}}")
  fit <- model_fit_summary(x$model)
  cli::cli_text("Adjusted R-squared: {round(fit$adj_r_squared, 3)} | {x$criterion} = {round(fit$criterion_aic, 1)}")
  cat("\n")
  print(as.data.frame(x$steps), row.names = FALSE, digits = 4)
  cat("\n")
  if (is.null(x$calibration)) {
    cli::cli_alert_warning(c(
      "These {x$criterion} differences are uncalibrated."
    ))
    cli::cli_text("{.emph A search over several candidates produces improvements of this size on permuted data. Re-run with {.code calibrate = TRUE} before reporting any of them.}")
  } else {
    cal <- x$calibration
    cli::cli_text("Calibrated against {cal$B} cluster permutations: empirical p = {.strong {signif(cal$p_value, 3)}} (5% threshold delta {round(cal$threshold_05, 1)}).")
    if (isTRUE(x$thresholds$delta < x$thresholds$delta_nominal)) {
      cli::cli_text("Acceptance threshold tightened from {x$thresholds$delta_nominal} to {round(x$thresholds$delta, 1)} by that calibration.")
    }
  }
  if (!is.null(x$conflict)) {
    cn <- x$conflict
    cli::cli_alert_danger("This selection kept nothing, and that disagrees with its own calibration.")
    cli::cli_text("{.field {cn$predictor}} improved {x$criterion} by {round(cn$delta, 1)} and the permutation null puts that at p = {signif(cn$calibration_p, 3)}, but it was rejected: {cn$reason}.")
    cli::cli_text("{.emph Look at the rejection before concluding there is no effect. If the destabilised term is one you are controlling for rather than interpreting, add it to the base terms.}")
  }
  invisible(x)
}

#' Fit statistics for a selected model
#'
#' For mixed models, marginal and conditional pseudo R-squared are reported in
#' place of adjusted R-squared.
#'
#' @param model A fitted model.
#' @return A list of fit statistics.
#' @export
model_fit_summary <- function(model) {
  if (inherits(model, "merMod")) {
    r2 <- pseudo_r2(model)
    return(list(adj_r_squared = r2$conditional, r_squared = r2$marginal,
                marginal_r2 = r2$marginal, conditional_r2 = r2$conditional,
                criterion_aic = stats::AIC(model), criterion_bic = stats::BIC(model),
                n = stats::nobs(model), f_statistic = NA_real_, p_value = NA_real_))
  }
  s <- summary(model)
  fstat <- s$fstatistic
  p <- if (!is.null(fstat)) {
    stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  } else NA_real_
  list(
    adj_r_squared = s$adj.r.squared, r_squared = s$r.squared,
    marginal_r2 = NA_real_, conditional_r2 = NA_real_,
    criterion_aic = stats::AIC(model), criterion_bic = stats::BIC(model),
    n = stats::nobs(model),
    f_statistic = if (!is.null(fstat)) unname(fstat[1]) else NA_real_,
    df1 = if (!is.null(fstat)) unname(fstat[2]) else NA_real_,
    df2 = if (!is.null(fstat)) unname(fstat[3]) else NA_real_,
    p_value = unname(p)
  )
}

pseudo_r2 <- function(model) {
  var_f <- stats::var(as.vector(lme4::fixef(model) %*% t(lme4::getME(model, "X"))))
  vc <- as.data.frame(lme4::VarCorr(model))
  var_r <- sum(vc$vcov[is.na(vc$var2) & vc$grp != "Residual"])
  var_e <- vc$vcov[vc$grp == "Residual"]
  list(marginal = var_f / (var_f + var_r + var_e),
       conditional = (var_f + var_r) / (var_f + var_r + var_e))
}

#' Re-estimate a selected model at the cluster level
#'
#' Refits the selected terms on the country-year means, which is the unit at
#' which the predictors vary, so the standard errors and degrees of freedom
#' count independent observations rather than sequences.
#'
#' It does **not** address the winner's curse, and it is worth being precise
#' about why: with equal cluster sizes and predictors that are constant within a
#' cluster, the size-weighted cluster-level least squares coefficient is
#' algebraically the same number as the sequence-level one. The estimate moves
#' only when cluster sizes are uneven. For an estimate that the selection did
#' not touch, use [holdout_estimate()].
#'
#' `weight` is the argument that actually changes the number. Weighting each
#' country-year by the number of sequences behind it reproduces the
#' sequence-level estimate exactly; setting `weight = FALSE` gives every
#' country-year one vote, which is a different estimand and usually the one you
#' want when sampling effort is as uneven as it is in GenBank. It is also the
#' likelier explanation of a gap between a sequence-level coefficient and a
#' country-year one.
#'
#' @param selection A `vs_selection`.
#' @param dataset The sequence-level `vs_dataset`, if the selection was run on
#'   an aggregated one and you want it refitted from the original rows.
#' @param weight Weight each cluster by the number of sequences behind it.
#'   `TRUE` reproduces the sequence-level coefficient and changes only the
#'   standard errors; `FALSE` gives each country-year equal weight.
#' @return A tibble from [tidy_coefficients()] for the refitted model, with the
#'   selected model attached as an attribute.
#'
#'   Its `unit` column is **mixed**: standardised predictors read `"standard
#'   deviation"` while `year` and the grouping contrasts read `"natural unit"`.
#'   Comparing two such tibbles row by row without looking at that column is an
#'   easy way to conclude that two identical estimates differ. [unscale()]
#'   handles the mixed case correctly - it converts only the standardised rows -
#'   but it should be applied to the whole tibble, never to a hand-picked row.
#' @examples
#' \donttest{
#' a <- example_analysis()
#' refit_at_cluster(a$selection, weight = FALSE)
#' }
#' @export
refit_at_cluster <- function(selection, dataset = NULL, weight = TRUE) {
  stopifnot(inherits(selection, "vs_selection"))
  ds <- dataset %||% selection$dataset
  if (!identical(ds$level, "cluster") || !weight) {
    base <- if (identical(ds$level, "cluster")) ds else ds
    ds <- aggregate_to_cluster(base, weight = weight)
  }
  fit <- fit_model(ds$data, ds$response, selection$terms, selection$random_term,
                   weights = model_weights(ds))
  out <- tidy_coefficients(fit, dataset = ds)
  attr(out, "model") <- fit
  attr(out, "note") <- "Refitted at the cluster level; the predictor set is still the one the search chose."
  out
}

#' Select on one half of the clusters, estimate on the other
#'
#' The coefficient of a predictor chosen because it fitted well is biased away
#' from zero, and no refit on the same data removes that. The standard remedy is
#' to spend the data twice over deliberately: split the clusters, let the search
#' see only the first half, and estimate the magnitude on the second half, which
#' had no say in which predictor won.
#'
#' The cost is power - half the clusters select and half estimate - so this is
#' the estimate to report once the search has already been calibrated, not a
#' replacement for the search. With `repeats > 1` the split is repeated and the
#' estimates averaged, which reduces the dependence on one arbitrary partition.
#'
#' @section It answers a different question from the full-data estimate:
#' This is the distinction to hold on to, because the two numbers are not two
#' attempts at the same quantity.
#'
#' The full-data coefficient answers *"the search picked this term; how big is
#' it?"* - and it is inflated, because the term was picked partly for being
#' large. The held-out estimate answers *"for terms a half-sample selects, what
#' does independent data say?"* - it is estimated on clusters that had no part
#' in choosing the term, but the term set is the split's, not the one you were
#' handed by the full-data run.
#'
#' Measured, conditioning on the full-data search having selected the predictor,
#' at a planted effect of 0.0015 under two noise levels:
#'
#' \preformatted{
#'   selected   full-data estimate      held-out estimate
#'    12/60     +136%  RMSE 2.2e-3      +66%  RMSE 1.4e-3
#'    22/60      +61%  RMSE 1.1e-3      +11%  RMSE 1.0e-3
#' }
#'
#' In both, the held-out estimate carried less bias and no more error. But it is
#' fitted on half the clusters, so it is the noisier estimator run for run, and
#' a careful reviewer measuring a third regime found it the worse of the two -
#' so treat the ordering as an empirical question about your data rather than a
#' rule. Both numbers are returned; compare them.
#'
#' `times_selected` is what tells you how much room there is for a correction at
#' all. Short of `n_splits` means selection was uncertain and the full-data
#' estimate is carrying the fluctuation that won the round. Equal to it means
#' the predictor wins on every split, there is no curse to remove, and the
#' held-out estimate is simply the noisier of the two.
#'
#' A term selected by *no* split does not appear in the table at all. That
#' filter is correlated with the estimate, so an absent term is a result -
#' nothing survived a half-sample - and not a missing value.
#'
#' @param dataset A `vs_dataset`.
#' @param base_terms,candidates,random_term,criterion Passed to
#'   [forward_select()].
#' @param prop Proportion of clusters used for selection.
#' @param repeats Number of independent splits to average over.
#' @param seed Optional seed.
#' @param ... Passed to [forward_select()].
#' @return A tibble with one row per term that at least one split selected:
#'   `times_selected` out of `n_splits`, the mean held-out estimate, the mean
#'   estimate from the selection half for comparison, and `unit`/`sd` saying
#'   what one unit of the predictor is - standardised, if [join_environment()]
#'   standardised it. Terms no split selected are absent rather than `NA`.
#' @examples
#' \donttest{
#' ds <- aggregate_to_cluster(example_analysis()$dataset)
#' holdout_estimate(ds, base_terms = c("country", "year"), repeats = 5)
#' }
#' @export
holdout_estimate <- function(dataset, base_terms = NULL, candidates = NULL,
                             random_term = NULL, criterion = c("AIC", "BIC"),
                             prop = 0.5, repeats = 10, seed = NULL, ...) {
  stopifnot(inherits(dataset, "vs_dataset"))
  criterion <- match.arg(criterion)
  if (inherits(base_terms, "vs_structure")) {
    random_term <- random_term %||% base_terms$random_term
    base_terms <- base_terms$base_terms
  }
  base_terms <- base_terms %||% intersect(c("country", "year"), names(dataset$data))
  candidates <- setdiff(candidates %||% dataset$predictors, base_terms)
  if (!is.null(seed)) set.seed(seed)

  d <- dataset$data
  key <- if (".cluster" %in% names(d)) as.character(d$.cluster) else as.character(seq_len(nrow(d)))
  lv <- unique(key)
  if (length(lv) < 8L) vs_abort("Splitting needs at least eight clusters; this dataset has {length(lv)}.")

  rows <- list()
  for (r in seq_len(repeats)) {
    pick <- sample(lv, max(2L, floor(prop * length(lv))))
    sel_ds <- est_ds <- dataset
    sel_ds$data <- d[key %in% pick, , drop = FALSE]
    est_ds$data <- d[!key %in% pick, , drop = FALSE]
    sel <- try(forward_select(sel_ds, base_terms = base_terms, candidates = candidates,
                              random_term = random_term, criterion = criterion,
                              calibrate = FALSE, verbose = FALSE, ...), silent = TRUE)
    if (inherits(sel, "try-error") || !length(sel$selected)) next
    held <- try(fit_model(est_ds$data, dataset$response, sel$terms, random_term,
                          weights = model_weights(est_ds)), silent = TRUE)
    if (inherits(held, "try-error")) next
    ch <- safe_coef(held)
    cs <- safe_coef(sel$model)
    for (tm in sel$selected) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        term = tm,
        estimate_selection_half = unname(cs[tm]) %||% NA_real_,
        estimate_heldout = unname(ch[tm]) %||% NA_real_
      )
    }
  }
  if (!length(rows)) {
    vs_warn("No split selected any predictor; nothing to estimate.")
    return(tibble::tibble(term = character(), times_selected = integer(),
                          n_splits = integer(), estimate_selection_half = numeric(),
                          estimate_heldout = numeric(), inflation = numeric(),
                          unit = character(), sd = numeric()))
  }
  all_rows <- dplyr::bind_rows(rows)
  out <- dplyr::summarise(
    dplyr::group_by(all_rows, term),
    times_selected = dplyr::n(),
    estimate_selection_half = mean(.data$estimate_selection_half, na.rm = TRUE),
    estimate_heldout = mean(.data$estimate_heldout, na.rm = TRUE),
    .groups = "drop"
  )
  out$n_splits <- repeats
  out$inflation <- out$estimate_selection_half / out$estimate_heldout
  ## carry the units, as tidy_coefficients() and refit_at_cluster() do; without
  ## them these estimates are silently in standard deviations
  sds <- scaling_sd(dataset, out$term)
  out$unit <- unname(ifelse(is.na(sds), "natural unit", "standard deviation"))
  out$sd <- unname(sds)
  out <- out[order(-out$times_selected), ]
  attr(out, "repeats") <- repeats
  attr(out, "prop") <- prop
  out
}
