# viroscape 0.5.2

Documentation only, from the 0.5.1 review, plus one column.

* `holdout_estimate()` now says what quantity it estimates and how that differs
  from the full-data coefficient. They are not two attempts at one number: the
  full-data estimate answers "the search picked this term, how big is it?", the
  held-out estimate answers "for terms a half-sample selects, what does
  independent data say?". Which has the lower error on a given dataset is an
  empirical question, and the function returns both so it can be checked.
* The claim that it "roughly halves" the winner's curse is gone. It reduced the
  inflation in every regime measured here (+136% to +66%, and +61% to +11%,
  conditioning on the full-data search having selected the predictor), but a
  careful reviewer measured a third regime where it was the worse of the two, so
  the docs now report the measurements and decline to state a rule.
* Documented that a term selected by no split is absent from the table rather
  than `NA`. That filter is correlated with the estimate, so absence is a
  result, not a missing value.
* New `n_splits` column, so `times_selected` can be read without going back to
  the call.

# viroscape 0.5.1

Documentation and mapping fixes from the 0.5.0 review. No behaviour changes
outside the nucleotide model mapping.

* `nt_distance_model()` no longer approximates models `ape::dist.dna()`
  implements exactly: `K81`/`K3P` and `T92` now map to themselves.
* The mapping rule is stated rather than left implicit, and one substitution
  changed because of it. **The base-frequency assumption is preserved first,
  and rate structure matched as closely as possible within that.** `SYM` has six
  rates and equal base frequencies; it previously mapped to `TN93`, which is
  richer in rates but estimates frequencies from the data. It now maps to `K81`,
  which keeps equal frequencies. The returned list carries a `family` field
  (`"equal"` or `"empirical"`) and the substitution message names which
  assumption was preserved.
* `holdout_estimate()` returns `unit` and `sd`, as `tidy_coefficients()` and
  `refit_at_cluster()` do. Without them its estimates were silently in standard
  deviations.
* `holdout_estimate()` documents when *not* to use it: `times_selected` is the
  diagnostic, and when it equals `repeats` the full-data estimate is the better
  one, because there is no curse to remove and the holdout only adds variance.
* `refit_at_cluster()` documents that its `unit` column is mixed - standardised
  predictors and natural-unit terms in one tibble - so `unscale()` is applied to
  the whole table, never to a hand-picked row.
* Corrected the winner's-curse figure in `forward_select()`. The previous "27%
  high" came from a single draw and did not survive replication. Replaced with a
  measured profile across effect sizes, which shows the inflation is large when
  selection is uncertain (+127% at an effect selected 6 times in 20) and small
  when it is not, and that `holdout_estimate()` roughly halves it rather than
  removing it.

# viroscape 0.5.0

Fixes from the 0.4.0 review, plus a defect that review turned up: the entire
nucleotide maximum-likelihood path was broken.

## Nucleotide distances actually work now

`distance_from_consensus(metric = "ml_nt")` failed with `'arg' should be one of
"JC69", "F81", "WAG", ...` for every nucleotide model except two, and
`select_substitution_model()` failed on its very first call for nucleotide data,
because its default `base_model = "JC"` is not a name `phangorn::dist.ml()`
accepts either.

The cause is a mismatch nothing in the package accounted for:
`phangorn::optim.pml()` fits GTR, HKY, SYM and the rest, so they are legitimate
entries in the model ranking, but `dist.ml()` implements only JC69 and F81 for
DNA - and its `bf` and `Q` arguments are ignored for nucleotides, so a GTR
distance cannot be recovered that way. There is no closed-form GTR distance
either.

* New `nt_distance_model()` maps a selected model onto the nearest distance that
  is actually implemented, routing to `ape::dist.dna()` where that is richer
  (K80, F84, TN93) and reporting the substitution rather than making it
  silently.
* `select_substitution_model()` maps its base model before building the
  starting tree.
* New `test-nucleotide.R`. The package had no nucleotide test at all, which is
  why a whole molecule type could be broken without failing anything.

## Post-selection estimates

* New `holdout_estimate()`: select on one half of the clusters, estimate on the
  other. This is the only thing here that addresses the winner's curse.
* New `refit_at_cluster()`, with an honest account of what it does. Its `weight`
  argument is the one that matters: size-weighted cluster means reproduce the
  sequence-level coefficient *exactly*, so `weight = TRUE` changes only the
  standard errors, while `weight = FALSE` gives each country-year one vote and
  is the likelier explanation of a gap between a sequence-level coefficient and
  a country-year one.
* `forward_select()` documents the winner's curse, and `report_analysis()` says
  in the report that the coefficients are post-selection.

## Smaller

* The `phylo_gls()` print header shows shrinkage to one decimal. 99.6% and 100%
  are different claims and the second reads as complete elimination.
* `fetch_environment()` retries a failed World Bank request three times with
  backoff, and when it does give up names the likely causes - proxy, firewall,
  timeout - the current timeout setting, and `read_environment()` as the offline
  route, instead of reporting only "cannot open the connection".

# viroscape 0.4.0

Third methods-review release. Both findings from the 0.3.0 review are fixed,
along with a Shiny layout regression introduced in 0.3.0.

## The stability guard no longer vetoes real effects

`forward_select()` was discarding a correctly specified predictor with a
calibrated, overwhelming effect because adding it moved the `year` coefficient -
while the calibration inside the very same object reported p = 0.005. Two rules
that are individually right were pulling against each other: `force_time = TRUE`
keeps `year` as a nuisance control, and the stability guard then treated
movement in it as evidence against the candidate.

* `stability_scope = "predictors"` (the default) now exempts **every** base
  term, not only the grouping factors. Base terms are what you control for, not
  what you interpret, and any real predictor correlated with time will move
  `year`. `stability_scope = "all"` restores the strict reading.
* Movement in exempt terms is still reported, as `base_max_shift` in the step
  log. The veto narrowed; the reporting did not.
* Detection of a planted country-year effect went from 1 in 30 runs to 10/10 at
  effect-to-noise 2.5 and 5, and 3/10 at 1, where a low rate is honest. The
  false-positive rate on pure noise is 1 in 25 with the offending calibration
  p = 0.035 - correctly sized against a nominal 5%, where the previous 0 in 25
  was the stability guard doing double duty.
* `allow_unstable_base` is deprecated in favour of `allow_collinear_base`. The
  old name was wrong: it only ever controlled the collinearity check.
* A selection that keeps nothing while its own calibration says the best
  candidate beat the permutation null now records that contradiction in
  `$conflict` and states it at `print()`, instead of quietly reporting nothing.

## phylo_gls() no longer over-claims

* The comparison table gains `estimate_clade`, `p_clade` and `shrinkage_clade`
  from a third fit with an explicit clade factor. A smooth Brownian or Pagel
  covariance approximates discrete block structure rather than absorbing it, so
  it can under-correct a lineage confound and still report significance; a
  categorical clade term does not. `print()` warns when the two disagree.
* `survives` is renamed `survives_pgls` and `shrinkage` to `shrinkage_pgls`, and
  the largest shrinkage is promoted into the print header, because it is the
  informative number.
* Documentation and the vignette now lead with the clade factor and treat
  `phylo_gls()` as the secondary check.

## Shiny layout

* Fixed the plots being squeezed into unreadably small panels, and the "invalid
  quartz() size" error before it. The panels scroll rather than fill, cards
  stack full width instead of splitting a narrow main column in two, plot
  heights are larger, and both sidebars are narrower. Measured at a 1440px
  viewport, the divergence plots went from 310x300 to 758x420 and the sampling
  plot to 1068x460.

# viroscape 0.3.0

Second methods-review release. 0.2.0 fixed the four defects that changed
results; this closes the last gap in that list and tightens the small edges
around it.

## Phylogenetic structure in the regression

* New `build_tree()`, `assign_clades()` and `add_clades()`. Clades are cut from
  the tree at a depth from the root, so groups are monophyletic by
  construction, and can be used as a model term, a random effect, or the
  clustering variable.
* New `phylo_gls()` fits the regression under a Brownian or Pagel covariance
  derived from the tree and reports every coefficient beside its ordinary least
  squares counterpart. On simulated data where an environmental predictor is
  associated with divergence only because lineages sit in particular countries,
  the association is significant at p = 6e-33 without the tree and p = 0.08
  with it, shrinking by 91%. Cluster-robust standard errors do not catch this;
  nothing before this release did.

## Calibration is now the default path, not the diligent one

* `forward_select(calibrate = TRUE)` is the default. The permutation null runs
  before the search and tightens `delta_threshold` to the calibrated 5%
  threshold when that is stricter, so the honest threshold is applied *during*
  selection instead of being available afterwards to whoever reads the warning.
  `calibrate = FALSE` restores the previous behaviour.
* `calibrate_selection()` gained `strata`, defaulting to the grouping terms
  already in `base_terms`. Permuting country-years *within* country holds the
  base model fixed and matches the question being asked; permuting across
  countries answered an easier one and gave a null that was too generous (10%
  false positives against a nominal 5%; within-country gives 3.6% with a
  uniform p-value distribution).

## Rates

* `divergence_rate(root = "best")` is the new default: the root position that
  minimises the residual variance of the root-to-tip regression, rather than
  the midpoint. On simulated constant-rate data with an outgroup present it
  recovers the realised rate to -3% against midpoint's -7%.
* `root = "outgroup"` added, with the outgroup dropped after rooting.
* A second guard warns when the regression explains less than 10% of the
  variance. Without a deep split to anchor it, root-to-tip distance is tip
  noise, and the function now says so instead of returning a slope.

## Host classification

* A bare genus - `Anas sp.`, `Anser sp.`, `Anatidae`, an unqualified
  `waterfowl` - is now `unknown` rather than `domestic_poultry`. It identifies a
  bird only as "some duck", which is less information than a vernacular `duck`,
  so it should not receive a more confident answer than the specific label.
* Fixed: the environment rule matched `water` inside `waterfowl`, so wild
  waterfowl were classified as environmental samples.

## Smaller things

* `detect_aligned()` now judges column identity against the chance identity
  implied by the data's own residue composition, instead of a fixed 0.5. The
  fixed threshold was unreachable at n = 2, where the commonest residue in a
  column is always at least half of it, and meant different things for
  nucleotides and amino acids.
* `parse_genbank()` given a file path says so, instead of reporting a missing
  `//` terminator.

# viroscape 0.2.0

Methods-review release. Four defects in 0.1.0 silently changed results; all
four are fixed, and the API changes below are deliberate.

## Breaking changes

* `distance_from_consensus(clock_correct = TRUE)` is now an error. Dividing
  divergence by elapsed years is not a rate: on data simulated under a
  perfectly constant rate it produced a strong, overwhelmingly significant
  *negative* association with year. Use `divergence_rate()`.
* `lrt_chain()` returns `nominal_p` and `below_alpha` in place of `p_value` and
  `significant`. The tests are post-selection and anti-conservative, and the
  column names now say so.
* `select_structure(force_time = TRUE)` is the default: `year` is retained
  whether or not it improves the criterion, because divergence accumulates with
  sampling date and year is therefore a confounder rather than a candidate.
* `run_analysis()` and `example_analysis()` model at the country-year level by
  default (`level = "cluster"`).

## Retrieval

* `build_entrez_query()` matches the subtype as an organism rather than with
  `[All Fields]`, which admitted viruses of other subtypes whose records merely
  mentioned H5 somewhere. Country and date are now recall filters only: NCBI has
  no indexed `/country` field, and `[PDAT]` is the deposit date, not the
  collection date. The date range is left open at the top, since a sample cannot
  be deposited before it was collected.
* `fetch_sequences()` downloads GenBank/GenPept flat records and filters exactly
  on the parsed `/country`, `/collection_date` and `/serotype` qualifiers.
* New `parse_genbank()`, `read_genbank()`, `parse_collection_year()` and
  `normalise_country()`.
* `filter_sequences(subtype = )` added.
* `classify_host()` defaults waterfowl to domestic unless explicitly marked
  wild, recognises Latin binomials, and takes a user-supplied `host_map`.
  `vs_host_rules()` exposes the default rules.
* The strain-name pattern allows single internal spaces, so `Hong Kong` and
  `Viet Nam` parse instead of silently becoming `NA`; unparsed deflines warn.
* `resolve_country()` uses the parsed place field only, no longer matching
  country names anywhere in the defline.
* `detect_aligned()` reports its evidence; equal lengths alone no longer pass
  silently as an alignment.

## Inference

* `join_environment(cluster = )`, `cluster_summary()` and
  `aggregate_to_cluster()`. Sequence-level selection on clustered data retained
  at least one pure-noise predictor in 100% of simulated null runs; aggregating
  to the country-year level takes that to 39%, and calibration to 5%.
* `calibrate_selection()` supplies the permutation null a delta has to be read
  against. `print.vs_selection()` says so when it is missing.
* `tidy_coefficients()` defaults to cluster-robust standard errors and reports
  the unit of each coefficient; `unscale()` converts per-SD coefficients back to
  natural units.
* `divergence_rate()` added.
* `report_analysis()` carries a banner when the analysis was built on the
  bundled simulated data, and states the modelling level, the ICC and design
  effect, the coefficient units and the calibration.
* `vs_cache_key()` exported, which was breaking two tests under `R CMD check`.
