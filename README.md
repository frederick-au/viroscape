# viroscape

An R platform for quantifying the relationship between ecological, land-use and
climate variables and viral genetic divergence.

It generalises the pipeline behind *Identifying Strong Abiotic Predictors for
Elevated Mutation Rates in Southeast Asian Viral Isolates of Avian Influenza A*
(Au & Buchan, 2026) into a reusable package, a `targets` pipeline and a Shiny
front-end, and extends it in the directions that paper identified as
limitations.

## What it does

```
retrieve  ->  align  ->  consensus  ->  substitution model  ->  divergence
                                                                    |
environmental predictors  ------------------------------------------+
                                                                    v
                          structural terms -> forward selection -> validation
```

Every stage is a function you can call on its own, and every stage is a
`targets` node so changing one setting reruns only what depends on it.

## Installation

```r
# install.packages("remotes")
remotes::install_local("path/to/viroscape")
```

Core dependencies are all on CRAN. MUSCLE alignment via Bioconductor is
optional:

```r
install.packages("BiocManager")
BiocManager::install(c("msa", "Biostrings"))   # to reproduce the published alignment
```

The package also works with `DECIPHER`, an external `mafft`/`muscle` binary, or
pre-aligned input. `available_aligners()` reports what your machine has.

## Quick start

```r
library(viroscape)

# bundled example: 468 simulated H5-like HA sequences across 12 countries,
# paired with real World Bank indicators
a <- example_analysis()
a

plot_selection_path(a$selection)
plot_coefficients(a)
```

Or run the whole thing live:

```r
a <- run_analysis(
  countries = names(sea_countries()),
  years     = c(2003, 2022),
  molecule  = "protein",
  metric    = "ml_aa"
)
```

Interactively:

```r
launch_app()
```

## The explorer

`launch_app()` opens a five-tab interface over the same functions:

* **Sampling** — how unevenly sequences are distributed across country and year,
  which is the context needed to read country coefficients honestly
* **Divergence** — swap the alignment backend, the substitution model criterion
  and the divergence metric, and see the ranking that results
* **Selection** — choose the country term, the candidate predictors and the
  thresholds, then watch forward selection accept and reject candidates with
  its reasons recorded
* **Diagnostics** — GVIF, the sequential likelihood ratio chain, and the
  coefficient stability table for the final step
* **Results** — coefficients with confidence intervals, and downloads

## What this adds to the original analysis

| Limitation in the paper | What the platform does |
| --- | --- |
| Only five countries, so a country random effect is not estimable | `fetch_environment()` pulls World Bank indicators for any set of countries, and `select_structure()` compares fixed against random and warns below ten groups |
| Amino acid distances only, so synonymous change is invisible | `distance_from_consensus()` supports nucleotide metrics (ML, TN93, K80) alongside the amino acid ones |
| No host metadata, so wild bird and poultry isolates cannot be separated | `parse_influenza_defline()` recovers host from the strain name and `classify_host()` groups it; `filter_sequences(hosts = ...)` stratifies |
| Country-year resolution cannot capture within-country variation | `join_environment()` accepts any join key, including `admin1` for subnational predictors |
| Collinearity handled by manual inspection | `gvif_table()` and `coef_stability()` run at every step, with reasons recorded in the step log |
| Analysis not reproducible without re-running scripts by hand | `_targets.R` orchestrates the whole pipeline with dependency-aware caching |

## Design decisions worth knowing

**BIC selects the substitution model, AIC drives forward selection.** BIC's
penalty scales with alignment length, so it is the more conservative choice
against an over-parameterised substitution matrix. Within a forward-selection
round the candidates are not nested in one another, so a likelihood ratio test
is not valid there and AIC is used instead; once the final model is fixed, the
chain *is* nested and `lrt_chain()` runs the sequential tests.

**GVIF flags, stability decides.** A high GVIF caused by structural correlation
with a grouping factor is not the same problem as redundancy between two
predictors. By default an exceedance is recorded but does not reject a
candidate on its own; what rejects it is a coefficient changing sign and losing
significance, or a significant non-grouping coefficient shifting by more than
the stability threshold *and* by more than one standard error. Pass
`gvif_action = "reject"` for the strict rule.

**No molecular clock by default.** The response is divergence from a regional
consensus and the aim is regression on environmental predictors, not tree
inference, so a strict clock would impose an assumption the analysis does not
need. `clock_correct = TRUE` divides by years elapsed if a per-year rate is
wanted.

## The example data

`example_sequences()` returns 468 **simulated** H5-like haemagglutinin proteins;
`example_environment()` returns **real** World Bank indicators for twelve
countries, 2003-2022. Divergence was generated from those real predictors with a
known ground truth:

| Component | Within-country | Between-country |
| --- | --- | --- |
| Agricultural land (%) | +2.6 | +5.0 |
| Livestock production index | -1.7 | -2.6 |
| Year | +1.15 substitutions/year | |

`example_analysis()` recovers agricultural land as a positive predictor and the
livestock index as a negative one, which is what the test suite asserts. The
sequences are not real isolates: use them to learn the interface, not to draw
biological conclusions.

## Reproducible pipeline

```r
targets::tar_make()          # runs everything, writes output/
targets::tar_visnetwork()    # see the dependency graph
targets::tar_read(selection)
```

Edit the configuration block at the top of `_targets.R` to switch to live data,
change the metric or move the thresholds.

## Caching

Live retrieval is cached on disk so a pipeline run is reproducible offline and
repeated Shiny interactions do not re-hit the network.

```r
vs_cache_dir()     # where it lives
vs_cache_list()    # what is in it
vs_cache_clear()   # start again
```

Set `VIROSCAPE_CACHE` or `options(viroscape.cache = )` to move it.

## Testing

```r
devtools::test()
```

The suite runs offline against small fixtures, plus an end-to-end check that the
planted predictors are recovered from the bundled example.
