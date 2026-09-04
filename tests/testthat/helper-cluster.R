## A dataset with the shape that breaks sequence-level inference: predictors
## that vary only by country-year, many near-identical sequences inside each
## country-year, and no environmental effect whatsoever.
noise_dataset <- function(n_ctry = 5, n_yr = 10, per = 20, n_pred = 7,
                          cluster_sd = 1, resid_sd = 0.27, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cl <- expand.grid(country = LETTERS[seq_len(n_ctry)],
                    year = seq_len(n_yr) + 2002)
  u <- stats::rnorm(nrow(cl), 0, cluster_sd)
  P <- matrix(stats::rnorm(nrow(cl) * n_pred), nrow(cl))
  colnames(P) <- paste0("x", seq_len(n_pred))
  cl <- cbind(cl, P)
  d <- cl[rep(seq_len(nrow(cl)), each = per), ]
  d$distance <- 0.05 + 0.01 * rep(u, each = per) +
    stats::rnorm(nrow(d), 0, resid_sd) * 0.01
  d$country <- factor(d$country)
  rownames(d) <- NULL
  d$.cluster <- interaction(d[c("country", "year")], drop = TRUE, sep = " / ")
  out <- structure(
    list(data = tibble::as_tibble(d), predictors = paste0("x", seq_len(n_pred)),
         response = "distance", keys = c("country", "year"), scaled = FALSE,
         cluster = c("country", "year"), weights = NULL, level = "sequence",
         metric = "ml_aa", model = "JTT"),
    class = "vs_dataset")
  cs <- cluster_summary(out)
  out$icc <- cs$icc
  out$design_effect <- cs$design_effect
  out
}

small_dataset <- function(n = 60) {
  aln <- align_sequences(small_example(n), method = "none")
  d <- distance_from_consensus(aln, metric = "ml_aa", model = "JTT")
  join_environment(d, example_environment())
}

small_analysis <- function(n = 60) {
  s <- small_example(n)
  aln <- align_sequences(s, method = "none")
  msel <- select_substitution_model(aln, models = c("JTT", "WAG"), criterion = "BIC")
  d <- distance_from_consensus(aln, metric = "ml_aa", model = msel)
  ds <- join_environment(d, example_environment())
  mds <- aggregate_to_cluster(ds)
  st <- select_structure(mds, criterion = "AIC")
  sel <- forward_select(mds, base_terms = st, calibrate = FALSE, verbose = FALSE)
  structure(
    list(sequences = s, alignment = aln, model_selection = msel,
         consensus = build_consensus(aln), distances = d,
         environment = example_environment(), dataset = ds,
         model_dataset = mds, level = "cluster", structure = st,
         selection = sel, calibration = NULL, lrt = lrt_chain(sel),
         coefficients = tidy_coefficients(sel), gvif = gvif_table(sel$model),
         fit = model_fit_summary(sel$model), sampling = sampling_summary(ds)),
    class = "vs_analysis")
}

toy_genbank <- function(path = tempfile(fileext = ".gb")) {
  writeLines(c(
    "LOCUS       AB123456                1701 bp    cRNA    linear   VRL 25-JUL-2016",
    "DEFINITION  Influenza A virus (A/duck/Viet Nam/HN-01/2005(H5N1)) HA gene for",
    "            haemagglutinin, complete cds.",
    "ACCESSION   AB123456",
    "VERSION     AB123456.1",
    "FEATURES             Location/Qualifiers",
    "     source          1..1701",
    '                     /organism="Influenza A virus"',
    '                     /strain="A/duck/Viet Nam/HN-01/2005"',
    '                     /serotype="H5N1"',
    '                     /host="Cairina moschata"',
    '                     /country="Viet Nam: Hanoi,',
    '                     Red River Delta"',
    '                     /collection_date="15-Mar-2005"',
    "     gene            1..1701",
    '                     /gene="HA"',
    "ORIGIN",
    "        1 atggagaaaa tagtgcttct",
    "//",
    "LOCUS       AB123457                1701 bp    cRNA    linear   VRL 25-JUL-2016",
    "DEFINITION  Influenza A virus (A/chicken/Cambodia/X1/2011(H5N1)) HA gene.",
    "ACCESSION   AB123457",
    "VERSION     AB123457.1",
    "FEATURES             Location/Qualifiers",
    "     source          1..1701",
    '                     /organism="Influenza A virus"',
    '                     /strain="A/chicken/Cambodia/X1/2011"',
    '                     /serotype="H9N2"',
    '                     /host="Gallus gallus"',
    '                     /country="Cambodia"',
    '                     /collection_date="2011"',
    "ORIGIN",
    "        1 atggagaaaa tagtgcttcg",
    "//"
  ), path)
  path
}
