## Sequences with real lineage structure, and a country-year predictor that is
## associated with divergence only because lineages sit in particular countries.
## Nothing in the environment causes anything here.
lineage_fixture <- function(n_lin = 4, per = 15, L = 200, seed = 99) {
  set.seed(seed)
  AA <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  anc <- sample(AA, L, TRUE)
  mutate <- function(s, n) {
    if (n > 0) for (i in sample(L, min(n, L))) s[i] <- sample(setdiff(AA, s[i]), 1)
    s
  }
  lin_seq <- lapply(seq_len(n_lin), function(g) mutate(anc, 30))
  lin_extra <- seq(4, 30, length.out = n_lin)
  countries <- LETTERS[1:5]

  rows <- list(); k <- 0
  for (g in seq_len(n_lin)) for (i in seq_len(per)) {
    k <- k + 1
    ## cycle, so a fixture with more lineages than countries still works
    ctry <- if (stats::runif(1) < 0.8) {
      countries[((g - 1) %% length(countries)) + 1]
    } else sample(countries, 1)
    rows[[k]] <- list(
      id = sprintf("SIM%05d.1", k),
      seq = paste(mutate(lin_seq[[g]], stats::rpois(1, lin_extra[g])), collapse = ""),
      country = ctry, year = sample(2003:2022, 1), lineage = sprintf("L%d", g)
    )
  }
  meta <- do.call(rbind, lapply(rows, function(r) {
    data.frame(id = r$id, country = r$country, year = r$year,
               lineage = r$lineage, stringsAsFactors = FALSE)
  }))
  fa <- tempfile(fileext = ".fasta")
  writeLines(unlist(lapply(rows, function(r) c(
    sprintf(">%s Influenza A virus (A/chicken/%s/X/%d(H5N1)) HA", r$id, r$country, r$year),
    r$seq))), fa)

  s <- suppressWarnings(read_sequences(fa))
  s$meta$country <- meta$country
  s$meta$year <- meta$year
  env <- expand.grid(country = countries, year = 2003:2022, stringsAsFactors = FALSE)
  env$x_confounded <- as.numeric(factor(env$country)) + stats::rnorm(nrow(env), 0, 0.3)
  env$x_null <- stats::rnorm(nrow(env))
  list(sequences = s, environment = env, truth = stats::setNames(meta$lineage, meta$id))
}

## A small nucleotide dataset. The package had no nucleotide test at all before
## 0.5.0, which is why the whole ml_nt path could be broken without noticing.
nucleotide_fixture <- function(n = 24, L = 300, seed = 12) {
  set.seed(seed)
  NT <- c("a", "c", "g", "t")
  anc <- sample(NT, L, TRUE)
  countries <- c("Cambodia", "Thailand", "Vietnam")
  rows <- lapply(seq_len(n), function(i) {
    s <- anc
    for (j in sample(L, stats::rpois(1, 12))) s[j] <- sample(setdiff(NT, s[j]), 1)
    list(id = sprintf("NT%05d.1", i), seq = paste(s, collapse = ""),
         country = countries[(i %% 3) + 1], year = 2003 + (i %% 20))
  })
  fa <- tempfile(fileext = ".fasta")
  writeLines(unlist(lapply(rows, function(r) c(
    sprintf(">%s Influenza A virus (A/duck/%s/X/%d(H5N1)) HA", r$id, r$country, r$year),
    r$seq))), fa)
  suppressWarnings(read_sequences(fa, molecule = "nucleotide"))
}
