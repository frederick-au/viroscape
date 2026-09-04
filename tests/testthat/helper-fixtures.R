## Small deterministic fixtures so the suite runs offline and quickly.

toy_alignment <- function(n = 12, len = 40, seed = 1) {
  set.seed(seed)
  aa <- c("A", "C", "D", "E", "F", "G")
  anc <- sample(aa, len, replace = TRUE)
  rows <- lapply(seq_len(n), function(i) {
    s <- anc
    k <- i %% 5
    if (k > 0) {
      pos <- sample(seq_len(len), k)
      for (p in pos) s[p] <- sample(setdiff(aa, s[p]), 1)
    }
    s
  })
  mat <- do.call(rbind, rows)
  rownames(mat) <- sprintf("seq%02d", seq_len(n))
  mat
}

toy_fasta <- function(mat, path = tempfile(fileext = ".fasta"),
                      countries = c("Cambodia", "Thailand", "Indonesia"),
                      years = 2004:2015) {
  set.seed(2)
  hosts <- c("chicken", "duck", "human")
  headers <- vapply(seq_len(nrow(mat)), function(i) {
    sprintf("EX%05d.1 Influenza A virus (A/%s/%s/AB/%d(H5N1)) hemagglutinin",
            i, sample(hosts, 1), sample(countries, 1), sample(years, 1))
  }, character(1))
  writeLines(as.vector(rbind(paste0(">", headers),
                             apply(mat, 1, paste, collapse = ""))), path)
  path
}

toy_environment <- function(countries = c("Cambodia", "Thailand", "Indonesia"),
                            years = 2004:2015, seed = 3) {
  set.seed(seed)
  grid <- expand.grid(country = countries, year = years,
                      stringsAsFactors = FALSE)
  ## the within-country component must not be a deterministic function of year,
  ## or it is not identifiable once year is in the model
  grid$agricultural_land_pct <- 30 + as.integer(factor(grid$country)) * 5 +
    (grid$year - 2004) * 0.4 + stats::rnorm(nrow(grid), 0, 3)
  grid$livestock_index <- 80 + (grid$year - 2004) * 1.5 +
    stats::rnorm(nrow(grid), 0, 2)
  grid$forest_area_pct <- 60 - as.integer(factor(grid$country)) * 4 +
    stats::rnorm(nrow(grid), 0, 1)
  tibble::as_tibble(grid)
}

small_example <- function(n = 60) {
  s <- viroscape::example_sequences()
  set.seed(7)
  keep <- sort(sample(seq_len(nrow(s$meta)), min(n, nrow(s$meta))))
  s$seqs <- s$seqs[keep]
  s$meta <- s$meta[keep, , drop = FALSE]
  s
}
