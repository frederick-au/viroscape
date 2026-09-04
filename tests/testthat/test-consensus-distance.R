test_that("majority rule picks the most common residue at each site", {
  mat <- rbind(c("A", "C", "D"), c("A", "C", "E"), c("A", "G", "D"))
  rownames(mat) <- c("a", "b", "c")
  aln <- structure(list(matrix = mat,
                        phy = phangorn::phyDat(mat, type = "AA"),
                        meta = tibble::tibble(id = rownames(mat)),
                        molecule = "protein", method = "none", n_sites = 3),
                   class = "vs_alignment")
  cons <- build_consensus(aln)

  expect_s3_class(cons, "vs_consensus")
  expect_equal(cons$sequence, c("A", "C", "D"))
  expect_equal(cons$support[1], 1)
})

test_that("the threshold rule emits an ambiguity character below support", {
  mat <- rbind(c("A", "C"), c("D", "C"), c("E", "C"))
  rownames(mat) <- c("a", "b", "c")
  aln <- structure(list(matrix = mat,
                        phy = phangorn::phyDat(mat, type = "AA"),
                        meta = tibble::tibble(id = rownames(mat)),
                        molecule = "protein", method = "none", n_sites = 2),
                   class = "vs_alignment")
  cons <- build_consensus(aln, method = "threshold", threshold = 0.5)
  expect_equal(cons$sequence[1], "X")
  expect_equal(cons$sequence[2], "C")
})

test_that("divergence increases with the number of substitutions", {
  mat <- toy_alignment(n = 10, len = 60)
  path <- toy_fasta(mat)
  seqs <- read_sequences(path)
  aln <- align_sequences(seqs, method = "none")
  d <- distance_from_consensus(aln, metric = "p_aa")

  expect_s3_class(d, "vs_distances")
  expect_equal(nrow(d), 10)
  expect_true(all(d$distance >= 0))

  ## FASTA records keep their input order, so substitution counts line up by
  ## position rather than by name
  n_sub <- apply(mat, 1, function(r) sum(r != mat[1, ]))
  ord <- match(d$id, names(seqs$seqs))
  expect_gt(stats::cor(d$distance, n_sub[ord]), 0.5)
})

test_that("a nucleotide metric is refused for a protein alignment", {
  aln <- align_sequences(read_sequences(toy_fasta(toy_alignment())), method = "none")
  expect_error(distance_from_consensus(aln, metric = "tn93"), "nucleotide")
})
