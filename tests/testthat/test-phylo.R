test_that("build_tree returns a rooted tree labelled by sequence id", {
  fx <- lineage_fixture(per = 8)
  aln <- align_sequences(fx$sequences, method = "none")
  tr <- build_tree(aln, metric = "p_aa")

  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, aln$meta$id)
  expect_true(ape::is.rooted(tr))
})

test_that("outgroup rooting drops the outgroup from the result", {
  fx <- lineage_fixture(per = 8)
  aln <- align_sequences(fx$sequences, method = "none")
  og <- aln$meta$id[1]
  tr <- build_tree(aln, metric = "p_aa", root = "outgroup", outgroup = og)

  expect_false(any(og %in% tr$tip.label))
  expect_equal(ape::Ntip(tr), nrow(aln$meta) - 1L)
  expect_error(build_tree(aln, metric = "p_aa", root = "outgroup"), "outgroup")
})

test_that("clades cut from the tree recover planted lineages", {
  fx <- lineage_fixture()
  aln <- align_sequences(fx$sequences, method = "none")
  cl <- assign_clades(aln, k = 4, metric = "p_aa")

  expect_length(levels(cl), 4)
  tab <- table(fx$truth, as.character(cl)[match(names(fx$truth), names(cl))])
  ## each true lineage lands in exactly one clade and each clade holds one lineage
  expect_true(all(rowSums(tab > 0) == 1))
  expect_true(all(colSums(tab > 0) == 1))
})

test_that("clade can be used as the clustering variable", {
  fx <- lineage_fixture()
  aln <- align_sequences(fx$sequences, method = "none")
  d <- distance_from_consensus(aln, metric = "p_aa")
  ds <- join_environment(d, fx$environment, scale = FALSE)
  dsc <- add_clades(ds, assign_clades(aln, k = 4, metric = "p_aa"), cluster = TRUE)

  expect_true("clade" %in% names(dsc$data))
  cs <- cluster_summary(dsc)
  expect_equal(cs$n_clusters, 4L)
  expect_gt(cs$icc, 0.5)
})

## The regression test for the last item on the reviewer's clustering list:
## an association that is really lineage structure must not survive the tree.
test_that("a lineage-confounded predictor does not survive phylo_gls", {
  fx <- lineage_fixture()
  aln <- align_sequences(fx$sequences, method = "none")
  d <- distance_from_consensus(aln, metric = "p_aa")
  ds <- join_environment(d, fx$environment, scale = FALSE)
  tr <- build_tree(aln, metric = "p_aa")

  ols <- stats::coef(summary(stats::lm(distance ~ year + x_confounded + x_null,
                                       data = ds$data)))
  expect_lt(ols["x_confounded", 4], 0.001)      # wildly significant without the tree

  pg <- phylo_gls(ds, tr, terms = c("year", "x_confounded", "x_null"))
  row <- pg$comparison[pg$comparison$term == "x_confounded", ]
  expect_gt(row$p_pgls, 0.01)                   # and not with it
  expect_gt(row$shrinkage_pgls, 0.5)            # coefficient collapses
  expect_false(row$survives_pgls)
})

test_that("phylo_gls refuses an aggregated dataset", {
  fx <- lineage_fixture(per = 8)
  aln <- align_sequences(fx$sequences, method = "none")
  ds <- join_environment(distance_from_consensus(aln, metric = "p_aa"),
                         fx$environment, scale = FALSE)
  expect_error(phylo_gls(aggregate_to_cluster(ds), aln), "one row per sequence")
})

test_that("phylo_gls reports a clade fit beside the phylogenetic one", {
  fx <- lineage_fixture()
  aln <- align_sequences(fx$sequences, method = "none")
  ds <- join_environment(distance_from_consensus(aln, metric = "p_aa"),
                         fx$environment, scale = FALSE)
  pg <- phylo_gls(ds, build_tree(aln, metric = "p_aa"),
                  terms = c("year", "x_confounded"), clade_k = 4)

  ## a smooth covariance approximates block structure; a clade factor absorbs
  ## it, so both are reported and the categorical one is the tie-breaker
  expect_true(all(c("p_clade", "shrinkage_clade", "survives_pgls") %in%
                    names(pg$comparison)))
  expect_false("survives" %in% names(pg$comparison))
  row <- pg$comparison[pg$comparison$term == "x_confounded", ]
  expect_gt(row$shrinkage_pgls, 0.5)
  expect_gt(row$p_clade, 0.05)
  expect_output(print(pg), "shrinkage")
})

test_that("the clade fit can be skipped", {
  fx <- lineage_fixture(per = 8)
  aln <- align_sequences(fx$sequences, method = "none")
  ds <- join_environment(distance_from_consensus(aln, metric = "p_aa"),
                         fx$environment, scale = FALSE)
  pg <- phylo_gls(ds, build_tree(aln, metric = "p_aa"),
                  terms = c("year", "x_null"), clades = NA)
  expect_false("p_clade" %in% names(pg$comparison))
})
