test_that("nucleotide model names map onto distances that exist", {
  ## dist.ml implements only JC69 and F81 for DNA; optim.pml fits far more, so
  ## a selected model has to be mapped rather than passed through
  expect_equal(nt_distance_model("JC")$model, "JC69")
  expect_true(nt_distance_model("JC")$exact)
  expect_equal(nt_distance_model("F81")$engine, "phangorn")
  expect_equal(nt_distance_model("K80")$engine, "ape")
  expect_equal(nt_distance_model("TN93")$model, "TN93")

  gtr <- nt_distance_model("GTR")
  expect_false(gtr$exact)
  expect_equal(gtr$model, "TN93")
})

test_that("models ape implements exactly are not approximated", {
  ## ape::dist.dna does implement these; mapping them to TN93 and calling it
  ## inexact threw away an exact answer that was available
  for (m in c("K81", "K3P", "T92", "F84", "K80", "TN93")) {
    spec <- nt_distance_model(m)
    expect_true(spec$exact, info = m)
  }
  expect_equal(nt_distance_model("K81")$model, "K81")
  expect_equal(nt_distance_model("T92")$model, "T92")
})

test_that("the substitution preserves the base-frequency assumption", {
  ## SYM has six rates and EQUAL base frequencies. TN93 is richer in rates but
  ## estimates frequencies from the data, which changes the estimator; K81
  ## keeps equal frequencies, so it is the closer neighbour on the axis that
  ## matters for a distance.
  sym <- nt_distance_model("SYM")
  expect_equal(sym$model, "K81")
  expect_equal(sym$family, "equal")
  expect_false(sym$exact)

  gtr <- nt_distance_model("GTR")
  expect_equal(gtr$family, "empirical")
  expect_equal(gtr$model, "TN93")

  ## every mapping stays inside its own frequency family
  for (m in c("JC", "K80", "K81", "SYM", "TVMe")) {
    expect_equal(nt_distance_model(m)$family, "equal", info = m)
  }
  for (m in c("F81", "HKY", "F84", "T92", "TN93", "GTR", "TVM")) {
    expect_equal(nt_distance_model(m)$family, "empirical", info = m)
  }
})

test_that("every mapped model produces a usable distance", {
  aln <- align_sequences(nucleotide_fixture(), method = "none")
  for (m in c("JC", "F81", "K80", "K81", "T92", "F84", "HKY", "TN93", "SYM", "GTR")) {
    d <- suppressMessages(distance_from_consensus(aln, metric = "ml_nt", model = m))
    expect_true(all(is.finite(d$distance)), info = m)
  }
})

test_that("the whole nucleotide pipeline runs", {
  s <- nucleotide_fixture()
  expect_equal(s$molecule, "nucleotide")
  aln <- align_sequences(s, method = "none")

  msel <- select_substitution_model(aln, models = c("JC", "F81", "GTR"),
                                    criterion = "BIC")
  expect_s3_class(msel, "vs_model_selection")
  expect_true(msel$best %in% c("JC", "F81", "GTR"))

  ## this is the call that failed in 0.4.0 with
  ## 'arg' should be one of "JC69", "F81", "WAG", ...
  d <- distance_from_consensus(aln, metric = "ml_nt", model = msel)
  expect_s3_class(d, "vs_distances")
  expect_true(all(is.finite(d$distance)))
  expect_gt(max(d$distance), 0)
})

test_that("every nucleotide metric produces finite distances", {
  aln <- align_sequences(nucleotide_fixture(), method = "none")
  for (m in c("ml_nt", "tn93", "k80", "raw")) {
    d <- distance_from_consensus(aln, metric = m, model = "GTR")
    expect_true(all(is.finite(d$distance)), info = m)
  }
  expect_error(distance_from_consensus(aln, metric = "ml_aa"), "requires protein")
})

test_that("trees and rates work on nucleotides", {
  aln <- align_sequences(nucleotide_fixture(n = 30), method = "none")
  tr <- build_tree(aln, metric = "ml_nt", model = "GTR")
  expect_s3_class(tr, "phylo")
  expect_equal(ape::Ntip(tr), 30L)

  r <- suppressWarnings(divergence_rate(aln, metric = "ml_nt", model = "GTR"))
  expect_s3_class(r, "vs_rate")
})
