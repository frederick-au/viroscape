test_that("model selection ranks every candidate and picks a winner", {
  aln <- align_sequences(read_sequences(toy_fasta(toy_alignment(n = 8, len = 80))),
                         method = "none")
  sel <- select_substitution_model(aln, models = c("JTT", "WAG", "LG"),
                                   criterion = "BIC")

  expect_s3_class(sel, "vs_model_selection")
  expect_equal(nrow(sel$table), 3)
  expect_true(sel$best %in% c("JTT", "WAG", "LG"))
  expect_equal(sel$table$model[1], sel$best)
  expect_true(!is.unsorted(sel$table$BIC))
  expect_equal(min(sel$table$delta), 0)
})

test_that("BIC penalises parameters more heavily than AIC", {
  aln <- align_sequences(read_sequences(toy_fasta(toy_alignment(n = 8, len = 80))),
                         method = "none")
  sel <- select_substitution_model(aln, models = c("JTT", "WAG"), criterion = "BIC")
  tab <- sel$table
  gap <- tab$BIC - tab$AIC
  expect_true(all(gap > 0))
})

test_that("candidate model lists differ by molecule", {
  expect_true("JTT" %in% candidate_models("protein"))
  expect_true("GTR" %in% candidate_models("nucleotide"))
  expect_false("JTT" %in% candidate_models("nucleotide"))
})
