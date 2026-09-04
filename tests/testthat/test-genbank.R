test_that("GenBank source qualifiers are parsed, including multi-line values", {
  p <- parse_genbank(readLines(toy_genbank()))

  expect_equal(nrow(p$meta), 2)
  expect_equal(p$meta$accession, c("AB123456.1", "AB123457.1"))
  expect_equal(p$meta$host_raw, c("Cairina moschata", "Gallus gallus"))
  expect_equal(p$meta$subtype, c("H5N1", "H9N2"))
  expect_equal(p$meta$year, c(2005L, 2011L))
  ## the country value wraps onto a second line and carries a sub-national part
  expect_equal(p$meta$country_raw[1], "Viet Nam: Hanoi, Red River Delta")
  expect_equal(unname(p$seqs[1]), "ATGGAGAAAATAGTGCTTCT")
})

test_that("collection dates parse in every layout GenBank allows", {
  expect_equal(
    parse_collection_year(c("2011", "Mar-2011", "15-Mar-2011", "2011-03",
                            "2011-03-15", "2011/2012", "", NA)),
    c(2011L, 2011L, 2011L, 2011L, 2011L, 2011L, NA, NA)
  )
})

test_that("country values are stripped of sub-national detail and normalised", {
  expect_equal(normalise_country(c("Viet Nam: Hanoi", "VIETNAM", "Cambodia",
                                   "Lao PDR", NA)),
               c("Vietnam", "Vietnam", "Cambodia", "Laos", NA))
})

test_that("a GenBank file becomes a vs_sequences with qualifier-derived metadata", {
  s <- read_genbank(toy_genbank(), molecule = "nucleotide")

  expect_s3_class(s, "vs_sequences")
  expect_equal(s$meta$country, c("Vietnam", "Cambodia"))
  expect_equal(s$meta$year, c(2005L, 2011L))
  ## Cairina moschata is the Muscovy duck: domestic, and a binomial
  expect_equal(as.character(s$meta$host_class),
               c("domestic_poultry", "domestic_poultry"))
})

test_that("the subtype guard drops the wrong subtype", {
  s <- read_genbank(toy_genbank(), molecule = "nucleotide")
  f <- filter_sequences(s, subtype = "H5", max_ambiguous = NULL)

  expect_equal(nrow(f$meta), 1)
  expect_equal(f$meta$subtype, "H5N1")
})

test_that("a file path is diagnosed rather than reported as a missing terminator", {
  expect_error(parse_genbank(toy_genbank()), "looks like a file path")
  expect_error(parse_genbank(c("LOCUS  X", "DEFINITION  y")), "record terminator")
})
