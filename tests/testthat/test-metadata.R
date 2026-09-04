test_that("influenza deflines are parsed into structured fields", {
  x <- c(
    "AB123456.1 Influenza A virus (A/chicken/Cambodia/X1/2011(H5N1)) hemagglutinin",
    "A/duck/Indonesia/BL/2004(H5N1)",
    "A/Thailand/KAN-1/2004(H5N1)"
  )
  out <- parse_influenza_defline(x)

  expect_equal(nrow(out), 3)
  expect_equal(out$accession[1], "AB123456.1")
  expect_true(is.na(out$accession[2]))
  expect_equal(out$host[1:2], c("chicken", "duck"))
  expect_equal(out$year, c(2011L, 2004L, 2004L))
  expect_equal(out$subtype[1], "H5N1")
})

test_that("strain names containing a space are parsed, not silently dropped", {
  x <- c("A/duck/Hong Kong/205/1977(H5N3)",
         "A/duck/Viet Nam/HN-01/2005(H5N1)",
         "A/chicken/Cambodia/X1/2011(H5N1)")
  out <- parse_influenza_defline(x)

  expect_equal(out$place, c("Hong Kong", "Viet Nam", "Cambodia"))
  expect_equal(out$year, c(1977L, 2005L, 2011L))
  expect_equal(out$country_guess, c("Hong Kong", "Vietnam", "Cambodia"))
})

test_that("a four-digit isolate number is not mistaken for the year", {
  out <- parse_influenza_defline("A/duck/Vietnam/1234/2005(H5N1)")
  expect_equal(out$year, 2005L)
})

test_that("unparsed deflines are reported rather than returned as quiet NA", {
  expect_warning(parse_influenza_defline("not a strain name at all"),
                 "did not match a strain name")
})

test_that("hosts collapse into ecological classes", {
  expect_equal(as.character(classify_host("chicken")), "domestic_poultry")
  expect_equal(as.character(classify_host("mallard")), "wild_bird")
  expect_equal(as.character(classify_host("swine")), "swine")
  expect_equal(as.character(classify_host("human")), "human")
  expect_equal(as.character(classify_host("nonsense-host")), "unknown")
})

test_that("waterfowl default to domestic unless marked wild", {
  expect_equal(as.character(classify_host(c("duck", "domestic duck",
                                            "Muscovy duck", "goose",
                                            "free-grazing duck"))),
               rep("domestic_poultry", 5))
  expect_equal(as.character(classify_host(c("wild duck", "migratory goose"))),
               c("wild_bird", "wild_bird"))
})

test_that("Latin binomials are recognised", {
  expect_equal(as.character(classify_host(c("Gallus gallus", "Cairina moschata",
                                            "Sus scrofa"))),
               c("domestic_poultry", "domestic_poultry", "swine"))
  expect_equal(as.character(classify_host("Anas platyrhynchos")), "wild_bird")
})

test_that("a bare genus is unknown, not domestic", {
  ## "Anas sp." says only "some duck, species not determined" - less information
  ## than a vernacular "duck", so it must not get a more confident answer
  expect_equal(as.character(classify_host(c("Anas sp.", "Anas spp", "anser sp.",
                                            "Anatidae", "waterfowl", "bird"))),
               rep("unknown", 6))
  expect_equal(as.character(classify_host("Anas platyrhynchos")), "wild_bird")
  expect_equal(as.character(classify_host("duck")), "domestic_poultry")
})

test_that("waterfowl is a bird, not an environmental sample", {
  ## the environment rule used to match "water" inside "waterfowl"
  expect_equal(as.character(classify_host(c("wild waterfowl", "waterfowl",
                                            "water sample", "lake water"))),
               c("wild_bird", "unknown", "environment", "environment"))
})

test_that("the host mapping can be replaced by the caller", {
  mine <- list(duck_of_interest = "duck", other = ".")
  expect_equal(as.character(classify_host(c("domestic duck", "chicken"),
                                          host_map = mine)),
               c("duck_of_interest", "other"))
})

test_that("the entrez query filters on subtype as an organism, not free text", {
  q <- build_entrez_query(c("Cambodia", "Thailand"), c(2003, 2010))

  expect_true(grepl('"H5N1 subtype"\\[Organism\\]', q))
  expect_false(grepl("H5N1\\[All Fields\\]", q))
  expect_true(grepl("Influenza A virus", q))
})

test_that("the entrez query does not use a Country field or an upper date bound", {
  q <- build_entrez_query(c("Cambodia", "Thailand"), c(2003, 2010))

  ## [Country] is not an indexed Entrez field; it is silently reinterpreted as
  ## free text, so the country term is written as free text deliberately and
  ## the real filter is applied after download
  expect_false(grepl("\\[Country\\]", q))
  expect_true(grepl("Cambodia", q))
  expect_true(grepl("Thailand", q))

  ## [PDAT] is the deposit date. The lower bound is lossless (nothing can be
  ## deposited before it was collected); an upper bound is not, so there is none
  expect_true(grepl('"2003"\\[PDAT\\] : "3000"\\[PDAT\\]', q))
  expect_false(grepl('"2010"\\[PDAT\\]', q))
})

test_that("alignment evidence is judged against chance, so n = 2 can be rejected", {
  set.seed(5)
  aa <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  ## with two sequences the most frequent residue in a column is always at least
  ## half of it, so any fixed identity threshold at or below 0.5 is unreachable
  pair <- ape::as.AAbin(lapply(1:2, function(i) sample(aa, 80, replace = TRUE)))
  d <- detect_aligned(pair)
  expect_false(isTRUE(d))
  expect_equal(attr(d, "evidence"), "no_evidence")

  same <- sample(aa, 80, replace = TRUE)
  near <- ape::as.AAbin(list(same, replace(same, 1:3, "W")))
  expect_true(isTRUE(detect_aligned(near)))
})

test_that("equal sequence lengths alone are not accepted as evidence of alignment", {
  set.seed(4)
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
          "F", "P", "S", "T", "W", "Y", "V")
  unaligned <- lapply(1:20, function(i) sample(aa, 60, replace = TRUE))
  d1 <- detect_aligned(ape::as.AAbin(unaligned))
  expect_false(isTRUE(d1))
  expect_equal(attr(d1, "evidence"), "no_evidence")

  gapped <- lapply(1:20, function(i) c(rep("-", 3), sample(aa, 57, replace = TRUE)))
  d2 <- detect_aligned(ape::as.AAbin(gapped))
  expect_true(isTRUE(d2))
  expect_equal(attr(d2, "evidence"), "gaps")
})

test_that("sequences read from FASTA carry country and host metadata", {
  path <- toy_fasta(toy_alignment())
  s <- read_sequences(path)

  expect_s3_class(s, "vs_sequences")
  expect_equal(length(s$seqs), 12)
  expect_true(all(s$meta$country %in% c("Cambodia", "Thailand", "Indonesia")))
  expect_false(anyNA(s$meta$year))
})
