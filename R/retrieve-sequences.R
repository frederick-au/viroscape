#' Regular expression matching an influenza strain name
#'
#' Strain names are `A/host/place/isolate/year`, and the place is frequently
#' two words (`Hong Kong`, `Viet Nam`, `New York`). The 0.1.0 pattern excluded
#' whitespace outright, so every such strain parsed to `NA` and was then dropped
#' silently by [join_environment()]. Single internal spaces are allowed here;
#' the year is anchored to `19xx`/`20xx` so a four-digit isolate number cannot
#' be mistaken for one.
#'
#' @return A single Perl-compatible regular expression.
#' @export
vs_strain_pattern <- function() {
  tok <- "[A-Za-z0-9._'+-]+(?: [A-Za-z0-9._'+-]+)*"
  sprintf("A/%s(?:/%s)*/(?:19|20)[0-9]{2}", tok, tok)
}

#' Parse influenza strain deflines into structured metadata
#'
#' Influenza sequence names follow the convention `A/host/place/isolate/year`
#' (the host field is omitted for human isolates). Recovering the host is what
#' allows analyses to be stratified by wild bird, domestic poultry and mammalian
#' isolates, which the country-year design of the original study could not do.
#'
#' A defline is a weak source of truth: it is free text, and the place field is
#' not always a country. Prefer [read_genbank()] or [fetch_sequences()], which
#' read the `source` feature's qualifiers instead, and treat this as the
#' fallback it is.
#'
#' @param x Character vector of FASTA headers or strain names.
#' @param host_map Passed to [classify_host()].
#' @param warn Warn about headers that could not be parsed. Silent parsing is
#'   how unparsed records used to disappear without trace.
#' @return A tibble with columns `accession`, `strain`, `host`, `place`,
#'   `year`, `subtype`, `host_class` and `country_guess`.
#' @examples
#' parse_influenza_defline(c(
#'   "AB123456.1 Influenza A virus (A/chicken/Cambodia/X1/2011(H5N1)) HA",
#'   "A/duck/Hong Kong/205/1977(H5N3)"
#' ))
#' @export
parse_influenza_defline <- function(x, host_map = NULL, warn = TRUE) {
  x <- as.character(x)
  acc <- sub("^>?\\s*([^ |]+).*$", "\\1", x)
  acc[!grepl("^[A-Z]{1,3}[0-9]{5,8}(\\.[0-9]+)?$", acc)] <- NA_character_

  strain <- rep(NA_character_, length(x))
  pat <- vs_strain_pattern()
  hit <- grepl(pat, x, perl = TRUE)
  strain[hit] <- regmatches(x[hit], regexpr(pat, x[hit], perl = TRUE))

  subtype <- rep(NA_character_, length(x))
  st_pat <- "H[0-9]{1,2}N[0-9]{1,2}|H[0-9]{1,2}Nx"
  has_st <- grepl(st_pat, x, perl = TRUE)
  subtype[has_st] <- regmatches(x[has_st], regexpr(st_pat, x[has_st], perl = TRUE))

  parts <- strsplit(sub("^A/", "", ifelse(is.na(strain), "", strain)), "/", fixed = TRUE)
  host <- vapply(parts, function(p) if (length(p) >= 4) p[1] else NA_character_, character(1))
  place <- vapply(parts, function(p) {
    if (length(p) >= 4) p[2] else if (length(p) >= 1) p[1] else NA_character_
  }, character(1))
  yr <- vapply(parts, function(p) {
    if (!length(p)) return(NA_character_)
    tl <- p[length(p)]
    if (grepl("^[0-9]{4}$", tl)) tl else NA_character_
  }, character(1))
  host[is.na(host) & !is.na(strain)] <- "human"
  place[!nzchar(place %||% "")] <- NA_character_

  if (warn && any(!hit)) {
    vs_warn(c(
      "{sum(!hit)} of {length(x)} defline{?s} did not match a strain name and will carry no host, place or year.",
      "i" = "First unparsed: {.val {utils::head(x[!hit], 1)}}"
    ))
  }

  tibble::tibble(
    accession = acc,
    strain = strain,
    host = tolower(host),
    place = place,
    year = suppressWarnings(as.integer(yr)),
    subtype = subtype,
    host_class = classify_host(host, host_map = host_map),
    country_guess = normalise_country(place)
  )
}

#' Build an NCBI Entrez query for influenza sequences
#'
#' Exposed separately so a query can be inspected (and pasted into the NCBI web
#' interface) before it is run.
#'
#' Two of the terms are deliberately imprecise, and this is the point. NCBI has
#' no indexed field for the `/country` source qualifier, and `[PDAT]` is the
#' date a record was *released*, not the date the sample was collected - in a
#' study whose central weakness is temporal confounding, filtering on deposit
#' date silently deletes records with a country-specific bias. Both terms are
#' therefore used only as recall filters to keep the download tractable, and are
#' applied exactly, after download, to the parsed `/country` and
#' `/collection_date` qualifiers by [fetch_sequences()]:
#'
#' * country enters as free text, which is over-inclusive rather than wrong;
#' * the date range is left open at the top (`years[1]` to `3000`), because a
#'   sample collected in `years[1]` or later cannot have been deposited before
#'   `years[1]`. The lower bound is therefore lossless; an upper bound is not.
#'
#' The subtype is matched as an organism (`"H5N1 subtype"[Organism]`) rather
#' than with `[All Fields]`, which matches any record whose text mentions the
#' subtype anywhere - including reference titles - and so admits viruses of
#' entirely different subtypes.
#'
#' @param countries Character vector of country names, or `NULL` for no country
#'   term at all.
#' @param years Length-2 integer vector giving the inclusive *collection* year
#'   range. Only the lower bound enters the query.
#' @param subtype Haemagglutinin subtype, e.g. `"H5"` for all H5Nx, or a full
#'   subtype such as `"H5N1"`.
#' @param segment Gene segment name, e.g. `"hemagglutinin"`.
#' @param molecule `"protein"` or `"nucleotide"`.
#' @param subtype_field `"organism"` uses the subtype taxonomy node;
#'   `"all_fields"` restores the 0.1.0 behaviour, which is imprecise and kept
#'   only as an escape hatch.
#' @return A single Entrez query string.
#' @export
build_entrez_query <- function(countries = NULL, years = c(2003, 2022),
                               subtype = "H5", segment = "hemagglutinin",
                               molecule = c("protein", "nucleotide"),
                               subtype_field = c("organism", "all_fields")) {
  molecule <- match.arg(molecule)
  subtype_field <- match.arg(subtype_field)

  sub_terms <- if (grepl("^H[0-9]{1,2}$", subtype)) {
    sprintf("%sN%d", subtype, 1:9)
  } else if (grepl("^H[0-9]{1,2}N[0-9]{1,2}$", subtype)) {
    subtype
  } else {
    vs_warn("{.arg subtype} {.val {subtype}} is not an H or HN code; falling back to a free-text term.")
    NA_character_
  }
  sub_field <- if (anyNA(sub_terms)) {
    sprintf('%s[All Fields]', subtype)
  } else if (subtype_field == "organism") {
    paste(sprintf('"%s subtype"[Organism]', sub_terms), collapse = " OR ")
  } else {
    paste(sprintf('%s[All Fields]', sub_terms), collapse = " OR ")
  }

  seg_field <- if (molecule == "protein") {
    sprintf('(%s[Protein Name] OR %s[Title])', segment, segment)
  } else {
    sprintf('(%s[Gene Name] OR %s[Title])', segment, segment)
  }

  q <- paste0('("Influenza A virus"[Organism]) AND (', sub_field, ') AND ', seg_field)
  if (length(countries)) {
    q <- paste0(q, ' AND (',
                paste(sprintf('"%s"[All Fields]', countries), collapse = " OR "), ')')
  }
  if (length(years)) {
    q <- paste0(q, sprintf(' AND ("%d"[PDAT] : "3000"[PDAT])', as.integer(years[1])))
  }
  q
}

#' Retrieve influenza sequences from NCBI
#'
#' Live retrieval with on-disk caching. Requires the `rentrez` package.
#'
#' GenBank/GenPept flat records are downloaded rather than FASTA, because
#' country, collection date and host are `source` feature qualifiers and are not
#' present in a defline. The Entrez query is a recall filter only (see
#' [build_entrez_query()]); the requested countries, collection years and
#' subtype are then applied exactly to the parsed qualifiers, and the number of
#' records each filter removes is reported.
#'
#' @inheritParams build_entrez_query
#' @param max_records Maximum number of records to download.
#' @param cache_dir Cache directory; see [vs_cache_dir()].
#' @param refresh Force a fresh download.
#' @param api_key Optional NCBI API key (raises the rate limit).
#' @param strict Apply the exact post-download filters. `FALSE` returns
#'   everything the query matched, which is useful for auditing what the recall
#'   filter actually pulled in.
#' @param host_map Passed to [classify_host()].
#' @return A `vs_sequences` object: a list with `seqs` (an `AAbin` or `DNAbin`
#'   list), `meta` (a tibble of per-sequence metadata) and `query`.
#' @export
fetch_sequences <- function(countries, years = c(2003, 2022), subtype = "H5",
                            segment = "hemagglutinin",
                            molecule = c("protein", "nucleotide"),
                            max_records = 5000,
                            cache_dir = vs_cache_dir(), refresh = FALSE,
                            api_key = Sys.getenv("ENTREZ_KEY", ""),
                            strict = TRUE, host_map = NULL) {
  molecule <- match.arg(molecule)
  need_pkg("rentrez", "download sequences from NCBI")
  query <- build_entrez_query(countries, years, subtype, segment, molecule)
  db <- if (molecule == "protein") "protein" else "nuccore"
  rettype <- if (molecule == "protein") "gp" else "gb"
  key <- vs_cache_key("sequences_gb", query, db, max_records)

  raw <- with_cache(key, cache_dir = cache_dir, refresh = refresh, expr = {
    if (nzchar(api_key)) rentrez::set_entrez_key(api_key)
    vs_inform("Querying NCBI {.field {db}} ...")
    srch <- rentrez::entrez_search(db = db, term = query,
                                   retmax = max_records, use_history = TRUE)
    if (srch$count == 0) vs_abort("NCBI returned no records for this query.")
    if (srch$count > max_records) {
      vs_warn(c(
        "The query matched {srch$count} record{?s} but {.arg max_records} is {max_records}.",
        "i" = "The download is truncated, and truncation is not random. Raise {.arg max_records} or narrow the query."
      ))
    }
    vs_inform("Found {srch$count} record{?s}; downloading {min(srch$count, max_records)}.")
    n <- min(srch$count, max_records)
    chunks <- seq(0, n - 1, by = 200)
    txt <- vapply(chunks, function(start) {
      rentrez::entrez_fetch(db = db, web_history = srch$web_history,
                            rettype = rettype, retmode = "text",
                            retstart = start, retmax = 200)
    }, character(1))
    paste(txt, collapse = "")
  })

  parsed <- parse_genbank(raw)
  out <- vs_sequences_from_genbank(parsed, molecule = molecule, host_map = host_map,
                                   source = "ncbi", query = query)
  if (strict) {
    out <- apply_retrieval_filters(out, countries = countries, years = years,
                                   subtype = subtype)
  }
  out
}

## Exact filters on parsed qualifiers, with a per-filter account of what went.
apply_retrieval_filters <- function(x, countries = NULL, years = NULL, subtype = NULL) {
  n0 <- nrow(x$meta)
  drops <- c()
  keep <- rep(TRUE, n0)

  if (length(countries)) {
    want <- normalise_country(countries)
    hit <- !is.na(x$meta$country) & x$meta$country %in% want
    drops <- c(drops, country = sum(keep & !hit))
    keep <- keep & hit
  }
  if (length(years) == 2L) {
    hit <- !is.na(x$meta$year) & x$meta$year >= years[1] & x$meta$year <= years[2]
    drops <- c(drops, collection_year = sum(keep & !hit))
    keep <- keep & hit
  }
  if (!is.null(subtype)) {
    hit <- subtype_matches(x$meta$subtype, subtype)
    drops <- c(drops, subtype = sum(keep & !hit))
    keep <- keep & hit
  }
  if (!any(keep)) {
    vs_abort("Every downloaded record was removed by the exact filters; check the query.")
  }
  if (sum(!keep)) {
    detail <- paste(sprintf("%s: %d", names(drops), drops), collapse = ", ")
    vs_inform("Exact filters removed {sum(!keep)} of {n0} downloaded record{?s} ({detail}).")
  }
  x$seqs <- x$seqs[keep]
  x$meta <- x$meta[keep, , drop = FALSE]
  x
}

subtype_matches <- function(x, subtype) {
  if (is.null(subtype)) return(rep(TRUE, length(x)))
  pat <- if (grepl("^H[0-9]{1,2}$", subtype)) {
    sprintf("^%s(N[0-9]{1,2}|Nx)?$", subtype)
  } else {
    sprintf("^%s$", subtype)
  }
  !is.na(x) & grepl(pat, x)
}

#' Read sequences from a local FASTA file
#'
#' @param path Path to a FASTA file (aligned or unaligned).
#' @param metadata Optional data frame of per-sequence metadata. Must contain a
#'   column matching the FASTA names (`accession` or `strain`); anything it does
#'   not supply is parsed from the deflines.
#' @param molecule `"protein"` or `"nucleotide"`.
#' @param country Optional country label applied to all sequences, or a vector
#'   the same length as the file.
#' @param host_map Passed to [classify_host()].
#' @return A `vs_sequences` object.
#' @export
read_sequences <- function(path, metadata = NULL,
                           molecule = c("protein", "nucleotide"),
                           country = NULL, host_map = NULL) {
  molecule <- match.arg(molecule)
  if (!file.exists(path)) vs_abort("File {.path {path}} does not exist.")
  lines <- readLines(path, warn = FALSE)
  headers <- sub("^>", "", grep("^>", lines, value = TRUE))
  if (!length(headers)) vs_abort("{.path {path}} contains no FASTA records.")

  seqs <- if (molecule == "nucleotide") {
    ape::read.FASTA(path, type = "DNA")
  } else {
    ape::read.FASTA(path, type = "AA")
  }
  names(seqs) <- headers

  meta <- parse_influenza_defline(headers, host_map = host_map)
  meta$label <- headers
  meta$id <- make.unique(ifelse(is.na(meta$accession),
                                substr(headers, 1, 40), meta$accession))
  names(seqs) <- meta$id

  meta$country <- resolve_country(meta$country_guess, country)
  if (!is.null(metadata)) meta <- merge_metadata(meta, metadata)

  al <- detect_aligned(seqs)
  structure(
    list(seqs = seqs, meta = meta, molecule = molecule,
         aligned = isTRUE(al), aligned_evidence = attr(al, "evidence"),
         query = NA_character_, source = path),
    class = "vs_sequences"
  )
}

## Country comes from the parsed place field only. 0.1.0 also grepped every
## country name against the whole defline, so a record merely *mentioning*
## another country could be assigned to it - the same free-text failure as the
## old [Country] Entrez term, but client side.
resolve_country <- function(country_guess, country = NULL) {
  if (!is.null(country) && length(country) == 1L) {
    return(rep(country, length(country_guess)))
  }
  if (!is.null(country) && length(country) == length(country_guess)) return(country)
  country_guess
}

#' Is a set of sequences already aligned?
#'
#' Equal lengths alone are weak evidence: full-length influenza HA records are
#' very often all exactly the same length while being entirely unaligned. Gap
#' characters are strong evidence. Where lengths agree but no gaps are present,
#' the decision is made on how much more similar the columns are than chance.
#'
#' "Than chance" rather than "above a fixed proportion" matters. A fixed
#' threshold on mean column identity is unreachable for small `n` - with two
#' sequences the most frequent residue in a column always accounts for at least
#' half of it - and it means something different for nucleotides, where chance
#' agreement is around 1 in 4, than for amino acids, where it is nearer 1 in 16.
#' Both are handled by comparing observed mean pairwise identity against the
#' identity implied by the residue composition of the data itself.
#'
#' @param seqs An `AAbin`/`DNAbin` list or matrix.
#' @param min_identity Required *excess over chance* identity, on a scale where
#'   0 is chance agreement and 1 is a perfectly conserved column.
#' @return `TRUE`/`FALSE`, with an `evidence` attribute of `"matrix"`,
#'   `"gaps"`, `"conserved_columns"`, `"single_sequence"`, `"unequal_lengths"`
#'   or `"no_evidence"`, and `identity`/`chance` attributes carrying the
#'   numbers behind the decision.
#' @examples
#' set.seed(1)
#' aa <- c("A", "C", "D", "E", "F", "G")
#' unaligned <- lapply(1:2, function(i) sample(aa, 60, replace = TRUE))
#' detect_aligned(ape::as.AAbin(unaligned))   # FALSE even at n = 2
#' @export
detect_aligned <- function(seqs, min_identity = 0.5) {
  ev <- function(v, e, id = NA_real_, ch = NA_real_) {
    structure(v, evidence = e, identity = id, chance = ch)
  }
  if (is.matrix(seqs)) return(ev(TRUE, "matrix"))
  lens <- lengths(seqs)
  if (length(unique(lens)) != 1L) return(ev(FALSE, "unequal_lengths"))
  ch <- as.character(seqs)
  m <- if (is.list(ch)) do.call(rbind, ch) else matrix(ch, nrow = length(seqs), byrow = TRUE)
  if (any(m %in% c("-", ".", "~"))) return(ev(TRUE, "gaps"))
  n <- nrow(m)
  if (n < 2L) return(ev(TRUE, "single_sequence"))

  ## mean probability that two sequences drawn at random agree at a column
  pid <- mean(apply(m, 2, function(col) {
    tb <- table(col)
    sum(tb * (tb - 1)) / (n * (n - 1))
  }))
  ## the same probability implied by the overall residue composition
  tb <- table(m)
  N <- sum(tb)
  chance <- if (N > 1) sum(tb * (tb - 1)) / (N * (N - 1)) else NA_real_
  excess <- (pid - chance) / (1 - chance)

  if (isTRUE(excess >= min_identity)) {
    ev(TRUE, "conserved_columns", pid, chance)
  } else {
    ev(FALSE, "no_evidence", pid, chance)
  }
}

merge_metadata <- function(meta, metadata) {
  metadata <- tibble::as_tibble(metadata)
  keys <- intersect(c("id", "accession", "strain", "label"), names(metadata))
  if (!length(keys)) {
    vs_warn("{.arg metadata} has no join column; ignoring it.")
    return(meta)
  }
  key <- keys[1]
  extra <- setdiff(names(metadata), names(meta))
  overlap <- intersect(setdiff(names(metadata), key), names(meta))
  joined <- dplyr::left_join(meta, metadata[c(key, extra, overlap)],
                             by = key, suffix = c("", ".supplied"))
  for (col in overlap) {
    sup <- paste0(col, ".supplied")
    if (sup %in% names(joined)) {
      joined[[col]] <- ifelse(is.na(joined[[sup]]), joined[[col]], joined[[sup]])
      joined[[sup]] <- NULL
    }
  }
  joined
}

#' @export
print.vs_sequences <- function(x, ...) {
  cli::cli_h3("vs_sequences")
  cli::cli_text("{length(x$seqs)} {x$molecule} sequence{?s}, aligned: {x$aligned}{if (!is.null(x$aligned_evidence)) paste0(' (', x$aligned_evidence, ')') else ''}")
  ctry <- sort(table(x$meta$country), decreasing = TRUE)
  if (length(ctry)) {
    cli::cli_text("Countries: {paste(sprintf('%s (%d)', names(ctry), ctry), collapse = ', ')}")
  }
  yr <- suppressWarnings(range(x$meta$year, na.rm = TRUE))
  if (all(is.finite(yr))) cli::cli_text("Years: {yr[1]}-{yr[2]}")
  hosts <- sort(table(x$meta$host_class), decreasing = TRUE)
  hosts <- hosts[hosts > 0]
  if (length(hosts)) {
    cli::cli_text("Hosts: {paste(sprintf('%s (%d)', names(hosts), hosts), collapse = ', ')}")
  }
  invisible(x)
}

#' Filter a vs_sequences object
#'
#' @param x A `vs_sequences` object.
#' @param countries Optional character vector to keep.
#' @param years Optional length-2 year range.
#' @param hosts Optional host classes to keep (see [classify_host()]).
#' @param subtype Optional subtype guard, e.g. `"H5"` keeps H5Nx only. Records
#'   whose subtype could not be determined are dropped, because an unverifiable
#'   subtype is exactly what a free-text query lets through.
#' @param min_length Drop sequences shorter than this many residues.
#' @param max_ambiguous Maximum proportion of ambiguous characters allowed.
#' @return A filtered `vs_sequences` object.
#' @export
filter_sequences <- function(x, countries = NULL, years = NULL, hosts = NULL,
                             subtype = NULL, min_length = NULL,
                             max_ambiguous = 0.05) {
  stopifnot(inherits(x, "vs_sequences"))
  keep <- rep(TRUE, nrow(x$meta))
  if (!is.null(countries)) keep <- keep & x$meta$country %in% countries
  if (!is.null(years)) keep <- keep & !is.na(x$meta$year) &
      x$meta$year >= years[1] & x$meta$year <= years[2]
  if (!is.null(hosts)) keep <- keep & as.character(x$meta$host_class) %in% hosts
  if (!is.null(subtype)) {
    hit <- subtype_matches(x$meta$subtype, subtype)
    if (any(keep & !hit)) {
      vs_inform("Subtype filter {.val {subtype}} removes {sum(keep & !hit)} sequence{?s} ({sum(keep & is.na(x$meta$subtype))} with no recorded subtype).")
    }
    keep <- keep & hit
  }
  lens <- lengths(x$seqs)
  if (!is.null(min_length)) keep <- keep & lens >= min_length
  if (!is.null(max_ambiguous)) {
    amb <- ambiguous_fraction(x$seqs, x$molecule)
    keep <- keep & amb <= max_ambiguous
  }
  if (!any(keep)) vs_abort("No sequences remain after filtering.")
  x$seqs <- x$seqs[keep]
  x$meta <- x$meta[keep, , drop = FALSE]
  x
}

ambiguous_fraction <- function(seqs, molecule) {
  chars <- as.character(seqs)
  vapply(chars, function(s) {
    s <- toupper(s)
    bad <- if (molecule == "nucleotide") {
      sum(!s %in% c("A", "C", "G", "T", "-"))
    } else {
      sum(s %in% c("X", "?", "*"))
    }
    bad / max(length(s), 1)
  }, numeric(1), USE.NAMES = FALSE)
}
