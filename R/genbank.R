#' Parse GenBank / GenPept flat files
#'
#' The package retrieves flat records rather than FASTA because the fields that
#' matter for this analysis - collection date, country of origin and host - live
#' in the `source` feature's qualifiers and are simply absent from a FASTA
#' defline. Parsing them is what allows the country and date filters to be
#' applied to what the record actually says, instead of to a free-text match
#' over the whole record (see [build_entrez_query()]).
#'
#' @param x A character vector of lines, or a single string containing an
#'   entire flat file.
#' @return A list with `meta` (one row per record) and `seqs` (a named
#'   character vector of sequence strings, upper case, gaps absent).
#' @export
parse_genbank <- function(x) {
  if (length(x) == 1L && !grepl("\n", x) && nzchar(x) && file.exists(x)) {
    vs_abort(c(
      "{.arg x} looks like a file path, not the contents of a flat file.",
      "i" = "Use {.code read_genbank({.val {x}})}, or pass {.code readLines({.val {x}})}."
    ))
  }
  lines <- if (length(x) == 1L && grepl("\n", x)) strsplit(x, "\n", fixed = TRUE)[[1]] else as.character(x)
  lines <- sub("\r$", "", lines)
  ends <- which(trimws(lines) == "//")
  if (!length(ends)) {
    vs_abort(c(
      "No GenBank records found: no {.code //} record terminator in {length(lines)} line{?s}.",
      "i" = "Records must be GenBank or GenPept flat files. FASTA goes to {.fn read_sequences}."
    ))
  }
  starts <- c(1L, utils::head(ends, -1L) + 1L)

  recs <- lapply(seq_along(ends), function(i) parse_genbank_record(lines[starts[i]:ends[i]]))
  recs <- Filter(Negate(is.null), recs)
  if (!length(recs)) vs_abort("No parseable GenBank records.")

  meta <- dplyr::bind_rows(lapply(recs, function(r) r$meta))
  seqs <- vapply(recs, function(r) r$sequence, character(1))
  names(seqs) <- meta$accession
  list(meta = meta, seqs = seqs)
}

parse_genbank_record <- function(rl) {
  if (!length(rl)) return(NULL)
  version <- grep("^VERSION", rl, value = TRUE)
  locus <- grep("^LOCUS", rl, value = TRUE)
  acc <- if (length(version)) {
    strsplit(trimws(sub("^VERSION", "", version[1])), "\\s+")[[1]][1]
  } else if (length(locus)) {
    strsplit(trimws(sub("^LOCUS", "", locus[1])), "\\s+")[[1]][1]
  } else NA_character_
  if (is.na(acc)) return(NULL)

  defn <- collapse_field(rl, "DEFINITION")
  quals <- source_qualifiers(rl)

  org <- quals[["organism"]] %||% NA_character_
  strain <- quals[["strain"]] %||% quals[["isolate"]] %||% extract_strain(c(defn, org))
  serotype <- quals[["serotype"]] %||% NA_character_
  if (is.na(serotype)) {
    st <- regmatches(paste(defn, org), regexpr("H[0-9]{1,2}N[0-9]{1,2}", paste(defn, org)))
    if (length(st)) serotype <- st
  }
  cdate <- quals[["collection_date"]] %||% NA_character_

  ori <- grep("^ORIGIN", rl)
  sq <- if (length(ori)) {
    body <- rl[(ori[1] + 1L):length(rl)]
    body <- body[trimws(body) != "//"]
    toupper(gsub("[^A-Za-z*-]", "", paste(body, collapse = "")))
  } else ""

  list(
    sequence = sq,
    meta = tibble::tibble(
      accession = acc,
      definition = defn %||% NA_character_,
      strain = strain %||% NA_character_,
      subtype = serotype,
      host_raw = quals[["host"]] %||% NA_character_,
      country_raw = quals[["country"]] %||% quals[["geo_loc_name"]] %||% NA_character_,
      collection_date = cdate,
      year = parse_collection_year(cdate),
      segment = quals[["segment"]] %||% NA_character_
    )
  )
}

collapse_field <- function(rl, key) {
  i <- grep(paste0("^", key), rl)
  if (!length(i)) return(NA_character_)
  i <- i[1]
  j <- i
  while (j + 1L <= length(rl) && grepl("^\\s", rl[j + 1L]) && !grepl("^\\s{5}\\S", rl[j + 1L])) {
    j <- j + 1L
  }
  trimws(gsub("\\s+", " ", paste(sub(paste0("^", key), "", rl[i:j]), collapse = " ")))
}

## Qualifiers of the first (source) feature. Continuation lines are folded in;
## a value is only closed on its terminating quote, so multi-line countries and
## hosts survive.
source_qualifiers <- function(rl) {
  fi <- grep("^FEATURES", rl)
  if (!length(fi)) return(list())
  start <- fi[1] + 1L
  ori <- grep("^ORIGIN", rl)
  stop_at <- if (length(ori)) ori[1] - 1L else length(rl)
  if (start > stop_at) return(list())
  block <- rl[start:stop_at]
  ## the source feature runs until the next feature key at the 5-space indent
  keys <- grep("^\\s{5}\\S", block)
  if (length(keys) >= 2) block <- block[seq_len(keys[2] - 1L)]

  out <- list()
  cur <- NULL
  buf <- character()
  flush <- function() {
    if (!is.null(cur) && is.null(out[[cur]])) {
      v <- paste(buf, collapse = " ")
      v <- gsub('^"|"$', "", trimws(v))
      out[[cur]] <<- trimws(gsub("\\s+", " ", v))
    }
  }
  for (ln in block) {
    t <- trimws(ln)
    if (startsWith(t, "/")) {
      flush()
      if (grepl("=", t, fixed = TRUE)) {
        cur <- sub("^/([^=]+)=.*$", "\\1", t)
        buf <- sub("^/[^=]+=", "", t)
      } else {
        cur <- sub("^/", "", t)
        buf <- "TRUE"
      }
    } else if (!is.null(cur) && nzchar(t)) {
      buf <- c(buf, t)
    }
  }
  flush()
  out
}

extract_strain <- function(x) {
  x <- paste(stats::na.omit(x), collapse = " ")
  m <- regmatches(x, regexpr(vs_strain_pattern(), x, perl = TRUE))
  if (length(m)) m else NA_character_
}

#' Collection year from a GenBank `/collection_date` qualifier
#'
#' Accepts every layout GenBank permits (`2011`, `Mar-2011`, `15-Mar-2011`,
#' `2011-03`, `2011-03-15`) and ranges such as `2011/2012`, for which the
#' earlier year is returned.
#'
#' @param x Character vector of `/collection_date` values.
#' @return Integer vector of years, `NA` where no year is present.
#' @export
parse_collection_year <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexpr("(19|20)[0-9]{2}", x))
  out <- rep(NA_integer_, length(x))
  has <- lengths(regmatches(x, gregexpr("(19|20)[0-9]{2}", x))) > 0
  out[has] <- as.integer(m)
  out
}

#' Read sequences from a GenBank or GenPept flat file
#'
#' @param path Path to a `.gb` / `.gp` flat file.
#' @param molecule `"protein"` or `"nucleotide"`.
#' @param host_map Passed to [classify_host()].
#' @return A `vs_sequences` object whose metadata comes from the record's
#'   `source` qualifiers rather than from the defline.
#' @export
read_genbank <- function(path, molecule = c("nucleotide", "protein"),
                         host_map = NULL) {
  molecule <- match.arg(molecule)
  if (!file.exists(path)) vs_abort("File {.path {path}} does not exist.")
  parsed <- parse_genbank(readLines(path, warn = FALSE))
  vs_sequences_from_genbank(parsed, molecule = molecule, host_map = host_map,
                            source = path)
}

vs_sequences_from_genbank <- function(parsed, molecule, host_map = NULL,
                                      source = NA_character_, query = NA_character_) {
  meta <- parsed$meta
  seqs <- parsed$seqs
  keep <- nzchar(seqs)
  if (!all(keep)) {
    vs_warn("Dropped {sum(!keep)} record{?s} with no sequence body.")
    meta <- meta[keep, , drop = FALSE]
    seqs <- seqs[keep]
  }
  if (!nrow(meta)) vs_abort("No usable records after parsing.")

  ## the defline parse is a fallback only: qualifiers win wherever they exist
  fallback <- parse_influenza_defline(
    ifelse(is.na(meta$strain), meta$definition %||% "", meta$strain),
    host_map = host_map, warn = FALSE
  )
  meta$host <- tolower(ifelse(is.na(meta$host_raw), fallback$host, meta$host_raw))
  meta$country <- normalise_country(meta$country_raw)
  meta$country[is.na(meta$country)] <- fallback$country_guess[is.na(meta$country)]
  meta$year[is.na(meta$year)] <- fallback$year[is.na(meta$year)]
  meta$subtype[is.na(meta$subtype)] <- fallback$subtype[is.na(meta$subtype)]
  meta$host_class <- classify_host(meta$host, host_map = host_map)
  meta$label <- meta$accession
  meta$id <- make.unique(meta$accession)

  bin <- as_ape_bin(seqs, molecule)
  names(bin) <- meta$id

  structure(
    list(seqs = bin, meta = meta, molecule = molecule,
         aligned = detect_aligned(bin),
         query = query, source = source),
    class = "vs_sequences"
  )
}

as_ape_bin <- function(seqs, molecule) {
  chars <- lapply(strsplit(toupper(seqs), "", fixed = TRUE), identity)
  if (molecule == "nucleotide") {
    ape::as.DNAbin(chars)
  } else {
    ape::as.AAbin(chars)
  }
}

#' Canonical country names for GenBank `/country` values
#'
#' GenBank writes `"Viet Nam: Hanoi"`; the World Bank and this package write
#' `"Vietnam"`. This strips the sub-national part after the colon and maps the
#' common spelling variants onto the names used by [country_to_iso3()].
#'
#' @param x Character vector of `/country` (or `/geo_loc_name`) values.
#' @return Character vector of canonical country names, `NA` where absent.
#' @export
normalise_country <- function(x) {
  x <- trimws(sub(":.*$", "", as.character(x)))
  x[!nzchar(x)] <- NA_character_
  alias <- c(
    "viet nam" = "Vietnam", "vietnam" = "Vietnam",
    "lao pdr" = "Laos", "lao people's democratic republic" = "Laos", "laos" = "Laos",
    "burma" = "Myanmar", "myanmar" = "Myanmar",
    "korea, south" = "South Korea", "republic of korea" = "South Korea",
    "south korea" = "South Korea", "korea" = "South Korea",
    "russian federation" = "Russia", "russia" = "Russia",
    "usa" = "United States", "united states of america" = "United States",
    "united states" = "United States",
    "uk" = "United Kingdom", "united kingdom" = "United Kingdom",
    "hong kong" = "Hong Kong", "china" = "China", "taiwan" = "Taiwan",
    "timor-leste" = "Timor-Leste", "east timor" = "Timor-Leste"
  )
  key <- tolower(x)
  out <- unname(alias[key])
  ## anything not in the alias table keeps its own spelling, title-cased so that
  ## "cambodia" and "Cambodia" do not become two countries
  keep <- is.na(out) & !is.na(x)
  out[keep] <- title_case(x[keep])
  out
}

title_case <- function(x) {
  gsub("(^|[ -])([a-z])", "\\1\\U\\2", tolower(x), perl = TRUE)
}
