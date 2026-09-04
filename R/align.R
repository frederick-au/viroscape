#' Which alignment backends are available
#'
#' The package deliberately does not hard-depend on any one aligner. MUSCLE via
#' the Bioconductor `msa` package reproduces the original study; `DECIPHER` is a
#' pure-R alternative; an external `mafft`/`muscle` binary is used if present;
#' and already-aligned input can be passed straight through.
#'
#' @return A tibble of backends with an `available` flag.
#' @export
available_aligners <- function() {
  ext <- c(mafft = unname(Sys.which("mafft")), muscle = unname(Sys.which("muscle")))
  tibble::tibble(
    backend = c("msa", "decipher", "mafft", "muscle", "none"),
    available = c(is_installed("msa"), is_installed("DECIPHER"),
                  nzchar(ext[["mafft"]]), nzchar(ext[["muscle"]]), TRUE),
    detail = c("Bioconductor msa (MUSCLE) - matches the published pipeline",
               "Bioconductor DECIPHER::AlignSeqs",
               ext[["mafft"]], ext[["muscle"]],
               "Pass through input that is already aligned")
  )
}

#' Align sequences
#'
#' @param x A `vs_sequences` object.
#' @param method Backend to use. `"auto"` picks the first available of `msa`,
#'   `decipher`, `mafft`, `muscle`, falling back to `"none"` if the input is
#'   already aligned.
#' @param gap_threshold Drop alignment columns that are gaps in more than this
#'   proportion of sequences. `NULL` disables trimming.
#' @param ... Passed to the underlying aligner.
#' @return A `vs_alignment` object wrapping a `phyDat` alignment plus metadata.
#' @export
align_sequences <- function(x, method = c("auto", "msa", "decipher", "mafft",
                                          "muscle", "none"),
                            gap_threshold = 0.9, ...) {
  stopifnot(inherits(x, "vs_sequences"))
  method <- match.arg(method)
  if (method == "auto") method <- pick_aligner(x)
  vs_inform("Aligning {length(x$seqs)} sequence{?s} with backend {.strong {method}}.")

  mat <- switch(method,
    msa = align_with_msa(x, ...),
    decipher = align_with_decipher(x, ...),
    mafft = align_external(x, "mafft", ...),
    muscle = align_external(x, "muscle", ...),
    none = align_passthrough(x)
  )

  if (!is.null(gap_threshold)) mat <- trim_gappy_columns(mat, gap_threshold)

  type <- if (x$molecule == "nucleotide") "DNA" else "AA"
  phy <- phangorn::phyDat(mat, type = type)

  structure(
    list(matrix = mat, phy = phy, meta = x$meta[match(rownames(mat), x$meta$id), , drop = FALSE],
         molecule = x$molecule, method = method, n_sites = ncol(mat)),
    class = "vs_alignment"
  )
}

pick_aligner <- function(x) {
  av <- available_aligners()
  for (b in c("msa", "decipher", "mafft", "muscle")) {
    if (av$available[av$backend == b]) return(b)
  }
  if (isTRUE(x$aligned)) {
    ev <- x$aligned_evidence %||% "unknown"
    if (!ev %in% c("gaps", "matrix")) {
      vs_warn(c(
        "Passing the input through unaligned: no backend is installed and the only evidence of alignment is {.val {ev}}.",
        "i" = "Equal sequence lengths alone do not mean aligned; full-length HA records are routinely all the same length.",
        "i" = "Install a backend, or pass {.code method = \"none\"} explicitly to silence this."
      ))
    }
    return("none")
  }
  vs_abort(c(
    "No alignment backend is available and the input is not aligned.",
    "i" = 'Install one with {.run BiocManager::install("DECIPHER")}, or supply a pre-aligned FASTA.'
  ))
}

align_with_msa <- function(x, ...) {
  need_pkg("msa", "align with MUSCLE")
  need_pkg("Biostrings", "align with MUSCLE")
  chars <- vapply(as.character(x$seqs), paste, character(1), collapse = "")
  ss <- if (x$molecule == "nucleotide") {
    Biostrings::DNAStringSet(chars)
  } else {
    Biostrings::AAStringSet(toupper(chars))
  }
  names(ss) <- names(x$seqs)
  al <- msa::msa(ss, method = "Muscle", ...)
  strings_to_matrix(as.character(methods_as_stringset(al)))
}

methods_as_stringset <- function(al) {
  if (is_installed("Biostrings")) {
    return(methods::as(al, "XStringSet"))
  }
  al
}

align_with_decipher <- function(x, ...) {
  need_pkg("DECIPHER", "align with DECIPHER")
  need_pkg("Biostrings", "align with DECIPHER")
  chars <- vapply(as.character(x$seqs), paste, character(1), collapse = "")
  ss <- if (x$molecule == "nucleotide") {
    Biostrings::DNAStringSet(chars)
  } else {
    Biostrings::AAStringSet(toupper(chars))
  }
  names(ss) <- names(x$seqs)
  al <- DECIPHER::AlignSeqs(ss, verbose = FALSE, ...)
  strings_to_matrix(as.character(al))
}

align_external <- function(x, cmd, extra_args = character(), ...) {
  bin <- unname(Sys.which(cmd))
  if (!nzchar(bin)) vs_abort("Could not find {.strong {cmd}} on the PATH.")
  infile <- tempfile(fileext = ".fasta")
  outfile <- tempfile(fileext = ".fasta")
  write_fasta(x$seqs, infile)
  args <- if (cmd == "mafft") {
    c("--auto", "--quiet", extra_args, shQuote(infile))
  } else {
    c("-align", shQuote(infile), "-output", shQuote(outfile), extra_args)
  }
  if (cmd == "mafft") {
    res <- system2(bin, args, stdout = outfile, stderr = FALSE)
  } else {
    res <- system2(bin, args, stdout = FALSE, stderr = FALSE)
  }
  if (!file.exists(outfile) || file.size(outfile) == 0) {
    vs_abort("{cmd} alignment failed (exit status {res}).")
  }
  aligned <- read_fasta_strings(outfile)
  strings_to_matrix(aligned)
}

align_passthrough <- function(x) {
  lens <- lengths(x$seqs)
  if (length(unique(lens)) != 1L) {
    vs_abort(c(
      "Sequences are not the same length, so they cannot be used unaligned.",
      "i" = "Lengths range from {min(lens)} to {max(lens)}."
    ))
  }
  chars <- vapply(as.character(x$seqs), paste, character(1), collapse = "")
  names(chars) <- names(x$seqs)
  strings_to_matrix(chars)
}

strings_to_matrix <- function(strings) {
  nm <- names(strings)
  split <- strsplit(toupper(unname(strings)), "", fixed = TRUE)
  mat <- do.call(rbind, split)
  rownames(mat) <- nm
  mat
}

write_fasta <- function(seqs, path) {
  chars <- vapply(as.character(seqs), paste, character(1), collapse = "")
  lines <- as.vector(rbind(paste0(">", names(seqs)), toupper(chars)))
  writeLines(lines, path)
  invisible(path)
}

read_fasta_strings <- function(path) {
  lines <- readLines(path, warn = FALSE)
  idx <- grep("^>", lines)
  nm <- sub("^>\\s*", "", lines[idx])
  starts <- idx + 1
  ends <- c(idx[-1] - 1, length(lines))
  out <- vapply(seq_along(idx), function(i) {
    if (starts[i] > ends[i]) return("")
    paste(lines[starts[i]:ends[i]], collapse = "")
  }, character(1))
  stats::setNames(toupper(out), nm)
}

trim_gappy_columns <- function(mat, threshold) {
  is_gap <- matrix(mat %in% c("-", ".", "?"), nrow = nrow(mat))
  gapfrac <- colMeans(is_gap)
  keep <- gapfrac <= threshold
  if (!any(keep)) vs_abort("Gap trimming removed every alignment column.")
  dropped <- sum(!keep)
  if (dropped) vs_inform("Trimmed {dropped} column{?s} with >{threshold*100}% gaps.")
  mat[, keep, drop = FALSE]
}

#' @export
print.vs_alignment <- function(x, ...) {
  cli::cli_h3("vs_alignment")
  cli::cli_text("{nrow(x$matrix)} sequence{?s} x {x$n_sites} site{?s} ({x$molecule}, backend: {x$method})")
  invisible(x)
}
