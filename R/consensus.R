#' Build a consensus sequence from an alignment
#'
#' The response variable in this framework is divergence from a regional
#' consensus, so the consensus is the reference point the whole analysis hangs
#' on. By default it is built from every sequence in the alignment by simple
#' majority rule, matching the published method. Building it from a subset (for
#' example wild bird isolates only) is often more interpretable, which is what
#' `subset` is for.
#'
#' @param alignment A `vs_alignment` object.
#' @param method `"majority"` takes the most common non-gap character at each
#'   site; `"threshold"` requires the plurality character to reach `threshold`
#'   and emits `ambiguity` otherwise.
#' @param threshold Proportion required by the `"threshold"` method.
#' @param ambiguity Character emitted where no residue meets the threshold.
#' @param subset Optional logical vector or sequence ids to build from.
#' @param ignore_gaps Exclude gaps when counting characters at each site.
#' @return A `vs_consensus` object: a character vector of the consensus with
#'   per-site support attached.
#' @export
build_consensus <- function(alignment, method = c("majority", "threshold"),
                            threshold = 0.5, ambiguity = NULL,
                            subset = NULL, ignore_gaps = TRUE) {
  stopifnot(inherits(alignment, "vs_alignment"))
  method <- match.arg(method)
  ambiguity <- ambiguity %||% if (alignment$molecule == "nucleotide") "N" else "X"

  mat <- alignment$matrix
  if (!is.null(subset)) {
    if (is.logical(subset)) mat <- mat[subset, , drop = FALSE]
    else mat <- mat[rownames(mat) %in% subset, , drop = FALSE]
    if (!nrow(mat)) vs_abort("`subset` selected no sequences.")
  }

  gaps <- c("-", ".", "?")
  res <- apply(mat, 2, function(col) {
    if (ignore_gaps) col <- col[!col %in% gaps]
    if (!length(col)) return(c(ambiguity, 0))
    tab <- sort(table(col), decreasing = TRUE)
    support <- tab[[1]] / sum(tab)
    chr <- names(tab)[1]
    if (method == "threshold" && support < threshold) chr <- ambiguity
    c(chr, support)
  })

  structure(
    list(sequence = unname(res[1, ]),
         support = as.numeric(res[2, ]),
         n_sequences = nrow(mat),
         molecule = alignment$molecule,
         method = method),
    class = "vs_consensus"
  )
}

#' @export
print.vs_consensus <- function(x, ...) {
  cli::cli_h3("vs_consensus")
  cli::cli_text("{length(x$sequence)} site{?s} from {x$n_sequences} sequence{?s} ({x$method} rule)")
  cli::cli_text("Mean per-site support: {round(mean(x$support), 3)}")
  preview <- paste(utils::head(x$sequence, 60), collapse = "")
  cli::cli_text("{.val {preview}}...")
  invisible(x)
}

#' @export
as.character.vs_consensus <- function(x, ...) paste(x$sequence, collapse = "")
