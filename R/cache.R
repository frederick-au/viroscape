#' Resolve the on-disk cache directory
#'
#' Live retrieval from NCBI and the World Bank is cached so that a pipeline run
#' is reproducible offline and so that repeated Shiny interactions do not
#' re-hit the network. Set the `VIROSCAPE_CACHE` environment variable or the
#' `viroscape.cache` option to override the default.
#'
#' @param path Optional explicit path. If supplied it is created and returned.
#' @return Path to the cache directory, created if necessary.
#' @export
vs_cache_dir <- function(path = NULL) {
  path <- path %||% getOption("viroscape.cache") %||%
    Sys.getenv("VIROSCAPE_CACHE", unset = "") 
  if (identical(path, "")) {
    path <- tools::R_user_dir("viroscape", which = "cache")
  }
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

#' Deterministic cache key for a set of inputs
#'
#' Exported because [with_cache()] takes a key and there is otherwise no
#' supported way to make one.
#'
#' @param ... Any objects. The key is a hash of their serialised form, so it is
#'   stable across sessions and sensitive to every argument.
#' @return A single character string.
#' @export
vs_cache_key <- function(...) {
  payload <- paste(vapply(list(...), function(x) {
    paste(utils::capture.output(utils::str(x, max.level = 3)), collapse = "|")
  }, character(1)), collapse = "||")
  bytes <- as.numeric(charToRaw(payload))
  ## two independent polynomial rolling hashes, kept inside double precision
  h1 <- 0; h2 <- 0
  for (b in bytes) {
    h1 <- (h1 * 31 + b) %% 4294967291
    h2 <- (h2 * 131 + b) %% 4294967279
  }
  paste0(format(as.hexmode(as.integer(h1 %% 2147483647)), width = 8),
         format(as.hexmode(as.integer(h2 %% 2147483647)), width = 8))
}

#' Cache the value of an expression on disk
#'
#' @param key Cache key, typically from [vs_cache_key()].
#' @param expr Expression to evaluate if the cache misses.
#' @param cache_dir Directory to store the cache in.
#' @param refresh If `TRUE`, ignore any existing cache entry and re-evaluate.
#' @param quiet Suppress the cache-hit message.
#' @return The value of `expr`, from cache when available.
#' @export
with_cache <- function(key, expr, cache_dir = vs_cache_dir(),
                       refresh = FALSE, quiet = FALSE) {
  file <- file.path(cache_dir, paste0(key, ".rds"))
  if (!refresh && file.exists(file)) {
    if (!quiet) vs_inform("Using cached result ({.path {basename(file)}}).")
    return(readRDS(file))
  }
  value <- force(expr)
  saveRDS(value, file)
  value
}

#' List or clear cached downloads
#'
#' @param cache_dir Cache directory.
#' @return For `vs_cache_list()`, a tibble of cache entries. `vs_cache_clear()`
#'   returns the number of files removed, invisibly.
#' @export
vs_cache_list <- function(cache_dir = vs_cache_dir()) {
  f <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  tibble::tibble(
    file = basename(f),
    size_kb = round(file.size(f) / 1024, 1),
    modified = file.mtime(f)
  )
}

#' @rdname vs_cache_list
#' @export
vs_cache_clear <- function(cache_dir = vs_cache_dir()) {
  f <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  n <- sum(file.remove(f))
  vs_inform("Removed {n} cached file{?s}.")
  invisible(n)
}
