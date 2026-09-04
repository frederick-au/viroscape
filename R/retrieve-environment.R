#' Environmental indicators available from the World Bank
#'
#' The eight indicators used in the source study, plus a handful of further
#' candidates. Automating retrieval is what makes it practical to extend the
#' analysis beyond five countries, which is the main limitation of the original
#' design: reliable random effects for country need roughly ten groups.
#'
#' @return A tibble of indicator codes and labels. The `default` column marks
#'   the eight indicators [fetch_environment()] retrieves when `variables` is
#'   not supplied, and is what the Shiny app pre-selects.
#' @export
wb_indicators <- function() {
  tibble::tribble(
    ~variable,               ~code,                  ~default, ~label,
    "agricultural_land_pct", "AG.LND.AGRI.ZS",       TRUE,     "Agricultural land (% of land area)",
    "forest_area_pct",       "AG.LND.FRST.ZS",       TRUE,     "Forest area (% of land area)",
    "population",            "SP.POP.TOTL",          TRUE,     "Population, total",
    "urban_population_pct",  "SP.URB.TOTL.IN.ZS",    TRUE,     "Urban population (% of total)",
    "fertilizer_kg_ha",      "AG.CON.FERT.ZS",       TRUE,     "Fertilizer consumption (kg/ha arable land)",
    "livestock_index",       "AG.PRD.LVSK.XD",       TRUE,     "Livestock production index (2014-2016 = 100)",
    "crop_index",            "AG.PRD.CROP.XD",       FALSE,    "Crop production index (2014-2016 = 100)",
    "cereal_yield",          "AG.YLD.CREL.KG",       TRUE,     "Cereal yield (kg per hectare)",
    "ghg_emissions",         "EN.GHG.ALL.MT.CE.AR5", TRUE,     "Total greenhouse gas emissions (Mt CO2e)",
    "rural_population_pct",  "SP.RUR.TOTL.ZS",       FALSE,    "Rural population (% of total)",
    "population_density",    "EN.POP.DNST",          FALSE,    "Population density (people per sq km)",
    "arable_land_pct",       "AG.LND.ARBL.ZS",       FALSE,    "Arable land (% of land area)",
    "gdp_per_capita",        "NY.GDP.PCAP.CD",       FALSE,    "GDP per capita (current US$)"
  )
}

#' Fetch environmental indicators from the World Bank API
#'
#' @param countries Character vector of country names (mapped via
#'   [country_to_iso3()]) or ISO3 codes directly.
#' @param years Length-2 inclusive year range.
#' @param variables Which of [wb_indicators()] to retrieve. Defaults to the
#'   eight candidates used in the source study.
#' @param cache_dir Cache directory.
#' @param refresh Force a fresh download.
#' @return A tibble in wide form: one row per country-year.
#' @export
fetch_environment <- function(countries, years = c(2003, 2022),
                              variables = NULL,
                              cache_dir = vs_cache_dir(), refresh = FALSE) {
  ind <- wb_indicators()
  variables <- variables %||% ind$variable[ind$default]
  if (!length(variables)) vs_abort("No indicators requested.")
  unknown <- setdiff(variables, ind$variable)
  if (length(unknown)) vs_abort("Unknown indicator{?s}: {.val {unknown}}.")
  ind <- ind[ind$variable %in% variables, ]

  iso <- ifelse(grepl("^[A-Z]{3}$", countries), countries, country_to_iso3(countries))
  if (anyNA(iso)) vs_abort("Could not resolve ISO3 codes for all countries.")
  names(iso) <- countries

  key <- vs_cache_key("worldbank", sort(iso), years, sort(variables))
  with_cache(key, cache_dir = cache_dir, refresh = refresh, expr = {
    vs_inform("Fetching {nrow(ind)} indicator{?s} for {length(iso)} countr{?y/ies} from the World Bank.")
    pieces <- lapply(seq_len(nrow(ind)), function(i) {
      raw <- wb_get(paste(iso, collapse = ";"), ind$code[i], years)
      if (!nrow(raw)) {
        vs_warn("No data returned for {.field {ind$variable[i]}}.")
        return(NULL)
      }
      raw$variable <- ind$variable[i]
      raw
    })
    long <- dplyr::bind_rows(pieces)
    if (!nrow(long)) vs_abort("The World Bank API returned no data.")
    wide <- tidyr::pivot_wider(long, id_cols = c("iso3", "country", "year"),
                              names_from = "variable", values_from = "value")
    wide$country <- names(iso)[match(wide$iso3, iso)]
    dplyr::arrange(wide, country, year)
  })
}

wb_get <- function(iso_string, code, years, tries = 3L, pause = 1) {
  url <- sprintf(
    "https://api.worldbank.org/v2/country/%s/indicator/%s?date=%d:%d&format=json&per_page=20000",
    iso_string, code, years[1], years[2]
  )
  ## A single readLines() with no retry turned every transient network blip into
  ## an aborted run, and reported it as "cannot open the connection", which says
  ## nothing about what to do next.
  txt <- NULL
  last <- NULL
  for (attempt in seq_len(tries)) {
    txt <- tryCatch(
      suppressWarnings(paste(readLines(url, warn = FALSE), collapse = "")),
      error = function(e) {
        last <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(txt) && nzchar(txt)) break
    if (attempt < tries) {
      vs_inform("World Bank request failed ({attempt}/{tries}); retrying in {pause * attempt}s.")
      Sys.sleep(pause * attempt)
    }
  }
  if (is.null(txt) || !nzchar(txt)) {
    vs_abort(c(
      "World Bank request failed after {tries} attempt{?s}: {last %||% 'empty response'}",
      "i" = "Check network access - a proxy or firewall will produce exactly this error.",
      "i" = "The timeout is {getOption('timeout')}s; raise it with {.code options(timeout = 120)} on a slow link.",
      "i" = "To work offline, download the indicators yourself and read them with {.fn read_environment}.",
      "x" = "URL: {url}"
    ))
  }
  parsed <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  if (length(parsed) < 2 || is.null(parsed[[2]]) || !length(parsed[[2]])) {
    return(tibble::tibble())
  }
  d <- parsed[[2]]
  tibble::tibble(
    iso3 = d$countryiso3code,
    country = d$country$value,
    year = as.integer(d$date),
    value = as.numeric(d$value)
  ) |> dplyr::filter(!is.na(value))
}

#' Read environmental data from a local file
#'
#' Use this when you have compiled predictors yourself, for example tree cover
#' loss from Global Forest Watch or land cover from the ESA CCI, neither of
#' which has an open API. Subnational data is supported: supply an `admin1`
#' column and pass `by = c("country", "admin1", "year")` to [join_environment()].
#'
#' @param path Path to a CSV.
#' @param country_col,year_col Names of the country and year columns.
#' @param ... Passed to [utils::read.csv()].
#' @return A tibble with `country` and `year` normalised.
#' @export
read_environment <- function(path, country_col = "country", year_col = "year", ...) {
  if (!file.exists(path)) vs_abort("File {.path {path}} does not exist.")
  d <- tibble::as_tibble(utils::read.csv(path, stringsAsFactors = FALSE, ...))
  if (!country_col %in% names(d)) vs_abort("Column {.field {country_col}} not found.")
  if (!year_col %in% names(d)) vs_abort("Column {.field {year_col}} not found.")
  names(d)[names(d) == country_col] <- "country"
  names(d)[names(d) == year_col] <- "year"
  d$year <- as.integer(d$year)
  d$country <- trimws(as.character(d$country))
  d
}
