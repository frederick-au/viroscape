#' Default rules mapping raw host labels to ecological classes
#'
#' An ordered list of regular expressions. Each element is one ecological
#' class; the first pattern that matches a host string wins, so the order is
#' part of the specification and not incidental.
#'
#' Two decisions are worth stating explicitly, because they are the ones that
#' change results.
#'
#' **Waterfowl default to domestic.** A bare `duck` or `goose` is classified
#' `domestic_poultry`, not `wild_bird`. In Southeast Asian H5N1 the free-grazing
#' domestic duck in rice-paddy systems is the established risk factor (Gilbert
#' et al. 2008, PNAS 105:4769), and it is the commonest unqualified waterfowl
#' label in GenBank. A record is only `wild_bird` if it carries an explicit wild
#' marker (`wild`, `migratory`, `feral`) or names a species that is unambiguously
#' wild. This inverts the behaviour of viroscape 0.1.0, which classified
#' `domestic duck` as `wild_bird`.
#'
#' **Latin binomials are recognised.** GenBank's `/host` qualifier frequently
#' carries a binomial rather than a common name, so `Gallus gallus`,
#' `Cairina moschata` and `Anas platyrhynchos` are matched directly.
#'
#' `Anas platyrhynchos` is genuinely ambiguous - the mallard is wild, but
#' domestic ducks are derived from it and some submitters use the binomial for
#' domestic birds. It is treated as wild here. If that is wrong for your data,
#' pass your own rules through the `host_map` argument of [classify_host()].
#'
#' **Bare genus names are `unknown`, not domestic.** `Anas sp.`, `Anser sp.`,
#' `Anatidae` and an unqualified `waterfowl` say only "some duck or goose,
#' species not determined". They carry less information than a vernacular
#' `duck`, not more, and they are the labels wild-bird surveillance tends to
#' produce - so classifying them as domestic would give the vaguer label a more
#' confident answer than the specific one, in the wrong direction. They are
#' excluded from host-stratified analyses instead of biasing them.
#'
#' @return A named list of regular expressions, in application order.
#' @seealso [classify_host()]
#' @export
vs_host_rules <- function() {
  list(
    ## a genus with no species epithet identifies a bird only as "some duck" or
    ## "some goose"; it is less specific than the vernacular label, so it does
    ## not inherit the vernacular label's domestic default
    unknown = paste0("^anas( ?sp?p?\\.?)?$|^anser( ?sp?p?\\.?)?$|^anatidae$|",
                     "^water ?(fowl|bird)s?$|^avian$|^bird$|^unknown$|^ *$"),
    human = "human|homo sapiens|h\\. sapiens|patient",
    ## \bwater\b, not "water": the unanchored form swallowed "waterfowl" and
    ## classified wild waterfowl as an environmental sample
    environment = "environment|\\bwater\\b|waste ?water|\\blake\\b|\\bpond\\b|market|sediment|soil|swab|air sample|faeces|feces",
    swine = "swine|\\bpig\\b|piglet|porcine|sus scrofa|boar",
    mammal = paste(
      "cat\\b|feline|felis|dog\\b|canine|canis|tiger|panthera|leopard|civet|mink",
      "neovison|mustela|ferret|fox\\b|vulpes|seal\\b|phoca|sea lion|otter|marten",
      "badger|raccoon|skunk|bear\\b|ursus|cattle|bovine|bos taurus|dairy|goat",
      "capra|sheep|ovis|horse|equine|equus|mouse|mus musculus|mice|\\brat\\b|rattus",
      sep = "|"),
    ## explicit wild markers, then unambiguously wild taxa
    wild_bird = paste(
      "\\bwild\\b|migrat|feral",
      "mallard|anas platyrhynchos|teal|anas crecca|pintail|anas acuta|wigeon",
      "widgeon|shoveler|gadwall|garganey|swan|cygnus|gull|larus|tern\\b|sterna",
      "shorebird|wader|sandpiper|calidris|egret|heron|ardea|crow|corvus|magpie",
      "falcon|falco|eagle|haliaeetus|aquila|owl\\b|sparrow|starling|sturnus",
      "grebe|cormorant|phalacrocorax|stork|ciconia|crane\\b|grus|pigeon|columba",
      sep = "|"),
    ## everything else kept in poultry, including bare duck / goose
    domestic_poultry = paste(
      "chicken|gallus|turkey|meleagris|quail|coturnix|guinea ?fowl|numida",
      "broiler|layer|poultry|domestic|backyard|muscovy|cairina|silkie|partridge",
      "pheasant|peafowl|ostrich|emu\\b|duck|duckling|goose|geese|gosling|anas|anser",
      sep = "|")
  )
}

#' Group raw influenza host labels into ecological classes
#'
#' @param host Character vector of host labels, from strain names or from a
#'   GenBank `/host` qualifier.
#' @param host_map Optional named list of regular expressions overriding
#'   [vs_host_rules()]. Names become the class labels; the first match wins.
#'   Anything unmatched is `unknown`.
#' @return Factor with the levels of `host_map` plus `unknown`.
#' @examples
#' classify_host(c("domestic duck", "wild duck", "Gallus gallus", "chicken"))
#' @export
classify_host <- function(host, host_map = NULL) {
  rules <- host_map %||% vs_host_rules()
  if (!is.list(rules) || is.null(names(rules)) || anyNA(names(rules)) ||
      !all(nzchar(names(rules)))) {
    vs_abort("{.arg host_map} must be a named list of regular expressions.")
  }
  h <- tolower(trimws(as.character(host)))
  out <- rep(NA_character_, length(h))
  for (cls in names(rules)) {
    hit <- is.na(out) & !is.na(h) & grepl(rules[[cls]], h, perl = TRUE)
    out[hit] <- cls
  }
  out[is.na(out)] <- "unknown"
  ## keep the historical level order so existing plots and factors are unchanged
  canonical <- c("wild_bird", "domestic_poultry", "swine", "mammal", "human",
                 "environment", "unknown")
  levs <- unique(c(intersect(canonical, c(names(rules), "unknown")),
                   setdiff(names(rules), canonical), "unknown"))
  factor(out, levels = levs)
}
