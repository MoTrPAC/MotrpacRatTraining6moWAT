#' @title Update Gene Ontology term descriptions
#'
#' @description Update entries in the \code{gs_description} column of
#'   \code{\link{msigdbr}} results to use terms from the
#'   \href{http://geneontology.org/}{Gene Ontology Consortium}.
#'
#' @param x object of class \code{data.frame} produced by
#'   \code{link[msigdbr]{msigdbr}} or \code{\link{msigdbr2}} containing columns
#'   \code{gs_description} and \code{gs_subcat}.
#' @param version character; specifies the version of \code{msigdbr} to use.
#'   Defaults to the current version.
#' @param capitalize logical; whether to capitalize the first letter of each
#'   description if the first word does not contain a mix of capital and
#'   lowercase letters. Improves appearance of plots, such as those produced by
#'   \code{\link{enrichmat}}.
#'
#' @returns Object of class \code{data.frame}. The same as \code{x}, but with
#'   updated descriptions.
#'
#' @details This function assumes that the phrase "GO-basic obo file released
#'   on" is present in the MSigDB release notes for that version and is followed
#'   by a date.
#'
#' @md
#'
#' @author Tyler Sagendorf
#'
#' @seealso
#' \href{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Release_Notes}{MSigDB
#' Release Notes}, \href{http://release.geneontology.org/}{Gene Ontology Data
#' Archive}
#'
#' @references Ashburner, M., et al. (2000). Gene ontology: tool for the
#'   unification of biology. The Gene Ontology Consortium. \emph{Nature
#'   genetics, 25}(1), 25--29. \url{https://doi.org/10.1038/75556}
#'
#'   Gene Ontology Consortium (2021). The Gene Ontology resource: enriching a
#'   GOld mine. \emph{Nucleic acids research, 49}(D1), D325--D334.
#'   \url{https://doi.org/10.1093/nar/gkaa1113}
#'
#' @importFrom data.table setDT setDF setkeyv `:=`
#' @importFrom ontologyIndex get_OBO
#' @importFrom utils head packageVersion
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd bfcrpath bfcnew
#'   bfccount
#'
#' @export update_GO_names
#'
#' @examples
#' x <- msigdbr2(species = "Homo sapiens",
#'               genes = "gene_symbol",
#'               gs_subcat = "GO:MF")
#' set.seed(9900)
#' idx <- sample(1:nrow(x), size = 6) # random indices to illustrate changes
#' x$gs_description[idx] # before
#'
#' y <- update_GO_names(x, capitalize = TRUE)
#' y$gs_description[idx] # after


update_GO_names <- function(x,
                            version = packageVersion("msigdbr"),
                            capitalize = FALSE)
{
  go_subcats <- c("GO:MF", "GO:CC", "GO:BP")

  setDT(x)

  if (any(x[["gs_subcat"]] %in% go_subcats)) {

    dir <- tools::R_user_dir(package = "MotrpacRatTraining6moWAT",
                             which = "cache")

    bfc <- BiocFileCache(cache = dir)
    rname <- sprintf("MSigDB_v%s_Release_Notes", version)

    q1 <- bfcquery(bfc, query = rname, field = "rname")

    if (bfccount(q1) == 1L) {
      file <- q1$fpath
    } else {
      message(sprintf("Searching MSigDB %s Release Notes for OBO file date:",
                      version))
      file <- obo_file(version = version)
      bfcadd(bfc, fpath = file, rname = rname, fname = "unique")
    }

    message(paste("Updating GO term descriptions with", file))

    obo_date <- sub(".*org/([^/]+)/.*", "\\1", file)
    obo_rname <- paste0("GO_OBO_data_", obo_date)

    q2 <- bfcquery(bfc, query = obo_rname, field = "rname")

    if (bfccount(q2) == 1L) {
      go.dt <- readRDS(file = bfcrpath(bfc, rnames = obo_rname))
    } else {
      message("Downloading OBO file to cache:")
      go_basic_list <- get_OBO(file,
                               propagate_relationships = "is_a",
                               extract_tags = "minimal")
      go.dt <- as.data.frame(go_basic_list)
      setDT(go.dt)
      go.dt <- go.dt[obsolete != TRUE & grepl("^GO", id), list(id, name)]
      saveRDS(go.dt, file = bfcnew(bfc, rname = obo_rname, fname = "unique"),
              compress = TRUE)
    }

    # Update gs_description column with names from OBO file
    x[go.dt, on = list(gs_exact_source = id), gs_description := i.name]

  } else {
    message("No Gene Ontology term descriptions to modify.")
  }

  if (capitalize) {
    message("Capitalizing entries of gs_description column.")
    x[, `:=` (gs_description = cap_names(gs_description))]
  }

  setDF(x)
  return(x)
}




## Helper functions ------------------------------------------------------------

# Get the path to the appropriate OBO file.
obo_file <- function(version = packageVersion("msigdbr"))
{
  # MSigDB release notes for appropriate version
  url <- sprintf(file.path("https://software.broadinstitute.org/cancer",
                           "software/gsea/wiki/index.php",
                           "MSigDB_v%s_Release_Notes"), version)

  x <- readLines(url)
  x <- paste(x, collapse = "")
  x <- gsub("\\\\n", "", x)

  phrase <- "GO-basic obo file released on"

  if (!grepl(phrase, x)) {
    stop(sprintf("Phrase '%s' not found in %s", phrase, url))
  }

  obo_date <- sub(sprintf(".*%s ([^ ]+).*", phrase), "\\1", x)
  obo_file <- sprintf("http://release.geneontology.org/%s/ontology/go-basic.obo",
                      obo_date)

  return(obo_file)
}


# Capitalize first letter of first word unless there is
# already a mix of uppercase and lowercase letters.
# x is a character vector.
cap_names <- function(x)
{
  first_word <- sub(" .*", "", x)
  idx <- first_word == tolower(first_word)

  x[idx] <- sub("(.)", "\\U\\1", x[idx], perl = TRUE)
  return(x)
}

