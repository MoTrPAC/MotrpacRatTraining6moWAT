#' @title Update Gene Ontology term descriptions
#'
#' @description Update entries in the \code{gs_description} column of
#'   \code{\link{msigdbr}} results to use terms from the
#'   \href{http://geneontology.org/}{Gene Ontology Consortium}.
#'
#' @param x object of class \code{data.frame} produced by
#'   \code{link[msigdbr]{msigdbr}} or \code{\link{msigdbr2}} containing columns
#'   \code{gs_description} and \code{gs_subcat}.
#' @param version character string specifying the version of \code{msigdbr} to
#'   use. Defaults to the current version.
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
#' @importFrom data.table setDT setDF setkeyv `:=`
#' @importFrom ontologyIndex get_OBO
#' @importFrom utils head packageVersion
#'
#' @export update_GO_names
#'
#' @author Tyler Sagendorf
#'
#' @seealso
#' \href{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Release_Notes}{MSigDB
#' Release Notes} \href{http://release.geneontology.org/}{Gene Ontology Data
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
#' @examples
#' x <- msigdbr2(species = "Homo sapiens",
#'               genes = "gene_symbol",
#'               gs_subcat = "GO:MF")
#' tail(x$gs_description) # before
#'
#' x <- update_GO_names(x)
#' tail(x$gs_description) # after

update_GO_names <- function(x,
                            version = packageVersion("msigdbr"),
                            capitalize = FALSE)
{
  go_subcats <- c("GO:MF", "GO:CC", "GO:BP")

  # setDT(x, key = "gs_exact_source")
  setDT(x)

  if (any(x[["gs_subcat"]] %in% go_subcats)) {
    file <- obo_file(version = version)
    message(paste("Updating GO term descriptions with", file))

    go_basic_list <- get_OBO(file,
                             propagate_relationships = "is_a",
                             extract_tags = "minimal")
    go.dt <- as.data.frame(go_basic_list)
    setDT(go.dt)
    go.dt <- go.dt[obsolete != TRUE & grepl("^GO", id), list(id, name)]
    # setkeyv(go.dt, cols = "id")

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
obo_file <- function(version = packageVersion("msigdbr")) {
  # MSigDB release notes for appropriate version
  path <- sprintf(file.path("https://software.broadinstitute.org/cancer",
                            "software/gsea/wiki/index.php",
                            "MSigDB_v%s_Release_Notes"), version)
  x <- readLines(path)
  x <- paste(x, collapse = "")
  x <- gsub("\\\\n", "", x)

  phrase <- "GO-basic obo file released on"

  if (!grepl(phrase, x)) {
    stop(sprintf("Phrase '%s' not found in %s", phrase, path))
  }

  obo_date <- sub(sprintf(".*%s ([^ ]+).*", phrase), "\\1", x)
  obo_file <- sprintf("http://release.geneontology.org/%s/ontology/go-basic.obo",
                      obo_date)
  return(obo_file)
}


# Capitalize first letter of first word unless there is
# already a mix of uppercase and lowercase letters.
# x is a character vector.
cap_names <- function(x) {
  first_word <- sub(" .*", "", x)
  idx <- first_word == tolower(first_word)

  x[idx] <- sub("(.)", "\\U\\1", x[idx], perl = TRUE)
  return(x)
}

