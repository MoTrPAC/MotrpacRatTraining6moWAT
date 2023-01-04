#' @title Get gene sets / pathways from MSigDB
#'
#' @description Wrapper around \code{\link[msigdbr]{msigdbr}}. Reformats pathway
#'   data for use with \code{\link{fgsea2}}.
#'
#' @param species character; scientific or common name of species. Default is
#'   "Homo sapiens". \code{\link[msigdbr]{msigdbr_species}} displays all
#'   possible choices.
#' @param genes character; type of gene ID to fetch. Options are "gene_symbol",
#'   "entrez_gene", and "ensembl_gene".
#' @param gs_subcat character vector of one or more MSigDB subcategories. See
#'   \code{\link[msigdbr]{msigdbr_collections}} for details.
#' @inheritParams update_GO_names
#'
#' @returns A \code{\link[data.table]{data.table}} with columns
#'   "gs_exact_source" (term ID), "gs_subcat" (subcategory), "gs_description"
#'   (term description), and a list column of gene IDs specified by
#'   \code{genes}.
#'
#' @importFrom data.table rbindlist `.SD` `:=` setnames setDF
#' @importFrom msigdbr msigdbr_collections msigdbr
#' @importFrom ontologyIndex get_OBO
#'
#' @export msigdbr2
#'
#' @references Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M.,
#'   Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database
#'   (MSigDB) hallmark gene set collection. \emph{Cell systems, 1}(6), 417--425.
#'   \url{https://doi.org/10.1016/j.cels.2015.12.004}
#'
#'   Dolgalev, I. (2022). msigdbr: MSigDB Gene Sets for Multiple Organisms in a
#'   Tidy Data Format. R package version 7.5.1,
#'   \url{https:://igordot.github.io/msigdbr}
#'
#' @examples
#' x <- msigdbr2(species = "rat",
#'               genes = "gene_symbol",
#'               gs_subcat = c("CP:REACTOME", "GO:MF"))

msigdbr2 <- function(species = "Homo sapiens",
                     genes = c("gene_symbol", "entrez_gene", "ensembl_gene"),
                     gs_subcat,
                     capitalize = FALSE)
{
  mcol <- msigdbr_collections()

  # Check that all gs_subcat are valid
  for (i in seq_along(gs_subcat)) {
    if (!(gs_subcat[i] %in% mcol$gs_subcat)) {
      stop(paste0(gs_subcat[i],
                  " is not a valid subcategory. See ?msigdbr_collections."))
    }
  }

  # Get categories
  gs_cat <- mcol[mcol$gs_subcat %in% gs_subcat, c("gs_subcat", "gs_cat")]
  gs_cat <- deframe(gs_cat)[gs_subcat]

  # Fetch MSigDB data
  paths <- lapply(seq_along(gs_cat), function(i) {
    msigdbr(species = species,
            category = gs_cat[i],
            subcategory = gs_subcat[i])
  })
  paths <- rbindlist(paths)
  setDT(paths)

  cols <- c("gs_subcat", "gs_exact_source", "gs_description", genes)
  paths <- unique(paths[, cols, with = FALSE])
  paths <- paths[, lapply(.SD, list),
                 by = c("gs_subcat", "gs_exact_source", "gs_description")]

  # Update GO descriptions
  paths <- update_GO_names(paths, capitalize = capitalize)

  setDF(paths)
  return(paths)
}

