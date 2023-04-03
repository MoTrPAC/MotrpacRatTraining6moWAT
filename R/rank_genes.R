#' @title Calculate GSEA ranking metric values
#'
#' @description Calculate GSEA ranking metric values for \code{\link{fgsea2}}.
#'
#' @param x object of class \code{data.frame}. Typically results from
#'   \code{\link{limma_full}}.
#' @param genes character; the name of a column in \code{x}. A ranking metric
#'   value will be calculated for each unique entry. Usually gene symbols,
#'   Entrez IDs, or Ensembl genes to work with the output of \code{msigdbr2}.
#' @param metric character; specifies how the ranking metric will be calculated
#'   from \code{x}. Can be the name of a single column.
#' @param contrast_column character; the name of a column in \code{x}. Entries
#'   in this column will be set as the names of the output list.
#'
#' @details The -log\eqn{_{10}}-transformed p-value signed by the fold-change is
#'   similar to the t-statistic, but more clearly separates extreme values in
#'   the tails from those close to 0 that are not as interesting.
#'
#' @returns A named vector (if no "contrast" column) or list of named vectors of
#'   sorted ranking metric values for each contrast in \code{x}.
#'
#' @importFrom data.table `.SD` `:=` setDT
#' @importFrom tibble deframe
#'
#' @export rank_genes

rank_genes <- function(x,
                       genes,
                       metric = "-log10(P.Value)*sign(logFC)",
                       contrast_column = "contrast")
{
  if (!(genes %in% colnames(x))) {
    stop(sprintf("'%s' is not a column in x.", genes))
  }
  if (!(contrast_column %in% colnames(x))) {
    warning(sprintf("'%s' is not a column in x. Output will be a vector.",
                    contrast_column))
    stats <- .rank_genes(x, genes = genes, metric = metric)
  } else {
    x <- split.data.frame(x = x, f = x[[contrast_column]], drop = TRUE)
    stats <- lapply(x, .rank_genes, genes = genes, metric = metric)
  }
  return(stats)
}


.rank_genes <- function(x, genes, metric) {
  setDT(x)
  # Subset to non-missing genes
  x <- subset(x, subset = !is.na(get(genes)))
  if (nrow(x) == 0) {
    stop(sprintf("All entries in '%s' are missing.", genes))
  }
  # Calculate ranking metric for each gene
  x[, rank_metric := eval(parse(text = metric))]
  cols <- c(genes, "rank_metric")
  x <- unique(x[, cols, with = FALSE])
  x <- x[, lapply(.SD, mean, na.rm = TRUE), by = genes] # mean by gene
  stats <- sort(deframe(x), decreasing = TRUE) # sorting doesn't actually matter
  names(stats) <- as.character(names(stats))
  return(stats)
}

