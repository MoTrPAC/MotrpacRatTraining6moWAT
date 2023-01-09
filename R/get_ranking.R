#' @title Calculate ranking metric values by contrast
#'
#' @description Calculate ranking metric values for \code{\link{fgsea2}}.
#'
#' @param x object of class \code{data.frame}. Typically results from
#'   \code{\link{limma_full}}.
#' @param genes string; the name of a column in \code{x}. A ranking metric value
#'   will be calculated for each unique entry. Can be gene symbols, Entrez IDs,
#'   or Ensembl genes.
#' @param metric character; specifies how the ranking metric will be calculated
#'   from \code{x}. Can be the name of a single column.
#' @param contrast_column character; the name of a column in \code{x}. Entries
#'   in this column will be set as the names of the output list.
#'
#' @details The -log\eqn{_{10}}-transformed p-value signed by the fold-change is
#'   similar to using t-statistics, but provides better separation between the
#'   extreme values in the tails and those close to 0 that are not as
#'   interesting. The log\eqn{_2} fold-change should not be used as the ranking
#'   metric, as it does not capture variability in the measurements.
#'
#' @returns A named vector (if no "contrast" column) or list of sorted ranking
#'   metric values for each contrast in \code{x}.
#'
#' @importFrom data.table `.SD` `:=` setDT
#' @importFrom tibble deframe
#'
#' @export get_ranking

get_ranking <- function(x,
                       genes,
                       metric = "-log10(P.Value)*sign(logFC)",
                       contrast_column = "contrast")
{
  stopifnot(is.data.frame(x))

  if (!(genes %in% colnames(x))) {
    stop(sprintf("'%s' is not a column in x.", genes))
  }
  if (!(contrast_column %in% colnames(x))) {
    warning(sprintf("'%s' is not a column in x. Output will be a vector.",
                    contrast_column))
    stats <- .get_ranking(x, genes = genes, metric = metric)
  } else {
    x <- split.data.frame(x = x, f = x[[contrast_column]], drop = TRUE)
    stats <- lapply(x, .get_ranking, genes = genes, metric = metric)
  }
  return(stats)
}


.get_ranking <- function(x, genes, metric) {
  setDT(x)
  # Subset to non-missing genes
  x <- subset(x, subset = !is.na(get(genes)))
  if (nrow(x) == 0) {
    stop(sprintf("All entries in '%s' are missing.", genes))
  }
  # Calculate ranking metric for each gene
  x[, rank.metric := eval(parse(text = metric))]
  cols <- c(genes, "rank.metric")
  x <- unique(x[, cols, with = FALSE])
  x <- x[, lapply(.SD, mean, na.rm = TRUE), by = genes] # mean by gene
  stats <- sort(deframe(x), decreasing = TRUE) # sorting doesn't actually matter
  names(stats) <- as.character(names(stats))
  return(stats)
}

