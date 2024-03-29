#' @title Over-representation analysis
#'
#' @description Custom wrapper for \code{\link[fgsea]{fora}}.
#'
#' @param pathways output of \code{link{msigdbr2}}.
#' @param genes named list of genes. These are usually clusters of related
#'   genes.
#' @param universe character vector of background genes. Usually, this is all
#'   unique elements of `genes`.
#' @param gene_column character; the name of the column in pathways containing
#'   the elements of each gene set.
#' @param minSize integer; minimum size of gene sets allowed for testing.
#' @param maxSize integer; maximum size of gene sets allowed for testing.
#' @param adjust.method character string; p-value correction method. Either
#'   "scale" (default), which will compute scaled p-values or one of
#'   \code{\link[stats]{p.adjust.methods}}.
#' @param adjust.globally logical; should p-values from all contrasts be
#'   adjusted together using \code{adjust.method}? Set to \code{FALSE} if the
#'   contrasts being tested are not closely related.
#'
#' @details If \code{adjust.method = "scale"}, the function will calculate the
#'   maximum overlap for each combination of \code{names(genes)} and (if
#'   \code{adjust.globally = TRUE}), the entries of the "gs_subcat" column of
#'   \code{pathways}. For each gene set, it will then calculate the ratio of the
#'   size of that set's overlap to one plus the maximum overlap. The addition of
#'   1 in the denominator is to penalize small maximum overlaps. The
#'   \code{log10(pval)} is multiplied by this overlap ratio and then
#'   back-transformed to obtain p-values adjusted by how well they describe each
#'   group defined by \code{genes}.
#'
#' @seealso \code{\link[fgsea]{fora}}
#'
#' @importFrom stats p.adjust p.adjust.methods
#' @importFrom data.table setDT rbindlist setorderv `:=`
#' @importFrom tibble deframe
#' @importFrom fgsea fora
#'
#' @export fora2
#'
#' @author Tyler Sagendorf

fora2 <- function(pathways,
                  genes,
                  universe,
                  gene_column = "entrez_gene",
                  minSize = 1,
                  maxSize = length(universe) - 1,
                  adjust.method = "scale",
                  adjust.globally = TRUE)
{
  # Check input
  adjust.method <- match.arg(adjust.method,
                             choices = c("scale", p.adjust.methods))

  # List of pathways to test
  setDT(pathways)
  paths <- deframe(pathways[, list(gs_exact_source, get(gene_column))])

  # Run FGSEA on each contrast
  res <- lapply(names(genes), function(group_i) {
    res_i <- fora(pathways = paths,
                  genes = genes[[group_i]],
                  universe = universe,
                  minSize = minSize,
                  maxSize = maxSize)
    res_i[["module"]] <- group_i
    return(res_i)
  })
  res <- rbindlist(res)
  res <- res[res[["overlap"]] != 0]
  setDT(res)

  pathways <- pathways[, list(gs_subcat, gs_exact_source, gs_description)]

  res <- merge(res, pathways, sort = FALSE,
               by.x = "pathway", by.y = "gs_exact_source")

  # p-value adjustment
  by <- "module"
  if (!adjust.globally) {
    by <- c(by, "gs_subcat")
  }

  # Transform p-values to account for overlap ratio
  if (adjust.method == "scale") {
    res[, maxOverlap := max(overlap), by = by]
    res[, overlapRatio := overlap / (maxOverlap + 1)]
    res[, padj := 10 ^ (log10(pval) * overlapRatio)]
  } else {
    res[, padj := p.adjust(pval, method = adjust.method), by = by]
  }
  res[, module := factor(module, levels = unique(module))]

  setorderv(res, cols = c(by, "padj"))


  return(res)
}

