#' @title Fast Gene Set Enrichment Analysis
#'
#' @description Custom wrapper for \code{\link[fgsea]{fgsea}} for use with
#'   output of \code{\link{msigdbr2}}.
#'
#' @param pathways output of \code{\link{msigdbr2}}.
#' @param stats output of \code{\link{get_ranking}}.
#' @param gene_column character string; the name of a column in \code{pathways}
#'   containing the genes in each pathway.
#' @param adjust.method character string; the p-value correction method. Can be
#'   abbreviated. See \code{\link[stats]{p.adjust.methods}}.
#' @param adjust.globally logical; should p-values from all contrasts be
#'   adjusted together using `adjust.method`? Set to `FALSE` if the contrasts
#'   being tested are not closely related.
#' @param seed numeric or `NULL`; passed to `set.seed`.
#' @param ... additional arguments passed to
#'   \code{\link[fgsea]{fgseaMultilevel}}.
#'
#' @returns A \code{\link[data.table]{data.table}} with results for each
#'   contrast. See \code{\link[fgsea]{fgsea}} for more details.
#'
#' @md
#'
#' @importFrom fgsea fgseaMultilevel
#' @importFrom data.table rbindlist `:=` setDT
#' @importFrom tibble deframe
#'
#' @export fgsea2
#'
#' @author Tyler Sagendorf

fgsea2 <- function(pathways,
                   stats,
                   gene_column = "entrez_gene",
                   adjust.method = "BH",
                   adjust.globally = TRUE,
                   seed = NULL,
                   ...)
{
  # Check input
  adjust.method <- match.arg(adjust.method, choices = p.adjust.methods)

  if (!is.list(stats)) {
    stop("`stats` must be a named list.")
  }

  # List of pathways to test
  setDT(pathways)
  paths <- deframe(pathways[, list(gs_exact_source, get(gene_column))])

  # Run FGSEA on each contrast
  res <- lapply(names(stats), function(contrast_i) {
    message(contrast_i)
    set.seed(seed = seed)
    res_i <- fgseaMultilevel(pathways = paths,
                             stats = stats[[contrast_i]],
                             ...)
    res_i[["contrast"]] <- contrast_i

    return(res_i)
  })
  res <- rbindlist(res)
  setDT(res)


  pathways <- pathways[, list(gs_subcat, gs_exact_source, gs_description)]

  res <- merge(res, pathways, sort = FALSE,
               by.x = "pathway", by.y = "gs_exact_source")

  # p-value adjustment
  by <- "gs_subcat"
  if (!adjust.globally) { by <- c(by, "contrast") }
  res[, padj := p.adjust(pval, method = adjust.method), by = by]

  return(res)
}

