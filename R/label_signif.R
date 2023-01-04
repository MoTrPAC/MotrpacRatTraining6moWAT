#' @title Label p-values according to significance
#'
#' @description Use asterisks to indicate level of statistical significance.
#'
#' @param pvals numeric; vector of p-values.
#'
#' @return character vector of labels: p > 0.05 = "", 0.01 <= p < 0.05 = "\*",
#'   0.001 <= p < 0.01 = "\*\*", p < 0.001 = "\*\*\*"
#'
#' @export label_signif

label_signif <- function(pvals)
{
  cut(pvals, breaks = c(0, 0.001, 0.01, 0.05, 1),
      labels = c("***", "**", "*", ""),
      include.lowest = TRUE, right = FALSE)
}

