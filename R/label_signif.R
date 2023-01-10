#' @title Label p-values according to significance
#'
#' @description Use asterisks to indicate level of statistical significance.
#'
#' @param pvals numeric; vector of p-values.
#'
#' @return character vector of labels: p \eqn{>} 0.05 = "", 0.01 \eqn{\leq} p
#'   \eqn{<} 0.05 = "*", 0.001 \eqn{\leq} p \eqn{<} 0.01 = "**", p \eqn{<} 0.001
#'   = "***"
#'
#' @noMd
#'
#' @export label_signif

label_signif <- function(pvals)
{
  cut(pvals, breaks = c(0, 0.001, 0.01, 0.05, 1),
      labels = c("***", "**", "*", ""),
      include.lowest = TRUE, right = FALSE)
}

