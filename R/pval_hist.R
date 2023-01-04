#' @title Create p-value Histograms
#'
#' @description Create p-value histograms with
#'   \code{\link[ggplot2]{ggplot2-package}}.
#'
#' @inheritParams plot_volcano
#' @param pval character; the name of a column in \code{x} containing p-values
#'   or adjusted p-values.
#' @param ylims numeric vector of length 2. Limits of y-axis.
#'
#' @returns A \code{ggplot2} object.
#'
#' @import ggplot2
#' @importFrom scales extended_breaks
#'
#' @export pval_hist

pval_hist <- function(x,
                      pval = "adj.P.Val",
                      ylims = c(0, NA)) {
  ggplot(x) +
    geom_histogram(aes(x = !!sym(pval)),
                   breaks = seq(0, 1, 0.05),
                   color = "black", fill = "lightgrey") +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       expand = expansion(5e-3)) +
    scale_y_continuous(name = "Count",
                       limits = ylims,
                       breaks = scales::extended_breaks(n = 5),
                       expand = expansion(mult = c(0.001, 0.05))) +
    theme_pub()
    # theme_bw() +
    # theme(text = element_text(size = 6.5*scale, color = "black"),
    #       line = element_line(size = 0.3*scale, color = "black"),
    #       panel.border = element_rect(size = 0.4*scale,
    #                                   fill = NULL,
    #                                   color = "black"),
    #       axis.line.y.right = element_blank(),
    #       axis.text = element_text(size = 5*scale),
    #       strip.background = element_blank(),
    #       strip.text = element_text(size = 6.5*scale,
    #                                 hjust = 0),
    #       panel.spacing = unit(0.08*scale, "in"),
    #       plot.title = element_text(size = 7*scale),
    #       panel.spacing.x = unit(0.3, "in"),
    #       plot.margin = unit(c(5, 10, 5, 5), "pt"))
}

