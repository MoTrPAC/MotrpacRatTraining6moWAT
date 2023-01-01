#' @title Plot module ORA results
#'
#' @description Dotplot of the top over-represented terms by WGCNA module.
#'
#' @param x output of \code{\link[fgsea]{fora}} or \code{\link{fora2}}.
#' @param n_terms integer; maximum number of top over-represented terms to
#'   display from each module.
#' @param mods integer; the module numbers to include in the plot. By default,
#'   terms from modules 1-5 (the 5 largest modules) are shown.
#' @param gs_subcat character; the gene set subcategory of terms to plot.
#'   Default is "GO:MF".
#' @param rel_heights numeric; vector of length 2 specifying the relative
#'   heights of the plot and color legend, respectively. Default is c(0.8, 0.2).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @md
#'
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom cowplot get_legend plot_grid
#' @importFrom scales breaks_extended
#' @importFrom data.table setDT setorder `:=`
#'
#' @export plot_ORA

plot_ORA <- function(x,
                     n_terms = 5, # number of top terms per module
                     mods = 1:5, # module numbers to plot
                     gs_subcat = "GO:MF",
                     rel_heights = c(0.8, 0.2)
)
{
  # TODO add input checks

  mods <- mods[mods %in% 1:nlevels(x$module)]

  setDT(x)
  x <- subset(x[, , with=FALSE], subset = (gs_subcat %in% gs_subcat) &
                (padj < 0.05) & (module %in% levels(module)[mods]))
  setorder(x, c("module", "padj", "pval"))
  x <- x[1:n_terms, , by = module]
  x[, row_labels := ifelse(nchar(row_labels) > 35 + nchar(pathway) + 5,
                           sprintf("%s...(%s)",
                                   substr(gs_description, 1, 30),
                                   gs_description),
                           gs_description)]
  x[, row_labels := factor(row_labels, levels = rev(unique(row_labels)))]


  p1 <- ggplot(x) +
    geom_point(aes(x = module, y = row_labels,
                   size = -log10(padj),
                   color = overlapRatio)) +
    scale_size_area(name = latex2exp::TeX("$-log_{10}$(scaled p-value)"),
                    breaks = scales::breaks_extended(),
                    max_size = 6) +
    scale_color_viridis_c(name = "Overlap Ratio",
                          direction = -1,
                          limits = c(0, 1),
                          breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = NULL) +
    guides(
      size = guide_legend(title.position = "top",
                          order = 1),
      color = guide_colorbar(
        title.position = "top",
        label.theme = element_text(angle = 45,
                                   hjust = 1, vjust = 1),
        title.theme = element_text(margin = margin(r = 5)),
        barheight = unit(5, "pt"),
        barwidth = unit(60, "pt"),
        frame.colour = "black",
        ticks.colour = "black",
        order = 2)) +
    scale_y_discrete(position = "right") +
    # theme_bw(base_size = 6) +
    theme_pub() +
    theme(#text = element_text(size = 6, color = "black"),
      line = element_line(size = 0.3, color = "black"),
      axis.ticks = element_blank(),
      # axis.text = element_text(size = 6, color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      # axis.title = element_text(size = 6.5, color = "black"),
      # legend.title = element_text(size = 6, color = "black"),
      # legend.text = element_text(size = 5, color = "black"),
      legend.position = "bottom",
      legend.box.just = "left",
      legend.direction = "horizontal",
      legend.margin = margin(t = 4, r = 5, b = 4, l = 8),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 10))

  p2 <- p1 + theme(legend.position = "none")
  le1 <- cowplot::get_legend(p1)
  p <- cowplot::plot_grid(p2, le1, nrow = 2, rel_heights = rel_heights)

  return(p)
}

