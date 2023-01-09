#' @title Plot module ORA results
#'
#' @description Dotplot of the top over-represented terms by WGCNA module.
#'
#' @param x output of \code{\link[fgsea]{fora}} or \code{\link{fora2}}.
#' @param n_terms integer; maximum number of top over-represented terms to
#'   display from each module.
#' @param mods integer; the module numbers to include in the plot. By default,
#'   terms from modules 1-5 (the 5 largest modules) are shown.
#' @param subset character; the gene set subcategory of terms to plot.
#'   Default is "GO:MF".
#' @param rel_heights numeric; vector of length 2 specifying the relative
#'   heights of the plot and color legend, respectively. Default is c(0.8, 0.2).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @importFrom latex2exp TeX
#' @importFrom cowplot get_legend plot_grid
#' @importFrom scales breaks_extended
#' @importFrom data.table setDT setorder `:=` `.SD`
#'
#' @export plot_ORA

plot_ORA <- function(x,
                     n_terms = 5, # number of top terms per module
                     mods = 1:5, # module numbers to plot
                     subset = "GO:MF",
                     rel_heights = c(0.8, 0.2))
{
  # TODO add input checks

  mods <- mods[mods %in% 1:nlevels(x$module)]

  setDT(x)
  x <- subset(x,
              subset = (gs_subcat %in% subset) &
                (padj < 0.05) &
                (module %in% levels(module)[mods]))
  # x <- subset(x[, , with = FALSE],
  #             subset = (gs_subcat %in% gs_subcat) &
  #               (padj < 0.05) & (module %in% levels(module)[mods]))
  setorder(x, module, padj, pval)
  x <- x[, head(.SD, n_terms), by = module]
  x[, row_labels := ifelse(nchar(gs_description) > 35 + nchar(pathway) + 5,
                           sprintf("%s...(%s)",
                                   substr(gs_description, 1, 30),
                                   pathway),
                           gs_description)]
  x[, row_labels := factor(row_labels, levels = rev(unique(row_labels)))]

  p1 <- ggplot(x) +
    geom_point(aes(x = module, y = row_labels,
                   size = -log10(padj),
                   color = overlapRatio)) +
    scale_size_area(name = latex2exp::TeX("$-log_{10}$(scaled p-value)"),
                    breaks = scales::breaks_extended(),
                    max_size = 3) +
    scale_color_viridis_c(name = "Overlap Ratio",
                          direction = -1,
                          limits = c(0, 1),
                          breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = NULL) +
    guides(
      size = guide_legend(title.position = "top",
                          title.vjust = -0.5,
                          keyheight = unit(5, "pt"),
                          keywidth = unit(5, "pt"),
                          title.theme = element_text(size = 6,
                                                     margin = margin(b = 0)),
                          order = 1),
      color = guide_colorbar(
        title.position = "top",
        label.theme = element_text(size = 5, angle = 45,
                                   hjust = 1, vjust = 1),
        title.theme = element_text(size = 6, margin = margin(r = 5)),
        barheight = unit(5, "pt"),
        barwidth = unit(60, "pt"),
        frame.colour = "black",
        ticks.colour = "black",
        order = 2)) +
    scale_y_discrete(position = "right") +
    theme_pub() +
    theme(line = element_line(linewidth = 0.3, color = "black"),
          axis.ticks = element_line(linetype = 0),
          axis.text.x = element_text(size = 6,
                                     angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y.right = element_text(size = 6),
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

