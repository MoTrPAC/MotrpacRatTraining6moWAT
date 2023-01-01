#' @title Visualize WGCNA module eigenfeatures
#'
#' @description Create scatterplots of each module eigenfeature (principal
#'   eigenvectors).
#'
#' @param x output of \code{\link{run_WGCNA}}.
#'
#' @import ggplot2
#' @importFrom scales breaks_extended
#' @importFrom ggbeeswarm position_beeswarm
#' @importFrom patchwork wrap_plots
#' @importFrom data.table `:=` `.N` `.SD` as.data.table
#'
#' @export plot_eigenfeature

plot_eigenfeature <- function(x)
{
  # Module sizes
  mod_size <- as.data.table(x$modules)
  mod_size <- mod_size[, list(n = .N), by = moduleID]

  # Module MEs
  x <- as.data.table(x$MEs)
  x <- subset(x, subset = moduleNum != 0) # remove grey module
  x <- merge(x, mod_size, by = "moduleID", all.x = TRUE, all.y = FALSE)
  # x <- subset(x, subset = moduleColor != "grey")
  x[, moduleID := droplevels(moduleID)]
  # x[, `:=` (sex = factor(sex, levels = sort(unique(sex)),
  #                        labels = c("F", "M")),
  #           timepoint = factor(timepoint,
  #                              levels = c("SED", paste0(2^(0:3), "W"))))]

  ## Combine plots
  plotlist <- lapply(levels(x$moduleID), function(mod_i) {
    xi <- subset(x, moduleID == mod_i)

    mod_i_size <- xi[["n"]][1] # module size

    set.seed(0)
    ggplot(xi, aes(x = timepoint, y = ME)) +
      stat_summary(fun.data = "mean_sdl",
                   fun.args = list(mult = 1),
                   mapping = aes(color = sex),
                   geom = "crossbar", width = 0.7,
                   fatten = 1, size = 0.5) +
      geom_point(size = 0.4, color = "black", shape = 16,
                 position = position_beeswarm(cex = 3, groupOnX = TRUE,
                                              dodge.width = 0.34)) +
      facet_wrap(~ sex, nrow = 1) +
      labs(title = sprintf("%s (n = %d)", mod_i, mod_i_size),
           y = NULL, x = NULL) +
      scale_color_manual(name = "Sex",
                         values = c("#ff6eff", "#5555ff")) +
      scale_y_continuous(breaks = breaks_extended()) +
      # theme_bw() +
      theme_pub() +
      theme(# axis.title = element_text(size = 6, color = "black"),
            # axis.text = element_text(size = 5, color = "black"),
            # axis.line = element_line(color = "black", size = 0.3),
            # axis.ticks = element_line(color = "black", size = 0.3),
            # panel.grid.major.y = element_line(size = 0.3),
            panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(),
            # plot.title = element_text(size = 6.5),
            axis.text.x = element_text(angle = 90,
                                       hjust = 1, vjust = 0.5,
                                       margin = margin(t = 0)),
            panel.spacing = unit(0, "points"),
            strip.text = element_blank(),
            # panel.border = element_blank(),
            # legend.title = element_text(color = "black", size = 6),
            # legend.text = element_text(color = "black", size = 5),
            plot.margin = unit(rep(3, 4), units = "pt")
      )
  })

  p <- wrap_plots(plotlist, ncol = 4, guides = "collect")

  return(p)
}

