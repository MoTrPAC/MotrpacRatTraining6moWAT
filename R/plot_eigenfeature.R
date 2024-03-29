#' @title Plot WGCNA module eigenfeatures
#'
#' @description Create scatterplots of each module eigenfeature (principal
#'   eigenvector).
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
#'
#' @author Tyler Sagendorf

plot_eigenfeature <- function(x)
{
  # Module sizes
  mod_size <- as.data.table(x$modules)
  mod_size <- mod_size[, list(n = .N), by = moduleID]

  # Module MEs
  x <- as.data.table(x$MEs)
  x <- subset(x, subset = moduleNum != 0) # remove grey module
  x <- merge(x, mod_size, by = "moduleID", all.x = TRUE, all.y = FALSE)
  x[, moduleID := droplevels(moduleID)]

  # Combine plots
  plotlist <- lapply(levels(x$moduleID), function(mod_i) {
    xi <- subset(x, moduleID == mod_i)

    mod_i_size <- xi[["n"]][1] # module size

    set.seed(0)
    ggplot(xi, aes(x = timepoint, y = ME)) +
      stat_summary(fun.data = ~ mean_sdl(.x),
                   mapping = aes(color = sex),
                   geom = "crossbar", width = 0.7,
                   fatten = 1, linewidth = 0.5) +
      geom_point(size = 0.4, color = "black", shape = 16,
                 position = position_beeswarm(cex = 3, dodge.width = 0.34)) +
      facet_wrap(~ sex, nrow = 1) +
      labs(title = sprintf("%s (n = %d)", mod_i, mod_i_size),
           x = NULL, y = NULL) +
      scale_color_manual(name = "Sex:",
                         values = c("#ff6eff", "#5555ff")) +
      scale_y_continuous(breaks = breaks_extended()) +
      theme_pub() +
      theme(panel.grid.major.y = element_line(linewidth = 0),
            panel.grid.minor.y = element_line(linewidth = 0),
            panel.grid.major.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 90,
                                       hjust = 1, vjust = 0.5,
                                       margin = margin(t = 0)),
            panel.spacing = unit(0, "points"),
            strip.text = element_blank(),
            plot.margin = unit(rep(3, 4), units = "pt")
      )
  })

  p <- wrap_plots(plotlist, ncol = 4, guides = "collect")

  return(p)
}

