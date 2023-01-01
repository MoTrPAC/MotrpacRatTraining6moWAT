#' @title Visualize clinical analyte data
#'
#' @description Generate a scatterplot of clinical analyte values.
#'
#' @param x a `data.frame` with required columns value, sex, and timepoint.
#' @param analyte character; the name of a column of `x`.
#'
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom ggnewscale new_scale_color
#'
#' @export plot_analyte


## Function to plot a specific analyte ----
plot_analyte <- function(x, analyte = "")
{
  ggplot(x, aes(x = timepoint, y = !!sym(analyte), color = sex)) +
    stat_summary(fun.data = ~ exp(mean_sdl(log(.x), mult = 1)),
                 # fun.args = list(mult = 1),
                 show.legend = FALSE,
                 geom = "crossbar", width = 0.8,
                 na.rm = TRUE, fatten = 1, size = 0.4) +
    scale_color_manual(values = c("#ff6eff", "#5555ff"),
                       breaks = c("Female", "Male")) +
    guides(color = NULL) +
    # new_scale_color() +
    geom_quasirandom(color = "black", width = 0.35, size = 0.4, na.rm = T) +
    facet_wrap(~ sex, scales = "free_x") +
    # scale_fill_manual(values = c("#ff6eff", "#5555ff")) +
    # scale_color_manual(name = NULL,
    #                    values = c("black", "red"),
    #                    breaks = c(FALSE, TRUE),
    #                    labels = c("Non-outlier", "Outlier")) +
    scale_x_discrete(name = NULL) +
    # theme_bw() +
    theme_pub() +
    theme(#text = element_text(size = 6.5, color = "black"),
          # line = element_line(size = 0.3, color = "black"),
          # axis.ticks = element_line(size = 0.3, color = "black"),
          panel.grid = element_blank(),
          # panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.text = element_text(size = 5,
          #                          color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          # axis.title = element_text(size = 6.5,
          #                           color = "black"),
          axis.title.y = element_text(margin = margin(r = 2.5,
                                                      unit = "pt")),
          # axis.line = element_line(size = 0.3, color = "black"),
          strip.background = element_blank(),
          # strip.text = element_text(size = 6.5),
          panel.spacing = unit(-1, "pt"),
          # plot.title = element_text(size = 7, color = "black"),
          # plot.subtitle = element_text(size = 6, color = "black"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          # legend.text = element_text(size = 5, color = "black"),
          legend.margin = margin(t = -3, b = -3, unit = "pt"))
}

