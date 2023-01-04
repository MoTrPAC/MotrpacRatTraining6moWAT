#' @title Visualize clinical analyte data
#'
#' @description Generate a scatterplot of clinical analyte values.
#'
#' @param x a \code{data.frame} with required columns value, sex, and timepoint.
#' @param analyte character; the name of a column of \code{x}.
#'
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#'
#' @export plot_analyte


## Function to plot a specific analyte ----
plot_analyte <- function(x, analyte = "")
{
  ggplot(x, aes(x = timepoint, y = !!sym(analyte), color = sex)) +
    stat_summary(fun.data = ~ exp(mean_se(log(.x))),
                 show.legend = FALSE,
                 geom = "crossbar", width = 0.8,
                 na.rm = TRUE, fatten = 1, size = 0.4) +
    scale_color_manual(values = c("#ff6eff", "#5555ff"),
                       breaks = c("Female", "Male")) +
    guides(color = NULL) +
    geom_quasirandom(color = "black", width = 0.35, size = 0.4, na.rm = T) +
    facet_wrap(~ sex, scales = "free_x") +
    scale_x_discrete(name = NULL) +
    theme_pub() +
    theme(panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.y = element_text(margin = margin(r = 2.5,
                                                      unit = "pt")),
          strip.background = element_blank(),
          panel.spacing = unit(-1, "pt"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.margin = margin(t = -3, b = -3, unit = "pt"))
}

