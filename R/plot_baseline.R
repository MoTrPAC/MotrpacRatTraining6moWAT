#' @title Plot baseline (pre-training) measures
#'
#' @description Plot baseline (pre-training) measures.
#'
#' @param x a \code{data.frame} with columns \code{"sex"}, \code{"age"},
#'   \code{"group"}, and a column specified by \code{response}.
#' @param response character; the name of a column in \code{x} to plot.
#' @param conf a \code{data.frame} with confidence interval data.
#' @param stats optional \code{data.frame} with results of statistical analyses
#'   filtered to the response of interest.
#' @param bracket.nudge.y numeric;
#'
#' @import ggplot2
#' @importFrom dplyr %>% rename filter mutate pull
#' @importFrom purrr map2_dbl
#' @importFrom rlang !! sym
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom ggbeeswarm position_beeswarm
#'
#' @export plot_baseline

plot_baseline <- function(x,
                          response,
                          conf,
                          stats,
                          bracket.nudge.y = 0)
{
  x <- x %>%
    dplyr::rename(response = !!sym(response)) %>%
    filter(!is.na(response)) %>%
    droplevels.data.frame()

  # Base plot
  p <- ggplot(data = x, aes(x = group, y = response))

  # Add
  if (missing(conf)) {
    p <- p +
      stat_summary(fun.data = ~ exp(mean_cl_normal(log(.x))),
                   mapping = aes(color = sex),
                   geom = "crossbar", fatten = 1, linewidth = 0.4)

  } else {
    conf <- filter(conf, group %in% x$group)

    if (any(grepl("^asymp", colnames(conf)))) {
      conf <- conf %>%
        mutate(lower.CL = ifelse(is.na(lower.CL), asymp.LCL, lower.CL),
               upper.CL = ifelse(is.na(upper.CL), asymp.UCL, upper.CL))
    }
    p <- p +
      geom_crossbar(aes(x = group, y = response_mean,
                        ymin = lower.CL, ymax = upper.CL,
                        color = sex),
                    data = conf, fatten = 1, linewidth = 0.4)
  }

  # Workaround to prevent warning about existing coordinate system
  cart <- coord_cartesian(clip = "off")
  cart$default <- TRUE

  p <- p +
    geom_point(color = "black", na.rm = TRUE, shape = 16,
               position = position_beeswarm(cex = 3.5, dodge.width = 0.7),
               size = 0.6) +
    scale_color_manual(values = c("#ff63ff", "#5555ff"),
                       breaks = c("Female", "Male")) +
    labs(x = NULL) +
    facet_grid(~ sex) +
    cart +
    theme_bw() +
    theme(text = element_text(size = 6.5, color = "black"),
          line = element_line(linewidth = 0.3, color = "black"),
          axis.ticks = element_line(linewidth = 0.3, color = "black"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 5,
                                   color = "black"),
          axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1,
                                     vjust = 0.5),
          axis.title = element_text(size = 6.5, margin = margin(),
                                    color = "black"),
          axis.line = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(size = 6.5, color = "black"),
          panel.spacing = unit(-1, "pt"),
          plot.title = element_text(size = 7, color = "black"),
          plot.subtitle = element_text(size = 6, color = "black"),
          legend.position = "none",
          strip.placement = "outside"
    )

  if (!missing(stats)) {
    # Add asterisks and brackets for significant comparisons
    stats <- stats %>%
      filter(signif != "") %>%
      mutate(
        groups = strsplit(as.character(contrast), split = " / "),
        y.position = map2_dbl(.x = groups, .y = sex, .f = function(.x, .y) {
          x %>%
            filter(group %in% levels(group)[1:which(levels(group) == .x[1])],
                   sex == .y) %>%
            pull(response) %>%
            max(na.rm = TRUE)
        }),
        group1 = sub(" .*", "", contrast),
        group2 = sub(".* ", "", contrast))


    if (nrow(stats) > 0) {
      p <- p +
        stat_pvalue_manual(stats,
                           bracket.nudge.y = bracket.nudge.y,
                           label.size = 3,
                           step.group.by = "sex",
                           vjust = 0.5,
                           step.increase = 0.075,
                           label = "signif")
    }
  }

  return(p)
}

