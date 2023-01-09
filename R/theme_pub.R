#' @title Custom ggplot2 theme
#'
#' @description Theme used extensively by the plotting functions of this
#'   package, as well as to generate final plots for the associated publication.
#'   Built to adhere to the Nature publishing guidelines.
#'
#' @import ggplot2
#'
#' @export theme_pub

theme_pub <- function() {
  theme_bw() +
    theme(
      text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 6, color = "black"),
      axis.text = element_text(size = 5, color = "black"),
      axis.text.y.right = element_text(size = 5, color = "black"),
      axis.line = element_line(size = 0.3, color = "black"),
      axis.ticks = element_line(size = 0.3, color = "black"),

      legend.title = element_text(size = 5.5, color = "black"),
      legend.text = element_text(size = 5, color = "black"),

      panel.border = element_blank(),
      panel.grid.major = element_line(size = 0.3),
      panel.spacing = unit(5, "pt"),

      plot.title = element_text(size = 6.5, color = "black"),
      plot.subtitle = element_text(size = 6, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),

      strip.text = element_text(size = 6.5, color = "black")
    )
}

