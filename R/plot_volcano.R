#' @title Create Volcano Plot
#'
#' @description Create a volcano plot.
#'
#' @param x a \code{matrix} or \code{data.frame} containing differential
#'   analysis results. Must include columns logFC, log10_pval, sign_logFC,
#'   contrast, and label.
#' @param pval_cutoff numeric; cutoff for p-values to be considered significant.
#'   Adds a dashed horizontal line to the plot.
#' @param colors character; length 3 vector of significantly negative,
#'   significantly positive, and non-significant (NS) logFC colors.
#'
#' @returns A \code{ggplot2} object.
#'
#' @md
#'
#' @import ggplot2
#' @importFrom scales extended_breaks
#' @importFrom data.table setDT `:=` `.N`
#'
#' @export plot_volcano

plot_volcano <- function(x,
                         pval_cutoff = 0.05,
                         colors = c("#3366ff", "darkred", "grey"))
{
  setDT(x)

  p <- ggplot(x, aes(x = logFC, y = log10_pval)) +
    geom_point(aes(color = sign_logFC),
               alpha = 0.5, size = 1, shape = 16) +
    scale_x_continuous(
      limits = range(x[["logFC"]], na.rm = T) * 1.1,
      expand = expansion(mult = 0),
      breaks = scales::extended_breaks(n = 5)
    ) +
    theme_pub() +
    theme(axis.line.y.right = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0,
                                    margin = margin(t=0, r=0, b=10, l=0)),
          legend.position = "none",
          panel.grid = element_line(linewidth = 0))

  # y-axis limits
  y_lims <- range_extend(x[["log10_pval"]], nearest = 1)
  # Use expand_limits instead
  # y_lims[1] <- 0
  # y_lims[2] <- max(3, y_lims[2])

  # Horizontal line and secondary axis for p-value cutoff
  p <- p +
    geom_hline(yintercept = -log10(pval_cutoff),
               linewidth = 0.3,
               lty = "dashed") +
    scale_y_continuous(
      breaks = scales::extended_breaks(n = 5),
      limits = y_lims,
      sec.axis = sec_axis(trans = ~ 10 ^ (-.),
                          breaks = pval_cutoff),
      expand = expansion(mult = c(0.01, 0.1))
    ) +
    expand_limits(y = c(0, 3)) +
    scale_color_manual(values = colors,
                       breaks = levels(x[["sign_logFC"]]))

  # add annotations
  p <- p +
    geom_label(data = unique(x[, list(contrast, label)]),
               aes(label = label, x = -Inf, y = Inf),
               size = 5 * 0.35, label.size = NA,
               label.padding = unit(4, "pt"),
               fill = alpha("white", 0.5),
               # label.r = unit(1.5, "pt"),
               hjust = 0.05, vjust = 0) +
    coord_cartesian(clip = "off")

  return(p)

  ## Label points (old code)
  # if (!missing(label)) {
  #   label_args <- list(mapping = aes(label = !!sym(label)),
  #                      na.rm = TRUE,
  #                      size = 5*0.352778, # convert points to mm
  #                      max.overlaps = Inf,
  #                      nudge_x = x[["nudge_x"]],
  #                      nudge_y = 0.1,
  #                      fill = alpha("white", 0.65),
  #                      force = 10, seed = 0,
  #                      color = "darkred",
  #                      min.segment.length = 0,
  #                      label.padding = 0.15)
  #   label_args <- modifyList(x = label_args, val = list(...),
  #                            keep.null = TRUE)
  #
  #   p <- p +
  #     do.call(what = geom_label_repel, args = label_args)
  # }
}

