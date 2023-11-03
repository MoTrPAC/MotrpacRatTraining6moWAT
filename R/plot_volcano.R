#' @title Create Volcano Plot
#'
#' @description Create a volcano plot from differential analysis results table.
#'
#' @param x a \code{matrix} or \code{data.frame} containing differential
#'   analysis results. Must include columns \code{logFC}, \code{pval}, contrast,
#'   and label.
#' @param pval_cutoff numeric; cutoff for adjusted p-values to be considered
#'   significant. Adds a dashed horizontal line to the plot.
#' @param colors character; length 3 vector of significantly negative,
#'   significantly positive, and non-significant (NS) logFC colors.
#'
#' @returns A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline sec_axis expansion
#'   scale_y_continuous scale_color_manual .pt
#' @importFrom scales extended_breaks
#'
#' @export plot_volcano
#'
#' @author Tyler Sagendorf

plot_volcano <- function(x,
                         pval_cutoff = 0.05,
                         colors = c("#3366ff", "darkred", "grey"))
{
  # Check input
  if (length(colors) != 3L) {
    stop("colors must be a vector of length 3.")
  }

  if (pval_cutoff < 0 | pval_cutoff > 1) {
    stop("pval_cutoff must be in [0, 1].")
  }

  p <- ggplot(x, aes(x = logFC, y = log10_pval)) +
    geom_point(aes(color = sign_logFC),
               alpha = 0.5, size = 1, shape = 16) +
    scale_x_continuous(
      limits = range(x[["logFC"]], na.rm = TRUE) * 1.1,
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

  # Horizontal line and secondary axis for p-value cutoff
  p <- p +
    geom_hline(yintercept = -log10(pval_cutoff),
               linewidth = 0.3, lty = "dashed") +
    scale_y_continuous(
      breaks = scales::extended_breaks(n = 5),
      limits = y_lims,
      sec.axis = sec_axis(trans = ~ 10 ^ (-.),
                          breaks = pval_cutoff),
      expand = expansion(mult = c(0.01, 0.1))
    ) +
    # Minimum y-axis limits
    expand_limits(y = c(0, 3)) +
    scale_color_manual(values = colors,
                       breaks = levels(x[["sign_logFC"]]))

  return(p)
}

# plot_volcano <- function(x,
#                          logFC = "logFC",
#                          pval = "adj.P.Val",
#                          pval_cutoff = 0.05,
#                          contrast = "contrast",
#                          label = "feature",
#                          nearest = 1,
#                          colors = c("#3366ff", "darkred", "grey"))
# {
#   # Check input
#   for (col_i in c(logFC, pval, contrast, label)) {
#     if (!col_i %in% colnames(x)) {
#       stop(sprintf("%s is not the name of a column in x.", col_i))
#     }
#   }
#
#   x <- x %>%
#     dplyr::select(label = !!sym(label),
#                   contrast = !!sym(contrast),
#                   logFC = !!sym(logFC),
#                   pval = !!sym(pval)) %>%
#     mutate(log10_pval = -log10(pval),
#            sign_logFC = sign(logFC) * (pval < pval_cutoff),
#            sign_logFC = factor(sign_logFC,
#                                levels = c(-1, 1, 0),
#                                labels = c("up", "down", "NS")))
#
#   # Counts of differential features by contrast
#   contr_summary <- x %>%
#     group_by(contrast, sign_logFC) %>%
#     summarise(n = n()) %>%
#     mutate(contr_label = paste(n, sign_logFC)) %>%
#     group_by(contrast) %>%
#     summarise(contr_label = paste(contr_label, collapse = ", "))
#
#   x <- left_join(x, contr_summary, by = "contrast")
#
#   p <- ggplot(x, aes(x = logFC, y = log10_pval)) +
#     geom_point(aes(color = sign_logFC),
#                alpha = 0.5, size = 1, shape = 16)
#
#   # Horizontal line and secondary axis for p-value cutoff
#   y_lims <- range_extend(x$log10_pval, nearest = 1)
#
#   p <- p +
#     geom_hline(yintercept = -log10(pval_cutoff),
#                linewidth = 0.3,
#                lty = "dashed") +
#     scale_y_continuous(breaks = extended_breaks(n = 5),
#                        limits = y_lims,
#                        sec.axis = sec_axis(trans = ~ 10 ^ (-.x),
#                                            breaks = pval_cutoff),
#                        expand = expansion(mult = c(0.01, 0.1))) +
#     expand_limits(y = c(0, 3))
#
#   # Modify appearance
#   p <- p +
#     scale_x_continuous(limits = range_extend(x$logFC, nearest = nearest),
#                        expand = expansion(mult = 0),
#                        breaks = extended_breaks(n = 5)) +
#     scale_color_manual(values = colors,
#                        breaks = levels(x$sign_logFC)) +
#     theme_pub() +
#     theme(axis.line.y.right = element_blank(),
#           strip.background = element_blank(),
#           strip.text = element_text(hjust = 0,
#                                     margin = margin(t=0, r=0, b=10, l=0)),
#           legend.position = "none",
#           panel.grid = element_line(linewidth = 0))
#
#   return(p)
# }

