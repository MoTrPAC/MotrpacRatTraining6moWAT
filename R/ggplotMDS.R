#' @title MDS plots with ggplot2
#'
#' @description A `ggplot2` version of \code{\link[limma]{plotMDS}} tailored to
#'   the MoTrPAC PASS1B data. Samples are labeled by the number of weeks of
#'   exercise training (0, 1, 2, 4, 8) and colored according to sex.
#'
#' @param x an \code{link[MSnbase]{MSnSet}} object.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @importFrom limma plotMDS
#' @importFrom MSnbase exprs
#'
#' @export ggplotMDS

ggplotMDS <- function(x) {
  mds.obj <- plotMDS(exprs(x), plot = F)
  var_expl <- mds.obj[["var.explained"]]
  var_expl <- sprintf("Leading logFC dim %g (%g%%)",
                      seq_along(var_expl), round(100 * var_expl))[1:2]

  point_label <- ifelse(x$timepoint == "SED", "0",
                        sub("W", "", x$timepoint))

  p <- ggplot(mapping = aes(x = mds.obj$x, y = mds.obj$y)) +
    labs(x = var_expl[1], y = var_expl[2]) +
    geom_text(aes(label = point_label, color = x$sex),
              size = 2.46*scale) +
    scale_color_manual(name = "Sex",
                       values = c("#ff6eff", "#3366ff"),
                       breaks = c("F", "M"),
                       labels = c("Female", "Male")) +
    labs(subtitle = "Samples labeled by weeks of ExT") +
    theme_pub() +
    # theme_bw() +
    theme(#text = element_text(size = 6.5*scale, color = "black"),
          # line = element_line(size = 0.3*scale, color = "black"),
          axis.line.y.right = element_blank(),
          # axis.text = element_text(size = 5*scale, color = "black"),
          # axis.title = element_text(size = 6.5*scale,
          #                           color = "black"),
          # strip.background = element_blank(),
          # panel.border = element_rect(fill = NULL, color = "black",
          #                             size = 0.3*scale),
          strip.text = element_text(#size = 6.5*scale,
                                    hjust = 0,
                                    margin = margin(t=0, r=0, b=5, l=0))#,
          # panel.spacing = unit(5, "pt"),
          # plot.title = element_text(size = 7, color = "black")
          )

  return(p)
}

