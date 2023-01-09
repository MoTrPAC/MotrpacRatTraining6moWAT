#' @title MDS plots with ggplot2
#'
#' @description A \code{ggplot2} version of \code{\link[limma]{plotMDS}}
#'   tailored to the MoTrPAC PASS1B data. Samples are labeled by the number of
#'   weeks of exercise training (0, 1, 2, 4, 8) and colored according to sex.
#'
#' @param object object of class \code{link[MSnbase:MSnSet-class]{MSnSet}}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @import ggplot2
#' @importFrom limma plotMDS
#' @importFrom MSnbase exprs
#'
#' @export ggplotMDS

ggplotMDS <- function(object) {
  mds.obj <- plotMDS(exprs(object), plot = F)
  var_expl <- mds.obj[["var.explained"]]
  var_expl <- sprintf("Leading logFC dim %g (%g%%)",
                      seq_along(var_expl), round(100 * var_expl))[1:2]

  point_label <- ifelse(object$timepoint == "SED", "0",
                        sub("W", "", object$timepoint))

  p <- ggplot(mapping = aes(x = mds.obj$x, y = mds.obj$y)) +
    labs(x = var_expl[1], y = var_expl[2]) +
    geom_text(aes(label = point_label, color = object$sex),
              size = 2.46) +
    scale_color_manual(name = "Sex",
                       values = c("#ff6eff", "#3366ff"),
                       breaks = c("Female", "Male")) +
    labs(subtitle = "Samples labeled by weeks of ExT") +
    theme_pub() +
    theme(axis.line.y.right = element_blank(),
          strip.text = element_text(hjust = 0,
                                    margin = margin(t=0, r=0, b=5, l=0)),
          legend.key.size = unit(6, "pt"),
          legend.margin = margin(r = 0, l = 0),
          panel.grid.minor = element_line(linewidth = 0))

  return(p)
}

