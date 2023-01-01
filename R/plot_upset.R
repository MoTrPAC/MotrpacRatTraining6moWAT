#' @title Create Annotated UpSet Plot
#'
#' @description Create an UpSet plot annotated with intersection and set sizes.
#'   Essentially a modified version of \code{\link[ComplexHeatmap]{UpSet}}.
#'
#' @inheritParams enrichmat
#' @param x list; named list of elements in each set. Passed to
#'   \code{\link[ComplexHeatmap]{make_comb_mat}}.
#' @param top_n_comb numeric; number of largest intersections to show. Default
#'   is 10.
#' @param comb_title character; title of top barplot annotation. Default is
#'   "Intersection Size".
#' @param set_title character; title of right barplot annotation. Default is
#'   "Set Size".
#' @param scale numeric; scaling factor for all plot elements. Improves
#'   appearance. Default is 2.
#' @param barplot_args list; additional arguments passed to
#'   \code{\link[ComplexHeatmap]{anno_barplot}}.
#' @param filename character; the filename used to save the UpSet plot. If
#'   provided, the plot will not be displayed.
#' @param height numeric; the save height.
#' @param width numeric; the save width.
#' @param units character; the units of `height` and `width`. Default is "in"
#'   (inches).
#' @param upset_args list; additional arguments passed to
#'   \code{\link[ComplexHeatmap]{UpSet}}.
#'
#' @md
#'
#' @details See \code{\link[ComplexHeatmap]{UpSet}} for details and sections
#'   mentioned in Arguments.
#'
#' @import ComplexHeatmap
#' @importFrom grid gpar unit
#' @importFrom utils modifyList
#' @importFrom grDevices dev.off
#'
#' @export plot_upset
#'
#' @references Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister,
#'   H. (2014). UpSet: Visualization of Intersecting Sets. *IEEE transactions on
#'   visualization and computer graphics, 20*(12), 1983--1992.
#'   \url{https://doi.org/10.1109/TVCG.2014.2346248}

plot_upset <- function(x,
                       top_n_comb = 10,
                       comb_title = "Intersection\nSize",
                       set_title = "Set Size",
                       scale = 2,
                       barplot_args = list(),
                       filename = character(0),
                       height = 2*scale,
                       width = 3.4*scale,
                       units = "in",
                       save_args = list(),
                       upset_args = list())
{
  m <- make_comb_mat(x)
  cs <- comb_size(m)
  ss <- set_size(m)
  m <- m[, order(-cs)[1:min(length(cs), top_n_comb)]]
  cs <- comb_size(m)
  row_order <- seq_along(ss)
  column_order <- order(-cs)

  ## Base UpSet plot ----
  upset_args <- modifyList(
    x = list(
      m = m,
      row_names_gp = gpar(fontsize = 7*scale),
      width = length(cs)*unit(12, "pt")*scale,
      height = length(ss)*unit(12, "pt")*scale,
      lwd = 1.8*scale,
      pt_size = unit(7*scale, "pt"),
      comb_order = column_order,
      set_order = row_order,
      row_labels = set_name(m)
    ),
    val = upset_args, keep.null = TRUE)

  up <- do.call(what = UpSet,
                args = upset_args)

  ## Modify barplot annotations ----
  barplot_args <- modifyList(
    x = list(
      height = 6*unit(12, "pt")*scale,
      width = 6*unit(12, "pt")*scale,
      gp = gpar(fill = "black", color = "black"),
      add_numbers = TRUE,
      border = FALSE,
      numbers_gp = gpar(fontsize = 6*scale),
      axis = FALSE
    ),
    val = barplot_args, keep.null = TRUE)

  barplot_args[["row_names_max_width"]] <- max_text_width(
    text = upset_args[["row_labels"]],
    gp = upset_args[["row_names_gp"]]
  )



  if (!("top_annotation" %in% names(upset_args))) {
    up@top_annotation <- HeatmapAnnotation(
      comb_size = do.call(what = anno_barplot,
                          args = c(list(x = cs), barplot_args)),
      annotation_name_gp = gpar(fontsize = 7*scale),
      annotation_name_rot = 0,
      annotation_label = comb_title,
      annotation_name_side = "left",
      height = unit(12, "pt")*6*scale
    )
  }

  if (!("right_annotation" %in% names(upset_args))) {
    up@right_annotation <- HeatmapAnnotation(
      set_size = do.call(what = anno_barplot,
                         args = c(list(x = ss), barplot_args)),
      annotation_name_gp = gpar(fontsize = 7*scale),
      annotation_name_rot = 0,
      annotation_label = set_title,
      which = "row",
      width = unit(12, "pt")*4*scale
    )
  }

  ## Save plot ----
  if (!identical(filename, character(0))) {
    save_heatmap(filename = filename,
                 height = height,
                 width = width,
                 units = units,
                 save_args = save_args)
    on.exit(expr = dev.off())
  }

  draw(up)
}

