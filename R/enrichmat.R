#' @title Enrichment Heatmap
#'
#' @description Creates a heatmap summarizing FGSEA results.
#'
#' @param x object of class \code{data.frame}; output of
#'   \code{\link[fgsea]{fgsea}} or \code{\link{fgsea2}}. Must contain a
#'   "contrast" column.
#' @param n_top integer; number of pathways to display. Defaults to the top 15
#'   pathways that are most significantly enriched across all contrasts.
#' @param top_pathways character; vector of specific pathways to display. If
#'   \code{NULL} (default), the \code{n_top} pathways will be displayed instead.
#' @param rownames_column character; the name of a column in \code{x} containing
#'   unique identifiers that will be used as the rownames in the heatmap.
#'   Default is "pathway".
#' @param NES_column similar to \code{rownames_column}. The name of a column
#'   containing the normalized enrichment scores (NES) that determines the
#'   heatmap body colors. Values between -1 and +1 (noise) will appear white.
#' @param padj_column similar to \code{rownames_column}. The name of a column
#'   containing the adjusted p-values that determine the area of each circle in
#'   the heatmap.
#' @param padj_legend_title character; title of the background fill legend.
#'   Defaults to \code{padj_column}.
#' @param padj_cutoff numeric; cutoff for terms to be statistically significant.
#'   If \code{plot_sig_only=TRUE}, only those pathways with at least one
#'   \code{padj_column} value less than this threshold may appear in the
#'   heatmap. Default is 0.05.
#' @param plot_sig_only logical; whether to plot only those \code{n_top} terms
#'   that have at least one \code{padj_column} value less than
#'   \code{padj_cutoff}.
#' @param padj_fill character; the background color used for values in
#'   \code{padj_column} that are less than \code{padj_cutoff}. Default is
#'   "lightgrey".
#' @param colors character; vector of length 2 specifying the colors for the
#'   largest negative and largest positive values, respectively. Defaults to
#'   blue (#3366ff) and red (specifically, darkred).
#' @param scale_by character; whether to scale the circles such that the
#'   most-significant term in each row (\code{scale_by="row"}), column
#'   (\code{scale_by="column"}), or overall (\code{scale_by="none"}) is of
#'   maximum area. Default is "row" to better visualize patterns across
#'   contrasts.
#' @param cell_size \code{unit} object; the size of each heatmap cell. Default
#'   is \code{unit(14, "points")}.
#' @param filename character; the filename used to save the heatmap. If
#'   \code{character(0)} (default), the heatmap will be displayed instead.
#' @param height numeric; height of the file in \code{units}.
#' @param width same as \code{height}.
#' @param units character; units that define \code{height} and \code{width}.
#'   Defaults to "in" (inches).
#' @param heatmap_args list; additional arguments passed to
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param padj_args list; additional arguments passed to
#'   \code{\link[ComplexHeatmap]{Legend}}. Modifies the adjusted p-value legend.
#' @param save_args list; additional arguments passed to the graphics device.
#'   See \code{\link[grDevices]{png}} for options.
#' @param draw_args list; additional arguments passed to
#'   \code{\link[ComplexHeatmap]{draw-HeatmapList-method}}.
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom data.table `:=` `.N` setorderv setDT dcast
#' @importFrom grid gpar unit grid.circle grid.rect convertUnit
#' @importFrom utils modifyList head
#' @importFrom grDevices bmp dev.off jpeg png tiff pdf
#' @importFrom latex2exp TeX
#'
#' @export enrichmat
#'
#' @seealso \code{\link[ComplexHeatmap]{ComplexHeatmap-package}}
#'
#' @author Tyler Sagendorf

enrichmat <- function(x,
                      n_top = 15,
                      top_pathways = NULL,
                      rownames_column = "pathway",
                      NES_column = "NES",
                      padj_column = "padj",
                      padj_legend_title = padj_column,
                      padj_cutoff = 0.05,
                      plot_sig_only = TRUE,
                      padj_fill = "grey",
                      colors = c("#3366ff", "darkred"),
                      scale_by = c("row", "column", "none"),
                      cell_size = unit(14, "points"),
                      filename = character(0),
                      height = 5,
                      width = 5,
                      units = "in",
                      heatmap_args = list(),
                      padj_args = list(),
                      save_args = list(),
                      draw_args = list())
{
  scale_by <- match.arg(scale_by, scale_by)

  # Filter to n_top
  cols_to_keep <- c(rownames_column, padj_column, NES_column, "contrast")

  # Check that all required columns are present
  col_present <- cols_to_keep %in% colnames(x)
  if (any(!col_present)) {
    missing_cols <- cols_to_keep[!col_present]
    stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
  }

  setDT(x)

  x[, `:=` (any_sig = any(get(padj_column) < padj_cutoff),
            criteria = max(-log10(get(padj_column)),
                           na.rm = TRUE)),
    by = rownames_column]

  setorderv(x, cols = c("contrast", "criteria"), order = c(1, -1))

  if (is.null(top_pathways)) {
    top_pathways <- head(unique(x[, get(rownames_column)]), n = n_top)
  }
  x <- subset(x, get(rownames_column) %in% top_pathways)

  # Only plot pathways that are significantly enriched at least once
  if (plot_sig_only) { x <- subset(x, subset = any_sig == TRUE) }

  x <- unique(x[, cols_to_keep, with = FALSE])
  n <- x[, .N, by = rownames_column][["N"]] # number of entries by row name

  # Multiple rownames_column entries per contrast
  if (!all(n == max(n))) {
    stop("rownames_column is not uniquely defined for each contrast.")
  }

  x <- dcast(x, get(rownames_column) ~ contrast,
             value.var = c(NES_column, padj_column))
  setDF(x, rownames = x[[1]])
  x <- as.matrix(x[, -1])

  if (nrow(x) == 0) {
    stop("There is nothing to plot. Consider setting plot_sig_only=FALSE.")
  }

  # Matrices of adjusted p-values and NES
  padj_mat <- x[, grepl(paste0("^", padj_column), colnames(x)), drop = FALSE]
  NES_mat <- x[, grepl(paste0("^", NES_column), colnames(x)), drop = FALSE]

  colnames(padj_mat) <- colnames(NES_mat) <- contrasts <-
    sub(paste0("^", padj_column, "_"), "", colnames(padj_mat))

  colorRamp2_args <- heatmap_color_fun(NES_mat, colors)

  # Create heatmap ------------------------------------------------
  # Arguments that will be passed to ComplexHeatmap::Heatmap
  base_heatmap_args <- list(
    matrix = NES_mat,
    col = do.call(what = circlize::colorRamp2, args = colorRamp2_args),
    heatmap_legend_param = list(
      title = NES_column,
      at = colorRamp2_args$breaks,
      border = "black",
      legend_height = max(cell_size * 5, unit(21.1, "mm")),
      grid_width = max(cell_size, unit(4, "mm"))
    ),
    border = TRUE,
    row_names_max_width = max_text_width(rownames(x)),
    column_names_max_height = max_text_width(contrasts),
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    clustering_method_rows = ifelse(anyNA(x), "average", "complete"),
    clustering_method_columns = ifelse(anyNA(x), "average", "complete"),
    na_col = "black",
    height = cell_size * nrow(NES_mat),
    width = cell_size * ncol(NES_mat),
    layer_fun = layer_fun
  )

  # Update with user-supplied arguments
  heatmap_args <-  modifyList(x = base_heatmap_args,
                              val = heatmap_args,
                              keep.null = TRUE)

  # Mark missing values
  heatmap_args$rect_gp <- gpar(col = NA, fill = heatmap_args$na_col)

  # Color function for circles and NES legend
  col_fun <- heatmap_args$col

  # If layer_fun is specified, set the environment
  if (!is.null(heatmap_args$layer_fun)) {
    environment(heatmap_args$layer_fun) <- environment()
  }

  # Create heatmap
  ht <- do.call(what = Heatmap, args = heatmap_args)

  # Legend for background fill
  # base args
  lt_args <- list(
    title = padj_legend_title,
    at = 1:2,
    labels = latex2exp::TeX(c("$< 0.05$", "$\\geq 0.05$")),
    legend_gp = gpar(fill = c(padj_fill, "white")),
    grid_height = heatmap_args$heatmap_legend_param$grid_width,
    grid_width = heatmap_args$heatmap_legend_param$grid_width,
    border = "black", nrow = 2, direction = "horizontal"
  )
  lt_args <- modifyList(x = lt_args, val = padj_args, keep.null = TRUE)

  lt <- do.call(what = Legend, args = lt_args)

  if (!identical(filename, character(0))) {
    on.exit(dev.off())
    base_save_args <- list(filename = filename,
                           height = height, width = width,
                           units = units)
    save_args <- modifyList(x = base_save_args,
                            val = save_args,
                            keep.null = TRUE)
    do.call(what = save_heatmap, args = save_args)
  }
  draw_args <- modifyList(
    x = list(object = ht,
             annotation_legend_list = list(lt),
             merge_legends = TRUE,
             legend_gap = unit(0.15, "in")),
    val = draw_args, keep.null = TRUE
  )

  do.call(what = draw, args = draw_args)
}



## Helper functions ------------------------------------------------------------

## Format cells of heatmap
layer_fun <- function(j, i, x, y, w, h, f)
{
  # Cell background
  grid.rect(x = x, y = y, width = w, height = h,
            gp = gpar(col = NA,
                      fill = ifelse(pindex(padj_mat, i, j) < padj_cutoff,
                                    padj_fill, # lightgrey
                                    ifelse(is.na(pindex(padj_mat, i, j)),
                                           NA, "white"))
            ))

  # Matrix of radii (optionally scaled to row or column max)
  rmat <- -log10(padj_mat)
  if (scale_by != "none") {
    margin <- 1 + (scale_by == "column")
    rmat <- sweep(rmat, MARGIN = margin,
                  apply(rmat, MARGIN = margin, max, na.rm = TRUE),
                  FUN = "/")
  } else {
    rmat <- rmat / max(rmat, na.rm = TRUE)
  }

  grid.circle(
    x = x, y = y,
    r = pindex(rmat, i, j) / 2 * cell_size,
    gp = gpar(col = ifelse(pindex(padj_mat, i, j) < padj_cutoff,
                           "black", NA),
                           fill = col_fun(pindex(NES_mat, i, j))
    ))
}


## Determine breaks and colors for heatmap legend
heatmap_color_fun <- function(NES_mat,
                              colors = c("#3240cd", "darkred"))
{
  # Color function
  NES_mat <- NES_mat[!is.na(NES_mat)]
  extended_range <- range_extend(NES_mat, nearest = 0.1)

  # if (all(sign(NES_mat) %in% c(0, +1))) {
  if (all(NES_mat >= -1)) {
    breaks <- c(-1, 1, extended_range[2])
    colors <- c("white", "white", colors[2])
  }
  # } else if (all(sign(NES_mat) %in% c(0, -1))) {
  else if (all(NES_mat <= +1)) {
    breaks <- c(extended_range[1], -1, 1)
    colors <- c(colors[1], "white", "white")
  } else {
    breaks <- c(extended_range[1], -1, 1, extended_range[2])
    colors <- c(colors[1], rep("white", 2), colors[2])
  }
  return(list(breaks = breaks, colors = colors))
}


# Save the heatmap
save_heatmap <- function(filename = "enrichmat%03d.png",
                         width = 480, height = 480,
                         units = "px", res = 500, ...)
{
  file_ext <- sub(".*\\.(.*)", "\\1", filename)
  save_fun <- switch(file_ext,
                     "png" = png,
                     "tif" = tiff,
                     "tiff" = tiff,
                     "jpg" = jpeg,
                     "bmp" = bmp,
                     "pdf" = pdf)

  default_args <- list(filename = filename,
                       file = filename,
                       width = width,
                       height = height,
                       units = units,
                       quality = 100,
                       res = res, dpi = 500,
                       compression = "lzw")
  save_args <- modifyList(default_args, val = list(...))
  save_args <- save_args[names(save_args) %in% names(formals(save_fun))]

  do.call(what = save_fun, args = save_args)
}

