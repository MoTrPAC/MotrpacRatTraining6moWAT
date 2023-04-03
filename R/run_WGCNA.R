#' @title Run WGCNA Pipeline
#'
#' @description Run the WGCNA pipeline used in the MoTrPAC PASS1B WAT
#'   manuscript.
#'
#' @param object object of class
#'   \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}. \code{exprs}
#'   should be a matrix of log\eqn{_2}-transformed values. Count data is not
#'   accepted.
#' @param power integer; optional soft power vector (length 1 or greater). If
#'   not specified, the lowest value of 1--20 that satisfies scale-free fit
#'   R\eqn{^2 \geq} \code{RsquaredCut} will be used.
#' @param RsquaredCut numeric; minimum acceptable scale-free fit R\eqn{^2}.
#'   Ignored if \code{power} is provided.
#' @param module_prefix character; appended to the module numbers to create the
#'   "moduleID" column.
#'
#' @returns Object of class \code{list} of length 2:
#'
#' \itemize{
#'   \item "modules": a \code{data.frame} with the following columns, in
#'   addition to all columns in \code{Biobase::fData(object)}.
#'     \describe{
#'       \item{moduleColor}{moduleColor; unique color assigned to each module.
#'       The "grey" module always contains features that are not co-expressed.}
#'       \item{moduleID}{factor; \code{module_prefix} followed by a unique
#'       module number. The "grey" module is always 0.}
#'     }
#'   \item "MEs": a \code{data.frame} with 3 variables, in addition to all
#'   variables in \code{pData(object)}.
#'     \describe{
#'       \item{\code{moduleID}}{factor; \code{module_prefix} followed by a
#'       unique module number. The "grey" module is always 0.}
#'       \item{\code{ME}}{numeric; module eigenfeature values. One per sample x
#'       module combination.}
#'       \item{\code{moduleNum}}{integer; unique module ID number. The "grey"
#'       module is always 0.}
#'     }
#' }
#'
#' @importFrom Biobase exprs pData fData
#' @importFrom WGCNA pickSoftThreshold adjacency TOMsimilarity labels2colors
#'   plotDendroAndColors moduleEigengenes mergeCloseModules standardColors bicor
#' @importFrom stats hclust as.dist cor
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom graphics par text
#' @importFrom tibble deframe
#' @importFrom data.table setDT setDF melt setorderv `:=`
#'
#' @export run_WGCNA
#'
#' @author Zhenxin Hou, Tyler Sagendorf
#'
#' @references Zhang, B., & Horvath, S. (2005). A general framework for weighted
#'   gene co-expression network analysis. \emph{Statistical applications in
#'   genetics and molecular biology, 4}, Article17.
#'   \url{https://doi.org/10.2202/1544-6115.1128}
#'
#'   Langfelder P and Horvath S, WGCNA: an R package for weighted correlation
#'   network analysis. BMC Bioinformatics 2008, 9:559
#'   \url{https://doi.org/10.1186/1471-2105-9-559}
#'
#'   Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
#'   Correlations and Hierarchical Clustering. \emph{Journal of Statistical
#'   Software, 46}(11), 1--17. \url{http://www.jstatsoft.org/v46/i11/}
#'
#'   Horvath, S. (2011). Weighted Network Analysis. \emph{Springer New York}.
#'   \url{https://doi.org/10.1007/978-1-4419-8819-5}
#'
#'   Langfelder P, Zhang B, Horvath wcfS (2016). \emph{dynamicTreeCut: Methods
#'   for Detection of Clusters in Hierarchical Clustering Dendrograms}. R
#'   package version 1.63-1,
#'   \url{https://CRAN.R-project.org/package=dynamicTreeCut}.

run_WGCNA <- function(object,
                      power = 1:20,
                      RsquaredCut = 0.90,
                      module_prefix = "") # "P", "M", "T"
{
  on.exit(invisible(gc()))

  # Transpose for WGCNA
  datExpr <- t(exprs(object))

  # # Check if correlations can be calculated for all pair of features
  # if (anyNA(datExpr %*% t(datExpr))) {
  #   stop(paste("Too many missing values to calculate correlations.",
  #              "Consider filtering to features with fewer than 50% missing ",
  #              "values or imputing missing values.", sep = "\t"))
  # }

  # Choose soft-thresholding power
  if (length(power) > 1) {
    message("Choosing soft-thresholding power...")
    # Choosing the soft-thresholding power
    sft <- pickSoftThreshold(data = datExpr,
                             powerVector = power,
                             verbose = 5,
                             RsquaredCut = RsquaredCut,
                             networkType = "signed",
                             corFnc = bicor,
                             corOptions = list(use = "pairwise.complete.obs"))

    # Plot the Scale-free topology fit index and Mean connectivity
    # sizeGrWindow(9, 5)
    par(mfrow = c(1, 2))
    cex1 <- 0.9

    plot(sft$fitIndices[, 1],
         -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
         xlab = "Soft Threshold (power)",
         ylab = "Scale Free Topology Model Fit, signed R^2",
         type = "n",
         main = "Scale independence")
    text(sft$fitIndices[, 1],
         -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
         labels = power, cex = cex1, col = "red")
    abline(h = 0.90, col = "red")

    plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
         xlab = "Soft Threshold (power)",
         ylab = "Mean Connectivity",
         type = "n",
         main = "Mean connectivity")
    text(sft$fitIndices[, 1], sft$fitIndices[, 5],
         labels = power, cex = cex1, col = "red")

    # Power estimate
    power <- sft$powerEstimate

    if (is.na(power)) {
      stop(paste("No values of powerVector satisfy SFT.R.sq >=",
                 RsquaredCut))
    } else {
      msg <- "The estimated soft-thresholding power is %d."
      message(sprintf(msg, power))
    }
  } else {
    msg <- "Using soft-thresholding power of %d."
    message(sprintf(msg, power))
  }

  message("Constructing adjacency matrix...")
  ## Adjacency matrix
  adjacency <- adjacency(datExpr = datExpr,
                         power = power,
                         type = "signed",
                         corFnc = bicor,
                         corOptions = list(use = "pairwise.complete.obs"))

  message("Constructing dissimilarity matrix from topological overlap...")
  dissTOM <- 1 - TOMsimilarity(adjMat = adjacency, TOMType = "unsigned")

  # Hierarchical clustering based on dissTOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  # Module identification using dynamic tree cut:
  dynamicMods <- cutreeDynamic(dendro = geneTree,
                               distM = dissTOM,
                               cutHeight = 0.99,
                               deepSplit = 2,
                               pamRespectsDendro = FALSE,
                               minClusterSize = 30)

  # Convert numeric labels to colors
  moduleColors <- labels2colors(labels = dynamicMods)

  dynamicMods <- factor(dynamicMods, levels = sort(unique(dynamicMods)))

  # Plot the dendrogram
  par(mfrow = c(1, 1))
  plot(geneTree, xlab = "", sub = "",
       main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)

  plotDendroAndColors(dendro = geneTree,
                      colors = moduleColors,
                      groupLabels =  "Module Color",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")

  # Calculate eigengenes --------------------
  MEList <- moduleEigengenes(expr = datExpr,
                             colors = moduleColors,
                             softPower = power) ## updated 2022-09-06
  MEs <- MEList$eigengenes

  # Eigengene dissimilarity clustering
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")

  # Plot the result
  plot(METree,
       main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  MEDissThres <- 0.15 # 0.85 correlation
  abline(h = MEDissThres, col = "green")

  # Merge highly correlated (1-MEDissThres) modules ----
  message("Merging highly-correlated modules...")
  # Call an automatic merging function
  merge <- mergeCloseModules(exprData = datExpr,
                             colors = moduleColors,
                             MEs = MEs,
                             cutHeight = MEDissThres,
                             verbose = 3,
                             relabel = TRUE)
  # The merged module colors
  moduleColors <- merge$colors
  # Eigengenes of the new merged modules:
  MEs <- merge$newMEs

  # Update dynamicMods
  color_lvls <- c("grey", standardColors(n = 50))
  color_lvls <- color_lvls[color_lvls %in% moduleColors]
  dynamicMods <- factor(moduleColors, levels = color_lvls)
  dynamicMods <- as.numeric(dynamicMods) - ("grey" %in% moduleColors)

  ## Reformat results ----
  modules <- data.frame(moduleColor = moduleColors,
                        moduleID = paste0(module_prefix, dynamicMods))
  color2id <- deframe(unique(modules))
  modules <- cbind(modules, fData(object))

  colnames(MEs) <- color2id[sub("^ME", "", colnames(MEs))]

  ME_long <- pData(object)
  id.vars <- colnames(ME_long)
  setDT(ME_long)
  ME_long <- cbind(ME_long, MEs)

  ME_long <- melt(ME_long, id.vars = id.vars,
                  variable.name = "moduleID",
                  value.name = "ME")
  ME_long[, `:=` (moduleNum = as.numeric(
    sub(paste0("^", module_prefix), "", moduleID)
  ))]
  setorderv(ME_long, cols = "moduleNum")
  ME_long[, `:=` (moduleID = factor(moduleID, levels = unique(moduleID)))]
  setDF(ME_long)
  MEs <- ME_long

  modules$moduleID <- factor(modules$moduleID, levels = levels(MEs$moduleID))

  return(list(modules = modules, MEs = MEs))
}

