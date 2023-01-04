#' @title LIMMA Convenience Wrapper
#'
#' @description A convenience wrapper for differential analysis with LIMMA.
#'   Performs moderated t-tests, F-tests, or linear regression.
#'
#' @param object Object of class \code{\link[MSnbase:MSnSet-class]{MSnSet}}. The
#'   \code{exprs} slot can be either a matrix of log\eqn{_2} relative abundances
#'   or counts. If the latter, the limma--edgeR pipeline
#'   (\url{https://doi.org/10.12688/f1000research.9005.3}) will automatically be
#'   used.
#' @param model.str character; formulation of the model (e.g. \code{"~ a + b"}).
#'   If \code{contrasts} are provided, this should be a no-intercept model (it
#'   should include 0 or -1 as terms). See \code{\link[stats]{lm}} for more
#'   details.
#' @param coef.str character; coefficient of interest. One of
#'   \code{colnames(pData(object))}.
#' @param contrasts character; optional contrasts of the form \code{"a-b"} to
#'   test.
#' @param trend logical; should an intensity-dependent trend be allowed for the
#'   prior variance? If \code{FALSE}, then the prior variance is constant.
#'   Alternatively, \code{trend} can be a row-wise numeric vector, which will be
#'   used as the covariate for the prior variance. See
#'   \code{\link[limma]{eBayes}} for more details.
#' @param robust logical; should the estimation of \code{df.prior} and
#'   \code{var.prior} be robustified against outlier sample variances? See
#'   \code{\link[limma]{eBayes}} for more details.
#' @param var.group character; the column in \code{pData(object)} indicating
#'   groups that will have different relative quality weights. The default
#'   (\code{character(0)}) will weight each sample equally (all samples given a
#'   weight of 1). Weights affect the logFC and P.Value columns of the results.
#' @param block \code{NULL} or character; name of a column in
#'   \code{pData(object)} specifying a blocking variable. Passed to
#'   \code{\link[limma]{duplicateCorrelation}}. See
#'   \href{https://support.bioconductor.org/p/125489/#125602}{When to use
#'   \code{duplicateCorrelation}}. Section 9.7 "Multi-level Experiments" of the
#'   \href{https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf}{LIMMA
#'   User's Guide} explains when to use \code{duplicateCorrelation}. Currently
#'   ignored if \code{exprs(object)} is a matrix of counts.
#' @param plot logical; whether to generate diagnostic plots. If \code{TRUE},
#'   generates a barplot of weights applied to each sample and a plot of the
#'   mean-variance trend using \code{\link[limma]{plotSA}}.
#' @param adjust.method method for p-value adjustment. Default is \code{"BH"}
#'   (Benjamini-Hochberg), which will control the false discovery rate (FDR).
#'   See \code{\link[stats]{p.adjust}} for details.
#' @param adjust.globally logical; should p-values from all contrasts be
#'   adjusted together using \code{adjust.method}? Set to \code{FALSE} if the
#'   contrasts being tested are not closely related. See the "Multiple Testing
#'   Across Contrasts" section of the LIMMA User's Guide (linked in "See Also")
#'   for more information.
#' @param ... additional arguments passed to \code{\link[limma]{plotSA}}.
#'
#' @return \code{data.frame}. Output of \code{\link[limma]{topTable}} with
#'   additional columns \code{feature}, \code{contrast}, and column(s) for the
#'   standard error (se) of the logFC. All columns from \code{fData} are also
#'   included.
#'
#'
#' @details An MDS plot (produced by \code{\link[limma]{plotMDS}}) is used to
#'   determine the appropriate value of \code{var.group}. If samples within
#'   phenotype groups tend to cluster well and have similar dispersions, the
#'   default \code{var.group = character(0)} is recommended. If samples within
#'   phenotype groups tend to cluster well, but different groups are more or
#'   less dispersed, it may be beneficial to set \code{var.group} to the name of
#'   the phenotype column. If one or more samples tend to cluster poorly with
#'   samples of the same phenotype, it may be beneficial to weight by sample.
#'   That is, set \code{var.group} to the name of the column in
#'   \code{pData(object)} that uniquely identifies each sample. If variation
#'   between samples or groups of samples is biological rather than technical,
#'   weighting is not recommended.
#'
#'   The plot of the mean-variance trend helps determine whether to fit a trend
#'   line to the prior variance (default is \code{trend = TRUE}). It also helps
#'   determine if the prior variance should be robustified against outlier
#'   sample variances (\code{robust = TRUE}). See \code{\link[limma]{eBayes}}
#'   for more details.
#'
#' @seealso
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf}{LIMMA
#' User's Guide}
#'
#' \href{https://online.stat.psu.edu/stat555/node/1/}{PennState STAT 555 -
#' Statistical Analysis of Genomics Data}
#'
#' \href{http://varianceexplained.org/statistics/interpreting-pvalue-histogram/}{How
#' to interpret a p-value histogram}
#'
#' Breheny, P., Stromberg, A., & Lambert, J. (2018). \emph{p}-Value Histograms:
#' Inference and Diagnostics. \emph{High-throughput, 7}(3), 23.
#' \url{https://doi.org/10.3390/ht7030023}
#'
#'
#' @references Smyth G. K. (2004). Linear models and empirical bayes methods for
#'   assessing differential expression in microarray experiments.
#'   \emph{Statistical applications in genetics and molecular biology, 3},
#'   Article3. \url{https://doi.org/10.2202/1544-6115.1027}
#'
#'   Smyth, G. K., Michaud, J., and Scott, H. (2005). The use of within-array
#'   replicate spots for assessing differential expression in microarray
#'   experiments. \emph{Bioinformatics 21}(9), 2067--2075.
#'   \url{http://bioinformatics.oxfordjournals.org/content/21/9/2067} Preprint
#'   with corrections: \url{http://www.statsci.org/smyth/pubs/dupcor.pdf}
#'
#'   Liu, R., Holik, A. Z., Su, S., Jansz, N., Chen, K., Leong, H. S., Blewitt,
#'   M. E., Asselin-Labat, M. L., Smyth, G. K., & Ritchie, M. E. (2015). Why
#'   weight? Modelling sample and observational level variability improves power
#'   in RNA-seq analyses. \emph{Nucleic acids research, 43}(15), e97.
#'   \url{https://doi.org/10.1093/nar/gkv412}
#'
#'   Law, C. W., Alhamdoosh, M., Su, S., Dong, X., Tian, L., Smyth, G. K., &
#'   Ritchie, M. E. (2016). RNA-seq analysis is easy as 1-2-3 with limma, Glimma
#'   and edgeR. \emph{F1000Research, 5}, ISCB Comm J-1408.
#'   \url{https://doi.org/10.12688/f1000research.9005.3}
#'
#'   Phipson, B., Lee, S., Majewski, I. J., Alexander, W. S., & Smyth, G. K.
#'   (2016). ROBUST HYPERPARAMETER ESTIMATION PROTECTS AGAINST HYPERVARIABLE
#'   GENES AND IMPROVES POWER TO DETECT DIFFERENTIAL EXPRESSION. \emph{The
#'   annals of applied statistics, 10}(2), 946--963.
#'   \url{https://doi.org/10.1214/16-AOAS920}
#'
#' @importFrom MSnbase exprs pData
#' @importFrom stats terms model.matrix p.adjust
#' @importFrom data.table rbindlist
#' @importFrom graphics barplot abline
#' @importFrom edgeR DGEList filterByExpr calcNormFactors cpm
#' @import limma
#' @import statmod
#'
#' @export limma_full
#'
#' @author Tyler Sagendorf, Vlad Petyuk

limma_full <- function(object,
                       model.str,
                       coef.str,
                       contrasts,
                       trend = TRUE,
                       robust = TRUE,
                       var.group = character(0),
                       block = NULL,
                       plot = FALSE,
                       adjust.method = "BH",
                       adjust.globally = TRUE,
                       ...)
{
  # Remove samples with all missing values
  na.idx.sample <- colMeans(is.na(exprs(object))) == 1
  if (sum(na.idx.sample) > 0) {
    message(paste(sprintf("There are %d samples with all missing values.",
                          sum(na.idx.sample)), "These will be removed.",
                  sep = "\n"))
    object <- object[, !na.idx.sample]
  }

  # Remove features with all missing values
  na.idx <- rowMeans(is.na(exprs(object))) == 1
  if (sum(na.idx) > 0) {
    message(paste(sprintf("There are %d features with all missing values.",
                          sum(na.idx)), "These will be removed.", sep = "\n"))
    object <- object[!na.idx, ]
  }

  model.formula <- eval(parse(text = model.str), envir = pData(object))

  if (!missing(contrasts)) {
    if (attr(terms(model.formula), which = "intercept") != 0) {
      stop(paste("Please specify a no-intercept model for model.str",
                 "if contrasts are provided. See ?lm for more details.",
                 sep = "\n"))
    }
  }

  # Create design matrix and remove coef.str from the beginning of the levels
  # (easier for user when specifying contrasts)
  design <- model.matrix(model.formula)
  if (!is.numeric(pData(object)[[coef.str]])) {
    colnames(design) <- sub(paste0("^", coef.str), "", colnames(design))
  }

  if (!is.null(nonEstimable(design))) {
    stop(paste("One or more coefficients are not estimable.",
               "The design matrix is not full rank."))
  }

  # Some samples may be NA, so they will be dropped in the design matrix.
  # Subset object to only those rows present in the design matrix.
  object <- object[, as.numeric(rownames(design))]

  # By default, sample-level weights are all 1
  weights <- rep(1, ncol(object))
  names(weights) <- seq_along(weights)

  # Separate processing if dealing with RNA-seq (count) data
  if (isTRUE(all.equal(exprs(object),
                       floor(exprs(object))))) {
    message("Count data detected. Using the RNA-Seq pipeline.")
    # Convert to DGEList
    dge <- DGEList(counts = exprs(object),
                   group = pData(object)[[coef.str]],
                   samples = pData(object),
                   genes = fData(object))

    ## Filter
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    message(
      sprintf("%d low count features were removed. %d features remaining.",
              sum(!keep), sum(keep))
    )

    if (!identical(var.group, character(0))) {
      if (length(keep) < 10) {
        message("Fewer than 10 features. var.group will be ignored.")
        var.group <- character(0)
      }
    }

    # Calculate normalization factors
    dge <- calcNormFactors(dge, method = "TMM")

    # MDS plot
    if (plot) {
      lcpm <- cpm(dge$counts, log = TRUE,
                  prior.count = 0.5)
      plotMDS(lcpm, label = dge$samples[[coef.str]])
    }

    # voom (convert to log2 CPM and calculate observational-level weights)
    if (!identical(var.group, character(0))) {
      y <- voomWithQualityWeights(dge, design = design,
                                  var.group = pData(object)[[var.group]],
                                  plot = plot, method = "genebygene")
      weights <- y$targets$sample.weights
      names(weights) <- seq_along(weights)
    } else {
      y <- voom(dge, design = design, plot = plot)
    }
    fit <- lmFit(y, design)
    trend <- FALSE
  } else {
    # Non RNA-seq data -----
    if (!identical(var.group, character(0))) {
      if (nrow(object) < 10) {
        message("Fewer than 10 features. var.group will be ignored.")
      } else {
        weights <- arrayWeights(object, design, method = "genebygene",
                                var.group = pData(object)[[var.group]])
      }
    }

    # Linear modeling
    if (!is.null(block)) {
      # dupecor currently only used for non count data
      block <- eval(parse(text = block), envir = pData(object))
      dupecor <- duplicateCorrelation(object, design, block = block,
                                      weights = weights)

      fit <- lmFit(object, design, weights = weights, block = block,
                   correlation = dupecor$consensus.correlation)
    } else {
      fit <- lmFit(object, design, weights = weights)
    }
    # If not enough data, set robust to FALSE
    if (nrow(object) < 10 & robust) {
      message("Fewer than 10 features. Using trend=FALSE, robust=FALSE instead.")
      trend <- robust <- FALSE
    }
  }

  # Contrasts
  if (!missing(contrasts)) {
    contrast.matrix <- makeContrasts(contrasts = contrasts,
                                     levels = design)
    fit <- contrasts.fit(fit, contrast.matrix)
  }

  # Empirical bayes moderation
  fit.smooth <- eBayes(fit, trend = trend, robust = robust)

  ## Diagnostic plots
  if (plot) {
    # Sample weights
    barplot(weights, space = 0,
            xlab = "Sample Index", ylab = "Weight",
            main = "Sample-specific weights")
    abline(1, 0, lty = "longdash", col = "red")

    # Plot of mean-variance trend
    plotSA(fit.smooth, ...)
  }

  # Moderated t-test for each contrast (if provided)
  # or moderated F-tests
  if (!missing(contrasts)) {
    # topTable for each contrast
    res <- lapply(contrasts, function(contrast_i) {
      x <- topTable(fit.smooth,
                    coef = contrast_i,
                    number = Inf,
                    p.value = 1,
                    sort.by = "none",
                    adjust.method = adjust.method)
      x$feature <- rownames(x)
      x$contrast <- contrast_i

      ## Column for logFC standard error
      se <- sqrt(fit.smooth$s2.post) * fit.smooth$stdev.unscaled
      se <- se[rownames(x), contrast_i, drop = FALSE]
      colnames(se) <- "se"
      x <- cbind(x, se)

      return(x)
    })
    res <- rbindlist(res)
    res$contrast <- factor(res$contrast, levels = contrasts)

    # Adjust p-values across all contrasts. Contrasts should be related.
    if (adjust.globally) {
      res$adj.P.Val <- p.adjust(res$P.Value, method = adjust.method)
    }
  } else {
    # If coef.str is a factor or character, do this
    if (!(coef.str %in% colnames(design))) {
      idx <- which(names(attr(design, "contrast")) == coef.str)
      idx <- attr(design, "assign") == idx
      coef.str <- colnames(design)[idx]
    }

    res <- topTable(fit.smooth,
                    coef = coef.str,
                    number = Inf,
                    p.value = 1,
                    sort.by = "none",
                    adjust.method = adjust.method)

    ## Column(s) for standard error of logFC (se)
    se <- sqrt(fit.smooth$s2.post) * fit.smooth$stdev.unscaled
    se <- se[rownames(res), coef.str, drop = FALSE]
    colnames(se) <- paste0("se.", colnames(se))
    res <- cbind(res, se)
  }

  return(res)
}


