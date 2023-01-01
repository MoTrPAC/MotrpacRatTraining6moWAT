#' @title Extended Range of Values
#'
#' @description `range_extend` returns a vector containing the minimum and
#'   maximum of all the given arguments rounded outward to the value provided by
#'   `nearest`.
#'
#' @param x any \code{\link[base]{numeric}} object.
#' @param nearest numeric; the range of \code{x} will be rounded out to the
#'   value specified by \code{nearest}. Default is 1.
#'
#' @returns A \code{\link[base]{numeric}} vector of length 2 specifying the
#'   range of values after rounding outward to the value provided by `nearest`.
#'
#' @seealso \code{\link[base]{range}}
#'
#' @export range_extend
#'
#' @author Tyler Sagendorf
#'
#' @examples
#' set.seed(0)
#' x <- runif(5, min = -10, max = 10)
#' range(x) # -4.689827  8.164156
#'
#' range_extend(x) # -5  9
#' range_extend(x, nearest = 0.1) # -4.7  8.2

range_extend <- function(x, nearest = 1) {
  x <- x / nearest
  xmin <- min(x, na.rm = T)
  xmax <- max(x, na.rm = T)
  c(floor(xmin), ceiling(xmax)) * nearest
}

