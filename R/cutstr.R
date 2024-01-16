#' @title Advanced Substrings of a Character Vector
#'
#' @param x a character vector.
#' @param split a length 1 character vector. Used to split \code{x}.
#' @param n integer; maximum length of output string.
#'
#' @returns A character vector of the same length as \code{x}.
#'
#' @examples
#' x <- "This is an example sentence"
#' # Return x as-is
#' cutstr(x)
#'
#' # Same as substr
#' cutstr(x, n = 15)
#'
#' # Break between words
#' cutstr(x, split = " ", n = 15)
#'
#' @export cutstr

cutstr <- function(x, split = "", n = Inf) {
  x <- strsplit(x, split = split)

  x <- sapply(x, function(.x) {
    keep <- cumsum(nchar(.x)) + nchar(split) * (seq_along(.x) - 1L) <= n
    .x <- paste(.x[keep], collapse = split)

    return(.x)
  })

  return(x)
}
