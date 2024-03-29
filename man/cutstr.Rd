% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutstr.R
\name{cutstr}
\alias{cutstr}
\title{Advanced Substrings of a Character Vector}
\usage{
cutstr(x, split = "", n = Inf)
}
\arguments{
\item{x}{a character vector.}

\item{split}{a length 1 character vector. Used to split \code{x}.}

\item{n}{integer; maximum length of output string.}
}
\value{
A character vector of the same length as \code{x}.
}
\description{
Advanced Substrings of a Character Vector
}
\examples{
x <- "This is an example sentence"
# Return x as-is
cutstr(x)

# Same as substr
cutstr(x, n = 15)

# Break between words
cutstr(x, split = " ", n = 15)

}
