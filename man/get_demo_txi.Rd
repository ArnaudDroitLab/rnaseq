% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demo.R
\name{get_demo_txi}
\alias{get_demo_txi}
\title{Get demo txi file}
\usage{
get_demo_txi(large = FALSE, txOut = FALSE, ercc = FALSE)
}
\arguments{
\item{large}{If \code{TRUE}, txi matrices will contain 8 samples, otherwise
they will contain 4. Default: \code{FALSE}}

\item{txOut}{Return counts and abundance at the transcript level. Default:
\code{FALSE}}

\item{ercc}{Add ERCC counts to the txi?}
}
\value{
A txi object
}
\description{
Get demo txi file
}
\examples{
txi <- get_demo_txi()

}
