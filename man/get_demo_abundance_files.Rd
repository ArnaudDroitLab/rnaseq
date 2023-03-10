% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demo.R
\name{get_demo_abundance_files}
\alias{get_demo_abundance_files}
\title{Get demo kallisto abundance files}
\usage{
get_demo_abundance_files(large = FALSE)
}
\arguments{
\item{large}{If \code{TRUE}, returns 8 filenames. Otherwise returns 4.}
}
\value{
A vector of kallisto abundance filenames
}
\description{
Get demo kallisto abundance files
}
\examples{
abundances <- get_demo_abundance_files()

}
