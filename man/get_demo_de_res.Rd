% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gprofiler2_analysis.R
\name{get_demo_de_res}
\alias{get_demo_de_res}
\title{Imports a standard resulting data.frame from DESeq2}
\usage{
get_demo_de_res(number = 1)
}
\arguments{
\item{nombre}{an \code{integer} that can take the value 1 or 2 specifying if you want to
import only one \code{data.frame} or a \code{list} of two data.frames
Default: 1}
}
\value{
if no argument given the object returned is a .csv file that contains a \code{data.frame} named df1 with:
df1 specifications:
25 rows with a padj < 0.05 (significant) & a abs(log2FoldChange) > log2(1.5)
25 rows with a padj < 0.05 (significant) & a abs(log2FoldChange) < log2(1.5)
25 rows with a padj > 0.05 (not significant) & a abs(log2FolChcange) > log2(1.5)
25 rows with a padj > 0.05 (not significant) & a abs(log2FoldChange) < log2(1.5)
if argument given (number=2) the object returned is a \code{list} of two .csv files that contain a \code{data.frame} named df1 and df2.
df2 specifications:
50 rows with a padj < 0.05 (significant) & a abs(log2FoldChange) > log2(2)
50 rows with a padj < 0.05 (significant) & a abs(log2FoldChange) < log2(2)
50 rows with a padj > 0.05 (not significant) & a abs(log2FolChcange) > log2(2)
50 rows with a padj > 0.05 (not significant) & a abs(log2FoldChange) < log2(2)
}
\description{
Imports a standard resulting data.frame from DESeq2
}
\examples{
de_res <- get_demo_de_res()
de_res <- get_demo_de_res(2)

}
