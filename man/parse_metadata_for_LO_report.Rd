% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_metadata.R
\name{parse_metadata_for_LO_report}
\alias{parse_metadata_for_LO_report}
\title{Parse metadata LO}
\usage{
parse_metadata_for_LO_report(
  metadata,
  pca_subset = "Cell",
  pca_batch_metadata = c("Cell", "Compound"),
  extra_count_matrix = NULL,
  report_title = "title"
)
}
\arguments{
\item{metadata}{(dataframe)}

\item{pca_subset}{(character) column of metadata, filter for automated PCA (ex.: Cell type)}

\item{pca_batch_metadata}{(character) extra columns for pca coloring (biological or batch effects)}

\item{extra_count_matrix}{(character)}

\item{report_title}{(character)}
}
\value{
a list of data.frame
}
\description{
Parse metadata and generate a draft for pca_info, volcano_info, report_info
files.
}
