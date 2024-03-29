% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomart.R
\name{get_ensembl_dataset_from_version}
\alias{get_ensembl_dataset_from_version}
\title{Get Ensembl dataset with default parameter}
\usage{
get_ensembl_dataset_from_version(ensembl_version = 103, species = "MUS")
}
\arguments{
\item{ensembl_version}{Ensembl version you want to get data for}

\item{species}{Species you want to extract ("MUS", "RNO" or "HUM")}
}
\description{
Get Ensembl dataset with default parameter
}
\examples{
   ensembl <- get_ensembl_dataset_from_version(103, "MUS")
}
