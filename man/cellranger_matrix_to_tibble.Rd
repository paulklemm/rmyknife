% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrnaseq.R
\name{cellranger_matrix_to_tibble}
\alias{cellranger_matrix_to_tibble}
\title{Helper to convert cellranger data into a tidy tibble}
\usage{
cellranger_matrix_to_tibble(mat, tidy)
}
\arguments{
\item{mat}{dgTMatrix object created from cellranger}

\item{tidy}{Return tidy tibble}
}
\value{
tibble of matrix object
}
\description{
Helper to convert cellranger data into a tidy tibble
}
