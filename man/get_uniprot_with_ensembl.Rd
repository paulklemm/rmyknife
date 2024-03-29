% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomart.R
\name{get_uniprot_with_ensembl}
\alias{get_uniprot_with_ensembl}
\title{Wrapper for get_uniprot that attaches the ensembl_gene_id and external_gene_name rows}
\usage{
get_uniprot_with_ensembl(ensembl, species = "MUS", ...)
}
\arguments{
\item{ensembl}{Ensembl database connection object}

\item{species}{Either "HUM" or "MUS". Leave empty for pulling all species}
}
\description{
Wrapper for get_uniprot that attaches the ensembl_gene_id and external_gene_name rows
}
\examples{
  gpcrs <-
    rmyknife::get_uniprot_with_ensembl(
      ensembl = rmyknife::get_ensembl_dataset_from_version(101, species = "MUS"),
      species = "MUS",
      search_string = "G-protein coupled receptor"
    ) \%>\%
    dplyr::filter(uniprot_reviewed == "reviewed") \%>\%
    dplyr::filter(!is.na(ensembl_gene_id)) \%>\%
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
}
