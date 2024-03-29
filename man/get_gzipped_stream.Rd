% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{get_gzipped_stream}
\alias{get_gzipped_stream}
\title{Get stream for gzipped content}
\usage{
get_gzipped_stream(url)
}
\arguments{
\item{url}{Url to call}
}
\value{
stream for piping to read functions
}
\description{
Get stream for gzipped content
}
\examples{
"https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession\%2Creviewed\%2Cid\%2Cprotein_name\%2Cgene_names\%2Corganism_name&format=tsv&query=\%28G-coupled\%20coupled\%20receptor\%29\%20AND\%20\%28model_organism\%3A10090\%29" \%>\%
   url() \%>\%
   gzcon() \%>\%
   readr::read_tsv()
}
