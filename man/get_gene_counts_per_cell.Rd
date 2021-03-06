% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrnaseq.R
\name{get_gene_counts_per_cell}
\alias{get_gene_counts_per_cell}
\title{Get cell counts for list of genes}
\usage{
get_gene_counts_per_cell(seurat_dat, genes, tidy = TRUE)
}
\arguments{
\item{seurat_dat}{Seurat object}

\item{genes}{List of genes to extract from Seurat object}

\item{tidy}{Return column "gene_id" or make each gene a separate column}
}
\value{
tibble of counts per cell
}
\description{
Get cell counts for list of genes
}
