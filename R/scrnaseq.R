#' Read Cell Ranger matrix and return a tibble
#' Code adapted from https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices#r-load-mat
#' @import Matrix magrittr tidyr Seurat hdf5r
#' @export
#' @param matrix_dir Path of directory containing the files "matrix.mtx.gz", "features.tsv.gz" and "barcodes.tsv.gz"
#' @param tidy Return tidy tibble
#' @return tibble of resulting matrix
read_cellranger_feature_bc_matrix <- function(matrix_dir, tidy = TRUE) {
  barcode_path <- file.path(matrix_dir, "barcodes.tsv.gz")
  features_path <- file.path(matrix_dir, "features.tsv.gz")
  matrix_path <- file.path(matrix_dir, "matrix.mtx.gz")
  mat <- Matrix::readMM(file = matrix_path)
  feature_names <-
    read.delim(
      features_path,
      header = FALSE,
      stringsAsFactors = FALSE
    )
  barcode_names <-
    read.delim(barcode_path,
      header = FALSE,
      stringsAsFactors = FALSE
    )
  colnames(mat) <- barcode_names$V1
  rownames(mat) <- feature_names$V1

  mat %>%
    cellranger_matrix_to_tibble(tidy = tidy) %>%
    return()
}

#' Helper to convert cellranger data into a tidy tibble
#' @param mat dgTMatrix object created from cellranger
#' @param tidy Return tidy tibble
#' @return tibble of matrix object
cellranger_matrix_to_tibble <- function(mat, tidy) {
  mat %<>%
    # Convert dgTMatrix to normal matrix
    as.matrix() %>%
    # Keep rownames
    tibble::as_tibble(rownames = NA)

  if(tidy) {
    # Make tidy dataframe
    mat %<>%
      tibble::rownames_to_column(var = "id") %>%
      tidyr::pivot_longer(
        -id,
        names_to = "cell_barcode",
        values_to = "count"
      )
  }

  return(mat)
}

#' Read Cell Ranger matrix and return a tibble
#' Code adapted from https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices#r-load-mat
#' @import Matrix magrittr tidyr
#' @export
#' @param h5_path Path of file
#' @param use_names Keep original ensembl IDs or use clear gene names
#' @param tidy Return tidy tibble
#' @return tibble of resulting matrix
read_cellranger_feature_bc_matrix_h5 <- function(h5_path, use_names = FALSE, tidy = TRUE) {
  Seurat::Read10X_h5(h5_path, use.names = use_names) %>%
    cellranger_matrix_to_tibble(tidy = tidy) %>%
    return()
}

#' Seurat clustering based on https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#' @param seurat_data Seurat dataset object
#' @param max_dimension_used_for_clustering Determines how many principal components are used for clustering
#' @param cluster_resolution See Seurat::FindClusters resolution parameter. "resolution: Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
#' @import magrittr Seurat
#' @export
#' @examples
#' cellatlas::get_campbell_count_data() %>%
#'   # transform into format that can be read by seurat
#'   dplyr::select(GeneID, CellID, count) %>%
#'   tidyr::spread(CellID, count) %>%
#'   tibble::column_to_rownames("GeneID") %>%
#'   # Initialize the Seurat object with the raw (non-normalized data)
#'   Seurat::CreateSeuratObject(
#'     counts = .,
#'     min.cells = 3,
#'     min.features = 200,
#'     project = "campbell"
#'   ) %>%
#'   # Get seurat clustering
#'   get_seurat_clustering()
get_seurat_clustering <-
  function(
    seurat_data,
    max_dimension_used_for_clustering = 10,
    cluster_resolution = 0.8
  ) {
  seurat_data %<>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData(., features = rownames(.)) %>%
    Seurat::RunPCA() %>%
    # Clustering
    Seurat::FindNeighbors(dims = 1:max_dimension_used_for_clustering) %>%
    Seurat::FindClusters(resolution = cluster_resolution) %>%
    Seurat::RunUMAP(dims = 1:max_dimension_used_for_clustering) %>%
    return()
}


#' Get cell counts for list of genes
#' @param seurat_dat Seurat object
#' @param genes List of genes to extract from Seurat object
#' @param tidy Return column "gene_id" or make each gene a separate column
#' @return tibble of counts per cell
#' @import magrittr Seurat tidyr dplyr tibble
#' @export
get_gene_counts_per_cell <- function(seurat_dat, genes, tidy = TRUE) {
  # Extract the genes from the list of cells
  dat <- Seurat::GetAssayData(seurat_dat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(-gene_id, names_to = "cell_id", values_to = "count") %>%
    dplyr::filter(gene_id %in% genes)
  
  if (!tidy)
    # Pivot wider again to make filtering easier
    dat %<>% tidyr::pivot_wider(names_from = gene_id, values_from = count)
  
  return(dat)
}
