#' Read in GTF file
#'
#' @param path Path to GTF file
#' @import magrittr rtracklayer tibble
#' @export
#' @return Tibble of GTF file
read_gtf <- function(path) {
  path %>%
    rtracklayer::import() %>%
    tibble::as.tibble() %>%
    return()
}

#' Read in data from cufflinks tsv
#'
#' @param path Path to cufflinks tsv file
#' @import readr magrittr
#' @export
#' @return Tibble of cufflinks data
read_cufflinks <- function(path) {
  path %>%
    readr::read_tsv() %>%
    return()
}

#' Read Cell Ranger matrix and return a tibble
#' Code adapted from https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices#r-load-mat
#' @import Matrix magrittr tidyr
#' @export
#' @param matrix_dir Path of files "matrix.mtx.gz", "features.tsv.gz" and "barcodes.tsv.gz"
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

  mat %<>%
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
