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

#' Print standard DT::datatable with extensions for exporting the data
#' @import DT magrittr
#' @export
#' @param dat Dataframe to print
dt_datatable <- function(dat) {
  DT::datatable(
    data = dat,
    extensions = "Buttons",
    options = list(
      dom = "Bfrtip",
      buttons = c("copy", "csv", "excel", "pdf", "print")
    )
  ) %>%
    return()
}