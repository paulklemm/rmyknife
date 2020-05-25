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

#' WriteXLS::WriteXLS() function with sane default parameter
#' @param dat data frame to export
#' @param ExcelFileName see WriteXLS::WriteXLS()
#' @param SheetNames see WriteXLS::WriteXLS()
#' @export
#' @import WriteXLS
write_xls <- function(dat, ExcelFileName, SheetNames = filename) {
  WriteXLS::WriteXLS(
    x = dat,
    ExcelFileName = file.path("low_salmon_gene_counts.xlsx"),
    AdjWidth = TRUE,
    AutoFilter = TRUE,
    BoldHeaderRow = TRUE,
    FreezeRow = 1,
    SheetNames = sheetnames
  )
}