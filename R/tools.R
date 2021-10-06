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
#' @param caption Caption of datatable
#' @param scroll_y Pixel size of data-table size in y
dt_datatable <- function(dat, caption = "", scroll_y = 300) {
  DT::datatable(
    caption = caption,
    data = dat,
    extensions = c("Scroller", "Buttons"),
    filter = list(position = "top"),
    options = list(
      dom = "Bfrtip",
      buttons = c("copy", "csv", "excel", "pdf", "print"),
      scrollX = TRUE,
       # Makes the table more responsive when it's really big
      scrollY = scroll_y,
      scroller = TRUE
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
write_xls <- function(
  dat,
  ExcelFileName,
  SheetNames = basename(ExcelFileName)
) {
  WriteXLS::WriteXLS(
    x = dat,
    ExcelFileName = ExcelFileName,
    AdjWidth = TRUE,
    AutoFilter = TRUE,
    BoldHeaderRow = TRUE,
    FreezeRow = 1,
    SheetNames = SheetNames
  )
}

#' Print message if rmyknife.verbose is TRUE
#' @param message Debug message to print
debug_message <- function(message) {
  if (isTRUE(getOption("rmyknife.verbose"))) {
    message(message)
  }
}

#' Convert a DESeq2 results object to a tibble
#' @param dat DESeq2::results object
#' @param row_name The column name for the DESeq2 results row names
#' @return tibble of DESeq2 result
#' @import tibble magrittr
#' @export
#' @examples
#'   deseq2_diff <- DESeq2::results(deseq2_ip_wt_vs_input_wt) %>%
#'     deseq2_result_to_tibble()
deseq2_result_to_tibble <- function(dat, row_name = "ensembl_gene_id") {
  dat %>%
    as.data.frame() %>%
    tibble::rownames_to_column(row_name) %>%
    tibble::as_tibble()
}

#' Downloads the file in the URL to a temp file and returns it's path
#' @param url File url
#' @export
#' @import tools curl
#' @return path to temp file
#' @examples
#'    get_tempfile_from_url("https://ars.els-cdn.com/content/image/1-s2.0-S0092867419300571-mmc2.xlsx")
get_tempfile_from_url <- function(url) {
  temp <- tempfile(fileext = glue::glue(".{tools::file_ext(url)}"))
  curl::curl_download(url, temp)
  return(temp)
}

#' Print params list with knitr
#' @param pars Parameter list to print
#' @param omit Character vector with names of parameters that are not printed
#' @import knitr tibble
#' @export
#' @examples
#'   print_params(rmd_params, "counts")
print_params <- function(pars, omit = c()) {
  pars[omit] <- "not displayed"
  tibble::tibble(
    name = names(pars),
    setting = pars[names(pars)]
  ) %>%
    knitr::kable()
}

#' Set default theme for ggplot2
#' @export
set_ggplot_defaults <- function() {
  ggplot2::theme_set(
    ggplot2::theme_minimal(
      base_size = 12
    )
  )
  ggplot2::theme_update(
    axis.ticks = ggplot2::element_line(color = "grey92"),
    axis.ticks.length = ggplot2::unit(.5, "lines"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(color = "grey30"),
    plot.title = ggplot2::element_text(size = 18, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 12, color = "grey30"),
    plot.caption = ggplot2::element_text(size = 9, margin = ggplot2::margin(t = 15)),
    plot.title.position = 'plot'
  )
}