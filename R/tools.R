#' Read in GTF file
#'
#' @param path Path to GTF file
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
#' @export
#' @return Tibble of cufflinks data
read_cufflinks <- function(path) {
  path %>%
    readr::read_tsv() %>%
    return()
}

#' Print standard DT::datatable with extensions for exporting the data
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

#' Volcano plot using data frame containing columns `log2FoldChange` and `padj``
#' 
#' @param dat data frame containing columns `log2FoldChange` and `padj`
#' @param min_padj Minimum p-value cutoff
#' @param min_log2fc Minimum log2fc cutoff
#' @param label_top_n Top-n genes to label
#' @export
plot_volcano <- function(
  dat,
  min_padj = 0.05,
  min_log2fc = 1,
  label_top_n = 10
) {
  # Attach significance
  dat <-
    dat %>%
    dplyr::mutate(
      significance =
        (
          (padj <= min_padj & log2FoldChange >= min_log2fc) |
          (padj <= min_padj & log2FoldChange <= -min_log2fc)
        ) %>%
          ifelse(., "significant", "not significant")
    )
  
  volcanoplot <-
    dat %>%
    # Remove entries that cannot be drawn
    dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = log2FoldChange,
        y = -log10(padj),
        color = significance
      )
    ) +
    ggplot2::geom_point(
      alpha = 0.3,
      size = 0.5
    ) +
    ggplot2::xlab(expression(log[2](fc))) +
    ggplot2::ylab(expression(-log[10](adjusted ~ p ~ value))) +
    ggplot2::labs(colour =
      paste0(
        "Significance\npadj = ",
        min_padj,
        "\nl2fc = ± ",
        min_log2fc,
        "(",
        round(2^min_log2fc, digits = 2),
        "-fold)"
      )
    ) +
    ggplot2::theme_minimal() +
    ggrepel::geom_text_repel(
      data = . %>%
        dplyr::arrange(padj) %>%
        head(label_top_n),
      mapping = ggplot2::aes(label = external_gene_name),
      size = 3
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(min_padj),
      linetype = "dotted"
    ) +
    ggplot2::geom_vline(
      xintercept = -min_log2fc,
      linetype = "dotted"
    ) +
    ggplot2::geom_vline(
      xintercept = min_log2fc,
      linetype = "dotted"
    )
  
  # Apply proper coloring
  significance_summary <-
    dat %>%
    dplyr::distinct(significance)
  if(nrow(significance_summary) > 1) {
    volcanoplot <-
      volcanoplot +
      ggplot2::scale_color_manual(values = c("grey", "blue"))
  } else if(dplyr::pull(significance_summary) == "significant") {
    volcanoplot <-
      volcanoplot +
      ggplot2::scale_color_manual(values = "blue")
  } else {
    volcanoplot <-
      volcanoplot +
      ggplot2::scale_color_manual(values = "grey")
  }
  return(volcanoplot)
}
