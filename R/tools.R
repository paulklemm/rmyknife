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

#' Volcano plot using data frame containing columns `log2FoldChange`, `padj` and `external_gene_name`
#' 
#' @param dat data frame containing columns `log2FoldChange`, `padj` and `external_gene_name`
#' @param min_padj Minimum p-value cutoff
#' @param min_log2fc Minimum log2fc cutoff
#' @param label_top_n Top-n genes to label
#' @param highlight Vector of `external_gene_name` to highlight
#' @export
plot_volcano <- function(
  dat,
  min_padj = 0.05,
  min_log2fc = 1,
  label_top_n = 10,
  highlight = c()
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
        " ratio)"
      )
    ) +
    ggplot2::theme_minimal() +
    ggrepel::geom_text_repel(
      data = . %>%
        # Highlighted genes are drawn separately, remove them here
        dplyr::filter(!(external_gene_name %in% highlight)) %>%
        dplyr::arrange(padj) %>%
        head(label_top_n),
      mapping = ggplot2::aes(label = external_gene_name),
      size = 3
    ) +
    # Highlight points
    ggplot2::geom_point(
      data = . %>%
        dplyr::filter(external_gene_name %in% highlight),
      alpha = 1,
      size = 0.8,
      color = "red"
    ) +
    # Highlight names
    ggrepel::geom_text_repel(
      data = . %>%
        dplyr::filter(external_gene_name %in% highlight),
      mapping = ggplot2::aes(label = external_gene_name),
      color = "red"
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

#' Get stream for gzipped content
#' @param url Url to call
#' @return stream for piping to read functions
#' @export
#' @examples
#' "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name&format=tsv&query=%28G-coupled%20coupled%20receptor%29%20AND%20%28model_organism%3A10090%29" %>%
#'    url() %>%
#'    gzcon() %>%
#'    readr::read_tsv()
get_gzipped_stream <- function(url) {
  url %>%
    url() %>%
    gzcon()
}

#' Make targets and load them into the global enrivonment
#' @param envir The environment to load the targets into
#' @param load_to_environment Whether to load the targets into the environment after building
#' @param workers Number of workers to use for make (requires crew package if > 1)
#' @return Number of outdated targets
#' @export
make <- function(
  envir = parent.frame(),
  load_to_environment = TRUE,
  workers = 10
) {

  library(magrittr)

   start_time <- Sys.time()  # Start timer
  
  # Set tar_options to use workers when workers > 1
  # See https://books.ropensci.org/targets/crew.html
  if (workers > 1) {
    targets::tar_option_set(
      controller = crew::crew_controller_local(
        workers = workers,
        seconds_idle = 3
      )
    )
  }

  # Check if make has been run before
  make_is_initialised <-
    exists(
      "make_initialised",
      envir = envir
    )
  
  # Make outdated targets
  outdated_targets <- targets::tar_outdated()
  outdated_targets_not_empty <- length(outdated_targets) > 0
  if (outdated_targets_not_empty) {
    paste0("Make outdated targets: ", paste(outdated_targets, collapse = ", ")) %>%
      message()
    targets::tar_make()
  }
  
  if (!make_is_initialised) {
    # Make is not initialised, we need to load all targets
    "Initialise make" %>%
      message()
    # Define global variable to indicate that make has been run before
    make_initialised <<- TRUE
    # Copy the entire tar_meta to the global environment
    tar_meta_local <<- targets::tar_meta()
    # Load the complete environment
    if (load_to_environment) {
      targets::tar_load(
        names = tidyselect::everything(),
        envir = envir
      )
    }
  } else if (load_to_environment) {
    # Get the current tar_meta table and compare the time column to the local one
    # If the time column is different, we need to load these targets
    tar_meta_global <- targets::tar_meta()
    # Get outdated targets in local environment
    outdated_targets_local <-
      dplyr::left_join(
        tar_meta_global %>%
          dplyr::filter(type == "stem") %>%
          dplyr::select(name, time),
        tar_meta_local %>%
          dplyr::filter(type == "stem") %>%
          dplyr::select(name, time),
        by = "name",
        suffix = c("_global", "_local")
      ) %>%
        dplyr::filter(
          # Target in global environment is newer than in local
          time_global > time_local |
          # We don't have the target in the local environment
          is.na(time_local)
        ) %>%
        dplyr::pull(name)

    # Message about loading targets
    paste0("Load outdated targets in local: ", paste(outdated_targets_local, collapse = ", ")) %>%
      message()

    # Load outdated targets
    targets::tar_load(
      tidyselect::all_of(outdated_targets_local),
      envir = envir
    )

    # Update the tar_meta_local table
    tar_meta_local <<- tar_meta_global
  }

  # Print timer
  end_time <- Sys.time()  # End timer
  total_time <- end_time - start_time
  total_seconds <- as.numeric(total_time, units = "secs")
  minutes <- floor(total_seconds / 60)
  seconds <- total_seconds %% 60
  message(sprintf("Total time taken for make: %d minutes and %.2f seconds", minutes, seconds))
  
  # Return number of outdated targets
  length(outdated_targets) %>%
    return()
}
