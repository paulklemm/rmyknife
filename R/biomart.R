#' Get ensembl host from ensembl_version without 'http://'
#' Example:
#'   get_ensembl_host_from_version(92)
#'   [1] "apr2018.archive.ensembl.org"
#'
#' @export
#' @import biomaRt magrittr stringr dplyr
#' @param ensembl_version Integer of ensembl_version
get_ensembl_host_from_version <- function(ensembl_version) {
  url <- biomaRt::listEnsemblArchives() %>%
    tibble::as.tibble() %>%
    dplyr::filter(name == paste0('Ensembl ', ensembl_version)) %>%
    .$url %>%
    stringr::str_remove('http://')
  if (identical(url, character(0))) {
    paste0('No host found for Ensembl version ', ensembl_version) %>%
      stop()
  } else {
    url %>% return()
  }
}

#' Get gene names from synonyms
#' @export 
#' @import tibble dplyr DBI org.Mm.eg.db magrittr
#' @param gene_name Vector of gene names
#' @param species Define species. Default is "MUS"
#' @param keep_missing If we cannot identify a gene, still keep the original name?
get_gene_name_from_synonym <- function(
    gene_name,
    species = "MUS",
    keep_missing = TRUE
  ) {
  # Setup the dataset based on species
  if (species == "MUS") {
    # Get a tibble of all gene name synonyms
    synonyms <- DBI::dbGetQuery(
      org.Mm.eg.db::org.Mm.eg_dbconn(),
      'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
    ) %>%
      tibble::as_tibble() %>%
      dplyr::select(alias_symbol, symbol) %>%
      dplyr::mutate(alias_symbol = tolower(alias_symbol))
  } else {
    paste0("Species ", species, " not supported") %>%
      stop()
  }

  # Convert vector to tibble
  gene_name_return <- tolower(gene_name) %>%
    tibble::tibble(gene_name = .)
  
  get_synonym <- function(gene) {
    synonym <- synonyms %>%
      dplyr::filter(alias_symbol == gene)
    
    # Handle multiple rows. Currently we just choose the top one
    if (nrow(synonym) > 1) {
      synonym %<>%
        # Oddly enough, a synonym can map to multiple genes. Therefore we will only keep the first
        # gene that end's up in the list to keep a n:n relationship of the function
        head(1)
    } else if (nrow(synonym) == 0) {
      # Handle no return
      if (keep_missing) {
        return(gene)
      } else {
        return(NA)
      }
    }
    # Return result as a single value
    synonym %>%
      dplyr::select(symbol) %>%
      dplyr::pull() %>%
      return()
  }

  # Filter the synonyms list for the ones we're interested in
  gene_name_return <- gene_name_return %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gene_name = get_synonym(gene_name))
  
  paste0("Returning ", nrow(gene_name_return), " gene names from ", nrow(gene_name_return), " input gene name(s) (", is.na(gene_name_return$gene_name) %>% sum(), " NA value(s))") %>%
    message()
  
  gene_name_return %>%
    dplyr::select(gene_name) %>%
    dplyr::pull() %>%
    return()
}

#' Wrapper function for getting Memoised function
#' Returns memoised function if option rmyknife.use_memoise is true
#' @import biomaRt memoise
get_memoised <- function(func) {
  mem_func <- func
  if (isTRUE(getOption("rmyknife.use_memoise"))) {
    mem_func <- memoise::memoise(
      func,
      cache = memoise::cache_filesystem(getOption("rmyknife.memoise_path"))
    )
  }
  return(mem_func)
}

#' Attach Biomart Gene identifier from gene name
#' @param dat input data frame containing Gene names
#' @param gene_name_var Gene name column
#' @param attributes biomart attributes to retrieve
#' @param ensembl_version integer of ensembl version
#' @param species Species identifier as string
#' @param verbose Print summary statistics to check for 1:1 or 1:N mappings
#'
#' @return Tibble of dat with attached biomaRt variables
#' @export
#' @examples
#'    tibble(
#'      gene_ids = c("Vamp8", "Mff", "Cers6", "Adcy3"),
#'      my_id = 1:4
#'    ) %>%
#'    attach_ensembl_gene_id_from_name(gene_name_var = "gene_ids", ensembl_version = 96, species = "MUS")
#'
attach_ensembl_gene_id_from_name <- function(
  dat,
  gene_name_var = "Gene",
  ensembl_version,
  species = "MUS",
  verbose = TRUE) {
  # Setup the dataset based on species
  if (species == "MUS") {
    ensembl_dataset <- "mmusculus_gene_ensembl"
  } else {
    paste0("Species ", species, " not supported") %>%
      stop()
  }
  ensembl <- rmyknife::get_ensembl_dataset_from_version(ensembl_version, ensembl_dataset)
  # BiomaRt call
  dat_result <- get_memoised(biomaRt::getBM)(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "external_gene_name",
    values = dat %>%
      dplyr::select_(gene_name_var) %>%
      dplyr::distinct_(gene_name_var) %>%
      as.list(),
    mart = ensembl
  ) %>% tibble::as_tibble()
  # Attach data to biomart output
  dat_result <- dat %>%
    dplyr::left_join(dat_result, by = setNames("external_gene_name", gene_name_var))
  # Print output statistics when verbose is true
  if (verbose) {
    # Count the number of NA's in the column
    na_count <- dat_result$ensembl_gene_id %>%
      is.na() %>%
      sum()
    paste0("Attaching Biomart gene information to input dataset (n = ", dat %>% nrow(), ", attached_n = ", dat_result %>% nrow(), ". NA-values in ensembl_gene_id = ", na_count,"). Species is ", ensembl_dataset, ".") %>%
      message()
  }
  dat_result %>%
    return()
}

#' Attach Biomart variables based on either gene or transcript IDs
#' @param dat input data frame containing either ensembl gene or transcript ids
#' @param ensembl_id_var name of variable containing the gene/transcript ids
#' @param attributes biomart attributes to retrieve
#' @param ensembl_version integer of ensembl version
#' @param verbose Print summary statistics to check for 1:1 or 1:N mappings
#'
#' @return Tibble of dat with attached biomaRt variables
#' @export
#' @examples
#'    test_data <- tibble(
#'      gene_ids = c("ENSMUSG00000102693", "ENSMUSG00000064842", "ENSMUSG00000102851", "ENSMUSG00000089699", "ENSMUSG00000103147", "ENSMUSG00000102348", "ENSMUSG00000102592", "ENSMUSG00000104238", "ENSMUSG00000102269", "ENSMUSG00000096126"),
#'      my_id = 1:10
#'    )
#'    attach_biomart(test_data, ensembl_id_var = "gene_ids", ensembl_version = 94)
#'
#'    test_data_transcripts <- tibble(
#'      transcripts = c("ENSMUST00000082423", "ENSMUST00000082422", "ENSMUST00000082421", "ENSMUST00000082420", "ENSMUST00000082419", "ENSMUST00000082418", "ENSMUST00000082417", "ENSMUST00000082416", "ENSMUST00000082415", "ENSMUST00000082414"),
#'      ids = 1:10
#'    )
#'    attach_biomart(test_data_transcripts, ensembl_id_var = "transcripts", ensembl_version = 92)
attach_biomart <- function(
    dat,
    ensembl_id_var = "ensembl_gene_id",
    attributes = c("description", "gene_biotype", "external_gene_name"),
    ensembl_version = 94,
    verbose = TRUE
  ) {
  # Check if we have genes or transcripts as input based on the first element
  type <- dat %>% dplyr::select_(ensembl_id_var) %>%
    head(1) %>%
    stringr::str_match(pattern = "ENS[a-zA-Z]{3}(\\w)") %>%
    # Identifier is the second entry
    .[2]
  verbose_id_text <- ""
  filter_type <- ""
  # Check gene type
  if (type == "G") {
    verbose_id_text <- "Identifier type is Genes"
    attributes %<>% c("ensembl_gene_id")
    filter_type <- "ensembl_gene_id"
  } else if (type == "T") {
    verbose_id_text <- "Identifier type is Transcripts"
    attributes %<>% c("ensembl_transcript_id")
    filter_type <- "ensembl_transcript_id"
    # type is neither G or T
  } else {
    paste0("Cannot identify type (gene or transcript) from gene identifier") %>%
      stop()
  }
  species_id <- dat %>% dplyr::select_(ensembl_id_var) %>%
    head(1) %>%
    stringr::str_match(pattern = "ENS([a-zA-Z]{3})") %>%
    .[2]
  # Setup the dataset based on species
  if (species_id == "MUS") {
    ensembl_dataset <- "mmusculus_gene_ensembl"
  } else {
    paste0("Species ", species_id, " not supported") %>%
      stop()
  }
  # Get Ensembl dataset
  ensembl <- rmyknife::get_ensembl_dataset_from_version(ensembl_version, ensembl_dataset)
  # BiomaRt call
  dat_result <- get_memoised(biomaRt::getBM)(
    attributes = attributes,
    filters = filter_type,
    values = dat %>%
      dplyr::select_(ensembl_id_var) %>%
      dplyr::distinct_(ensembl_id_var) %>%
      as.list(),
    mart = ensembl
  ) %>% tibble::as.tibble()
  # Attach data to biomart output
  dat_result <- dat %>%
    dplyr::left_join(dat_result, by = setNames(filter_type, ensembl_id_var))
  # Print output statistics when verbose is true
  if (verbose) {
    paste0("Attaching Biomart gene information to input dataset (n = ", dat %>% nrow(), ", attached_n = ", dat_result %>% nrow(), "). Species is ", ensembl_dataset, ".") %>%
      message()
  }
  dat_result %>%
    return()
}

#' Get Ensembl dataset with default parameter
#' @param ensembl_version Ensembl version you want to get data for
#' @param ensembl_dataset Species you want to extract
#' @import biomaRt magrittr
#' @export
#'
#' @examples
#'    ensembl <- get_ensembl_dataset_from_version(94, "mmusculus_gene_ensembl")
get_ensembl_dataset_from_version <- function(
    ensembl_version = 94,
    ensembl_dataset = "mmusculus_gene_ensembl"
  ) {
  get_memoised(biomaRt::useMart)(
    host = rmyknife::get_ensembl_host_from_version(ensembl_version),
    biomart = 'ENSEMBL_MART_ENSEMBL',
    dataset = ensembl_dataset
  ) %>%
    return()
}

#' Get genes associated with GO-term based on BiomaRt
#' @param go_accession ID of GO term
#' @param ensembl Biomart connection
#' @import biomaRt tibble
get_genes_of_goterm_helper <- function(
  go_accession,
  ensembl
) {
  get_memoised(biomaRt::getBM)(
    attributes = c("ensembl_gene_id"),
    filters = "go_parent_term",
    values = c(go_accession),
    mart = ensembl
  ) %>%
    as_tibble() %>%
    return()
}

#' Get genes associated with GO-term based on BiomaRt
#' @param go_accession ID of GO term
#' @param ensembl Biomart connection
#' @export
#' @import dplyr biomaRt magrittr
#' @return Tibble with all genes associated with GO-term
#'
#' @examples
#'    get_genes_of_goterm(go_accession = "GO:0006811", ensembl = rmyknife::get_ensembl_dataset_from_version(94, "mmusculus_gene_ensembl"))
get_genes_of_goterm <- function(
  go_accession,
  ensembl,
  verbose = TRUE
) {
  # go_accession <- "GO:0006811"; ensembl <- rmyknife::get_ensembl_dataset_from_version(94, "mmusculus_gene_ensembl") ; verbose <- TRUE
  go_terms <- get_genes_of_goterm_helper(go_accession, ensembl)
  if (verbose) {
    go_name <- get_goterm_name_from_id(go_accession)
    # Print out verbose message
    paste0("Get genes of GO term ", go_accession, " (", go_name, "): ", go_terms %>% nrow(), " genes found") %>%
      message()
  }
  go_terms %>%
    return()
}

#' Get GO name.
#' @import biomaRt dplyr GO.db magrittr AnnotationDbi
#' @export
#' @param go_accession GO term id, e.g. "GO:0032680"
#'
#' @examples
#'    get_goterm_name_from_id(go_accession = "GO:0032680")
get_goterm_name_from_id <- function(go_accession) {
  # Convert GOterm DB to a tibble we can filter on
  goterm_name <- AnnotationDbi::Term(GO.db::GOTERM) %>%
    as.list() %>%
    tibble::as_tibble() %>%
    tidyr::gather(key = "accession", value = "name") %>%
  # Filter for term
  dplyr::filter(accession == go_accession)
  if (goterm_name %>% nrow() != 1) {
    paste0("Could not find GO-term ", go_accession, " in provided ensembl data") %>%
      stop()
  }
  goterm_name %>%
    .$name %>%
    return()
}

# #' @param dat
# #' @param ensembl_id_var
# #' @param attributes
# #' @param ensembl_version
# #' @param verbose
# #'
# #' @return Tibble of dat with attached biomaRt variables
# #' @export
# #' @examples
# #'    test_data <- tibble(
# #'      gene_ids = c("ENSMUSG00000000001","ENSMUSG00000000003","ENSMUSG00000000028","ENSMUSG00000000031","ENSMUSG00000000037","ENSMUSG00000000049","ENSMUSG00000000056","ENSMUSG00000000058","ENSMUSG00000000078","ENSMUSG00000000085"),
# #'      my_id = 1:10
# #'    )
# #'    get_goterm_genes_from_biomart(test_data, ensembl_id_var = "gene_ids", ensembl_version = 94)
# #'
# #'    test_data_transcripts <- tibble(
# #'      transcripts = c("ENSMUST00000082423", "ENSMUST00000082422", "ENSMUST00000082421", "ENSMUST00000082420", "ENSMUST00000082419", "ENSMUST00000082418", "ENSMUST00000082417", "ENSMUST00000082416", "ENSMUST00000082415", "ENSMUST00000082414"),
# #'      ids = 1:10
# #'    )
# #'    get_goterm_genes_from_biomart(test_data_transcripts, ensembl_id_var = "transcripts", ensembl_version = 92)
# attach_goterms_from_biomart <- function(dat, ensembl_id_var, ensembl_version, verbose) {
#   attributes <- c("name_1006", "definition_1006")
#   go_terms <-
#     attach_biomart(dat = dat, ensembl_id_var = ensembl_id_var, attributes = attributes, ensembl_version = ensembl_version, verbose = verbose) %>%
#     return()
#     # # Select only the GO term attributes
#     # subset(select = attributes) %>%
#     # # dplyr::distinct() %>%
#     # return()
# }
