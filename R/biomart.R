#' Get ensembl host from ensembl_version without 'http://'
#' Example:
#'   get_ensembl_host_from_version(103)
#'   [1] "apr2018.archive.ensembl.org"
#'
#' @export
#' @import biomaRt magrittr stringr dplyr
#' @param ensembl_version Integer of ensembl_version
get_ensembl_host_from_version <- function(ensembl_version) {
  url <- biomaRt::listEnsemblArchives() %>%
    tibble::as_tibble() %>%
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
      # There is a duplicated `_id` column which we need to take care of before we can select
      tibble::as_tibble(.name_repair = "unique") %>%
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

#' Wrapper function for getting memoised function
#' Returns memoised function if option rmyknife.use_memoise is true
#' The memoised function will store it's value to getOption("rmyknife.memoise_path")
#' @import biomaRt memoise
#' @export
#' @examples
#'   get_memoised(biomaRt::getBM)
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
#' @param ensembl ensembl database connection object
#' @param verbose Print summary statistics to check for 1:1 or 1:N mappings
#'
#' @return Tibble of dat with attached biomaRt variables
#' @export
#' @examples
#'    tibble(
#'      gene_ids = c("Vamp8", "Mff", "Cers6", "Adcy3"),
#'      my_id = 1:4
#'    ) %>%
#'    attach_ensembl_gene_id_from_name(gene_name_var = "gene_ids")
#'
attach_ensembl_gene_id_from_name <- function(
  dat,
  gene_name_var = "Gene",
  ensembl = get_ensembl_dataset_from_version(103, species = "MUS"),
  verbose = TRUE) {
  # BiomaRt call
  dat_result <- get_memoised(biomaRt::getBM)(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "external_gene_name",
    values = dat %>%
      dplyr::select(gene_name_var) %>%
      dplyr::distinct() %>%
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
    paste0("Attaching Biomart gene information to input dataset (n = ", dat %>% nrow(), ", attached_n = ", dat_result %>% nrow(), ". NA-values in ensembl_gene_id = ", na_count,"). Species is ", ensembl@dataset, ".") %>%
      message()
  }
  dat_result %>%
    return()
}

#' Attach Biomart Gene identifier from entrez id
#' @param dat input data frame containing entrez id's
#' @param entrez_id_var Entrez id column in character type.
#' @param ensembl ensembl database connection object
#' @param verbose Print summary statistics to check for 1:1 or 1:N mappings
#'
#' @return Tibble of dat with attached biomaRt variables
#' @export
#' @examples
#'    tibble(
#'      entrez_gene_ids = c("4496", "4494", "4495", "1544"),
#'      my_id = 1:4
#'    ) %>%
#'    rmyknife::attach_ensembl_gene_id_from_entrez_id(entrez_id_var = "entrez_gene_ids")
#'
attach_ensembl_gene_id_from_entrez_id <- function(
  dat,
  entrez_id_var = "entrez_id",
  ensembl = get_ensembl_dataset_from_version(103, species = "MUS"),
  verbose = TRUE) {
  # BiomaRt call
  dat_result <- rmyknife::get_memoised(biomaRt::getBM)(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "entrezgene_id",
    values = dat %>%
      dplyr::select(entrez_id_var) %>%
      dplyr::distinct() %>%
      as.list(),
    mart = ensembl
  ) %>% tibble::as_tibble() %>%
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id))

  # Attach data to biomart output
  dat_result <- dat %>%
    dplyr::left_join(
      dat_result,
      by = setNames("entrezgene_id", entrez_id_var)
    )
  # Print output statistics when verbose is true
  if (verbose) {
    # Count the number of NA's in the column
    na_count <- dat_result$ensembl_gene_id %>%
      is.na() %>%
      sum()
    paste0("Attaching Biomart gene information to input dataset (n = ", dat %>% nrow(), ", attached_n = ", dat_result %>% nrow(), ". NA-values in ensembl_gene_id = ", na_count,"). Species is ", ensembl@dataset, ".") %>%
      message()
  }
  dat_result %>%
    return()
}

#' Attach Biomart variables based on either gene or transcript IDs
#' @param dat Input data frame containing either ensembl gene or transcript ids
#' @param ensembl_id_var Name of variable containing the gene/transcript ids
#' @param attributes Biomart attributes to retrieve
#' @param ensembl Biomart mart. Providing this overrides ensembl_version
#' @param ensembl_version Only required if ensembl not provided. Integer of ensembl version
#' @param verbose Print summary statistics to check for 1:1 or 1:N mappings
#'
#' @return Tibble of dat with attached biomaRt variables
#' @export
#' @examples
#'    test_data <- tibble(
#'      gene_ids = c("ENSMUSG00000102693", "ENSMUSG00000064842", "ENSMUSG00000102851", "ENSMUSG00000089699", "ENSMUSG00000103147", "ENSMUSG00000102348", "ENSMUSG00000102592", "ENSMUSG00000104238", "ENSMUSG00000102269", "ENSMUSG00000096126"),
#'      my_id = 1:10
#'    ) %>%
#'    attach_biomart(
#'      ensembl_id_var = "gene_ids",
#'      ensembl = get_ensembl_dataset_from_version(103)
#'
#'    test_data <- tibble(
#'      gene_ids = c("ENSMUSG00000102693", "ENSMUSG00000064842", "ENSMUSG00000102851", "ENSMUSG00000089699", "ENSMUSG00000103147", "ENSMUSG00000102348", "ENSMUSG00000102592", "ENSMUSG00000104238", "ENSMUSG00000102269", "ENSMUSG00000096126"),
#'      my_id = 1:10
#'    ) %>%
#'    attach_biomart(ensembl_id_var = "gene_ids", ensembl_version = 103)
#'
#'    test_data_transcripts <- tibble(
#'      transcripts = c("ENSMUST00000082423", "ENSMUST00000082422", "ENSMUST00000082421", "ENSMUST00000082420", "ENSMUST00000082419", "ENSMUST00000082418", "ENSMUST00000082417", "ENSMUST00000082416", "ENSMUST00000082415", "ENSMUST00000082414"),
#'      ids = 1:10
#'    )
#'    attach_biomart(test_data_transcripts, ensembl_id_var = "transcripts", ensembl_version = 103)
attach_biomart <- function(
    dat,
    ensembl_id_var = "ensembl_gene_id",
    attributes = c("description", "gene_biotype", "external_gene_name"),
    ensembl = NULL,
    ensembl_version = 103,
    verbose = TRUE
  ) {
  # Check if we have genes or transcripts as input based on the first element
  type <- dat %>%
    dplyr::select(ensembl_id_var) %>%
    head(1) %>%
    stringr::str_match(pattern = "ENS[a-zA-Z]{3}(\\w)") %>%
    # Identifier is the second entry
    .[2]
  # If type is NA, maybe we have human genes that look like ENSG00000140505
  if (is.na(type)) {
    type <- dat %>%
      dplyr::select(ensembl_id_var) %>%
      head(1) %>%
      stringr::str_match(pattern = "ENS(\\w)") %>%
      # Identifier is the second entry
      .[2]
  }

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
  first_gene_id <- dat %>%
    dplyr::select(ensembl_id_var) %>%
    head(1)
  
  species_id <- first_gene_id %>%
    stringr::str_match(pattern = "ENS([a-zA-Z]{3})") %>%
    .[2]
  # Test for human gene type. names look like ENSG00000140505
  if (is.na(species_id)) {
    if (first_gene_id %>% stringr::str_detect(pattern = "^ENS[GT](\\d*)$")) {
      species_id <- "HUM"
    }
  }
  verbose_mart_text <- "Using user-provided ensembl mart."
  # Get mart based on information deduced from IDs
  if (is.null(ensembl)) {
    verbose_mart_text <- glue::glue("Using mart version {ensembl_version} for species {species_id}.")
    # Get Ensembl dataset
    ensembl <- get_ensembl_dataset_from_version(
      ensembl_version = ensembl_version,
      species = species_id
    )
  }
  # BiomaRt call
  dat_result <- get_memoised(biomaRt::getBM)(
    attributes = attributes,
    filters = filter_type,
    values = dat %>%
      dplyr::select(ensembl_id_var) %>%
      dplyr::distinct() %>%
      as.list(),
    mart = ensembl
  ) %>% tibble::as_tibble()
  # Force the keys join keys to be of same type. This is rarely required
  # For example it can happen with retired genes, that no result is reported with getBM
  # leading to a ensembl_gene_id column of type lgl (e.g. gene ENSMUSG00000089672)
  mode(dat_result[filter_type][[1]]) <- mode(dat[ensembl_id_var][[1]])
  # Attach data to biomart output
  dat_result <- dat %>%
    dplyr::left_join(dat_result, by = setNames(filter_type, ensembl_id_var))
  # Print output statistics when verbose is true
  if (verbose) {
    paste0(verbose_mart_text, " Attaching Biomart gene information to input dataset (n = ", dat %>% nrow(), ", attached_n = ", dat_result %>% nrow(), "). Species is ", ensembl@dataset, ".") %>%
      message()
  }
  dat_result %>%
    return()
}

#' Get Ensembl dataset with default parameter
#' @param ensembl_version Ensembl version you want to get data for
#' @param species Species you want to extract ("MUS" or "HUM")
#' @import biomaRt magrittr
#' @export
#'
#' @examples
#'    ensembl <- get_ensembl_dataset_from_version(103, "MUS")
get_ensembl_dataset_from_version <- function(
    ensembl_version = 103,
    species = "MUS"
  ) {
  # Setup the dataset based on species
  if (species == "MUS") {
    ensembl_dataset <- "mmusculus_gene_ensembl"
  } else if (species == "HUM") {
    ensembl_dataset <- "hsapiens_gene_ensembl"
  } else {
    paste0("Species ", species, " not supported") %>%
      stop()
  }

  # Check if we should use the backup version or not
  if (isTRUE(getOption("rmyknife.use_biomart_mirror"))) {
    debug_message("Using biomart mirror host in `get_ensembl_dataset_from_version`")
    ensembl <- biomaRt::useMart(
      host = getOption("rmyknife.biomart_mirror_host"),
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = ensembl_dataset
    )
  # Use the standard version
  } else {
    debug_message("Using standard biomart host in `get_ensembl_dataset_from_version`")
    ensembl <- biomaRt::useMart(
      host = rmyknife::get_ensembl_host_from_version(ensembl_version),
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = ensembl_dataset
    )
  }

  return(ensembl)
}

#' Get genes associated with GO-term based on BiomaRt
#' This will get the genes of the selected GO-term and
#' *all* child-terms.
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
    tibble::as_tibble() %>%
    return()
}

#' Get genes associated with GO-term based on BiomaRt
#' @param go_accession ID of GO term
#' @param ensembl Biomart connection
#' @param verbose Print summary statistic of the query
#' @export
#' @import dplyr biomaRt magrittr
#' @return Tibble with all genes associated with GO-term
#'
#' @examples
#'    get_genes_of_goterm(go_accession = "GO:0006811", ensembl = rmyknife::get_ensembl_dataset_from_version(103, species = "MUS"))
get_genes_of_goterm <- function(
  go_accession,
  ensembl,
  verbose = TRUE
) {
  # go_accession <- "GO:0006811"; ensembl <- rmyknife::get_ensembl_dataset_from_version(103, "mmusculus_gene_ensembl") ; verbose <- TRUE
  go_terms <- get_genes_of_goterm_helper(go_accession, ensembl)
  if (verbose) {
    go_name <- get_goterm_name_from_id(go_accession, ensembl)
    # Print out verbose message
    paste0("Get genes of GO term ", go_accession, " (", go_name, "): ", go_terms %>% nrow(), " genes found") %>%
      message()
  }
  go_terms %>%
    return()
}

#' Get genes associated with GO-term as vector based on GO.db
#' @param go_accession ID of GO term
#' @param species Define species, either "HUM" or "MUS"
#' @param verbose Print summary statistic of the query
#' @param memoised Use memoised function (search for local caches of given parameters)
#' @import GO.db tibble AnnotationDbi org.Hs.eg.db org.Mm.eg.db magrittr dplyr
#' @export
#' @examples
#'   get_genes_of_goterm_godb("GO:1900746")
get_genes_of_goterm_godb <- function(
  go_accession,
  species = "MUS",
  verbose = TRUE,
  memoised = TRUE
) {
  # Get memoised version of helper function if required
  if (memoised) {
    get_genes_of_goterm_godb_func <- get_memoised(get_genes_of_goterm_godb_helper)
  } else {
    get_genes_of_goterm_godb_func <- get_genes_of_goterm_godb_helper
  }

  get_genes_of_goterm_godb_func(
    go_accession = go_accession,
    species = species,
    verbose = verbose
  ) %>%
    return()
}

#' Get genes associated with GO-term as vector based on GO.db
#' @param go_accession See get_genes_of_goterm_godb
#' @param species See get_genes_of_goterm_godb
#' @param verbose See get_genes_of_goterm_godb
#' @import GO.db tibble AnnotationDbi org.Hs.eg.db org.Mm.eg.db magrittr dplyr
#' @examples
#'   get_genes_of_goterm_godb_helper("GO:1900746")
get_genes_of_goterm_godb_helper <- function(
  go_accession,
  species = "MUS",
  verbose = TRUE
) {
  # Code adapted from https://davetang.org/muse/2011/05/20/extract-gene-names-according-to-go-terms/
  # Check for species
  if (species == "MUS") {
    go2allegs <- org.Mm.eg.db::org.Mm.egGO2ALLEGS
    symbol <- org.Mm.eg.db::org.Mm.egSYMBOL
  } else if (species == "HUM") {
    go2allegs <- org.Hs.eg.db::org.Hs.egGO2ALLEGS
    symbol <- org.Hs.eg.db::org.Hs.egSYMBOL
  } else {
    paste0("Species ", species, " not supported") %>%
      stop()
  }
  # Get GOterm genes
  goterm_genes <- AnnotationDbi::get(go_accession, go2allegs) %>%
    AnnotationDbi::mget(symbol) %>%
    unlist() %>%
    tibble::enframe(name = NULL) %>%
    dplyr::distinct() %>%
    dplyr::pull()
  # Verbose output
  if (verbose) {
    go_name <- get_goterm_name_from_id_godb(go_accession)
    # Print out verbose message
    paste0("Get genes of GO term ", go_accession, " (", go_name, "): ", goterm_genes %>% length(), " genes found") %>%
      message()
  }
  # Return goterm genes
  return(goterm_genes)
}

#' DEPRECATED. Get GO name using GO.db
#' @import biomaRt dplyr GO.db magrittr AnnotationDbi
#' @export
#' @param go_accession GO term id, e.g. "GO:0032680"
#' @return GO-term name
#'
#' @examples
#'    get_goterm_name_from_id_godb(go_accession = "GO:0032680")
get_goterm_name_from_id_godb <- function(go_accession) {
  # Convert GOterm DB to a tibble we can filter on
  goterm_name <-
    AnnotationDbi::Term(GO.db::GOTERM) %>%
    as.list() %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(cols = tidyr::everything(), names_to = "accession", values_to = "name") %>%
    # Filter for term
    dplyr::filter(accession == go_accession)
  if (goterm_name %>% nrow() != 1) {
    paste0("Could not find GO-term ", go_accession, " in GO.db") %>%
      stop()
  }
  goterm_name %>%
    .$name %>%
    return()
}

#' Get GO name using biomaRt
#' @import biomaRt tibble magrittr
#' @export
#' @param go_accession GO term id, e.g. "GO:0032680"
#' @param ensembl ensembl biomaRt object
#' @return GO-term name
#'
#' @examples
#'    get_goterm_name_from_id(go_accession = "GO:0032680")
get_goterm_name_from_id <- function(
  go_accession,
  ensembl = get_ensembl_dataset_from_version(103)
) {
  # BiomaRt call
  dat_result <-
    get_memoised(biomaRt::getBM)(
      attributes = c("go_id", "name_1006"),
      mart = ensembl
    ) %>%
    tibble::as_tibble() %>%
    dplyr::filter(go_id == go_accession)
  if (dat_result %>% nrow() != 1) {
    paste0("Could not find GO-term ", go_accession, " in provided ensembl mart") %>%
      stop()
  }
  dat_result %>%
    .$name_1006 %>%
    return()
}

#' Get promotor sequence upstream of a gene flank
#' @param ensembl_gene_ids genes to get sequence from
#' @param ensembl biomaRt connection
#' @param upstream_bases How many bases upstream of the gene start to return
#' @return tibble with columns gene_flank and ensembl_gene_id
#' @export
#' @import tibble
#' @examples
#'    get_promotor_sequence(
#'      ensembl_gene_ids = c("ENSMUSG00000102693", "ENSMUSG00000064842", "ENSMUSG00000102851"),
#'      ensembl = rmyknife::get_ensembl_dataset_from_version(103)
#'    )
#'    get_promotor_sequence(
#'      ensembl_gene_ids = c("ENSG00000140505", "ENSG00000205358", "ENSG00000125144", "ENSG00000198417", "ENSG00000205364", "ENSG00000169715", "ENSG00000187193", "ENSG00000125148"),
#'      ensembl = rmyknife::get_ensembl_dataset_from_version(103, species = "HUM")
#'    )
get_promotor_sequence <- function(
  ensembl_gene_ids,
  ensembl = get_ensembl_dataset_from_version(103),
  upstream_bases = 1000
) {
  get_memoised(biomaRt::getSequence)(
    id = ensembl_gene_ids,
    mart = ensembl,
    type = "ensembl_gene_id",
    # We probably need "coding_gene_flank"
    # Source: https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#given-a-set-of-entrezgene-identifiers-retrieve-100bp-upstream-promoter-sequences
    seqType = "coding_gene_flank",
    upstream = upstream_bases
  ) %>%
    tibble::as_tibble() %>%
    return()
}

#' Wrapper function for get_promotor_sequence. Attaches the result to a data frame
#' @param dat tibble with ensembl_gene_id variable
#' @param ensembl_id_var ensembl id variable. Default is "ensembl_gene_id"
#' @param ensembl ... from get_promotor_sequence
#' @param upstream_bases ... from get_promotor_sequence
#' @return original tibble with column gene_flank added
#' @export
#' @import tibble
#' @examples
#'   tibble::tibble(EnsemblIDs = c("ENSMUSG00000102693", "ENSMUSG00000064842", "ENSMUSG00000102851")) %>%
#'     get_promotor_sequence_tibble(
#'       ensembl_id_var = "EnsemblIDs",
#'       ensembl = get_ensembl_dataset_from_version(103, species = "MUS")
#'     )
get_promotor_sequence_tibble <- function(
  dat,
  ensembl_id_var = "ensembl_gene_id",
  ensembl = get_ensembl_dataset_from_version(103),
  upstream_bases = 1000
) {
  dplyr::left_join(
    dat,
    get_promotor_sequence(
      ensembl_gene_ids = dplyr::select(dat, ensembl_id_var) %>% dplyr::pull(),
      upstream_bases = upstream_bases,
      ensembl = ensembl
    ),
    by = setNames("ensembl_gene_id", ensembl_id_var)
  ) %>%
    return()
}

#' Get orgdb for species as a string
#' @param species Either MUS or HUM
#' @return OrgDB database string
#' @export
#' @import magrittr
get_orgdb_for_species <- function(species = "MUS") {
  if (species == "MUS") {
    return("org.Mm.eg.db")
  } else if (species == "HUM") {
    return("org.Hs.eg.db")
  } else {
    paste0("Species ", species, " not supported") %>%
      stop()
  }
}

#' Gets Entrez IDs and add it to dat
#'
#' @export
#' @import magrittr clusterProfiler org.Mm.eg.db dplyr tibble
#' @param dat Data frame with ensembl identifier
#' @param ensembl_id_name Name of column containing the ensembl identifier
#' @param keep_only_rows_with_entrez Only keep rows for which entrez IDs could be found
#' @param drop_duplicates Often there is a n:1 relationship between Entrez-IDs and Ensembl-IDs. If this value is true, only keep the first hit
#' @param species Either MUS or HUM
#' @examples
#'   tibble::tibble(
#'     ensembl_gene_id = c("ENSMUSG00000000001","ENSMUSG00000000003","ENSMUSG00000000028","ENSMUSG00000000031","ENSMUSG00000000037","ENSMUSG00000000049","ENSMUSG00000000056","ENSMUSG00000000058","ENSMUSG00000000078","ENSMUSG00000000085"),
#'   ) %>%
#'     rmyknife::ensembl_to_entrez()
ensembl_to_entrez <- function(
  dat,
  ensembl_id_name = "ensembl_gene_id",
  keep_only_rows_with_entrez = TRUE,
  drop_duplicates = TRUE,
  species = "MUS"
) {
  org_db <- get_orgdb_for_species(species)
  # Get the Entrez IDs
  ens_to_ent <- dat[ensembl_id_name][[1]] %>%
    clusterProfiler::bitr(
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = org_db
    ) %>%
    dplyr::rename(EntrezID = ENTREZID)
  # Drop duplicates if required
  if (drop_duplicates) {
    ens_to_ent <- ens_to_ent %>% dplyr::filter(!duplicated(ENSEMBL))
  }
  if (keep_only_rows_with_entrez) {
    dat %>%
      dplyr::right_join(
        ens_to_ent,
        by = setNames("ENSEMBL", ensembl_id_name)
      ) %>%
      tibble::as_tibble() %>%
      return()
  } else {
    dat %>%
      dplyr::full_join(
        ens_to_ent,
        by = setNames("ENSEMBL", ensembl_id_name)
      ) %>%
      tibble::as_tibble() %>%
      return()
  }
}

#' Add Gene Symbol from EntrezID
#'
#' @export
#' @import clusterProfiler dplyr magrittr
#' @param dat data frame with EnzrezID
#' @param species Either "MUS" or "HUM"
#' @examples
#'   tibble::tibble(
#'     ensembl_gene_id = c("ENSMUSG00000000001","ENSMUSG00000000003","ENSMUSG00000000028","ENSMUSG00000000031","ENSMUSG00000000037","ENSMUSG00000000049","ENSMUSG00000000056","ENSMUSG00000000058","ENSMUSG00000000078","ENSMUSG00000000085"),
#'   ) %>%
#'     rmyknife::ensembl_to_entrez() %>%
#'     rmyknife::attach_gene_symbol_from_entrez()
attach_gene_symbol_from_entrez <- function(
  dat,
  species = "MUS"
) {
  if (!('EntrezID' %in% colnames(dat))) {
    stop('EntrezID not available in data frame')
  }
  org_db <- get_orgdb_for_species(species)
  dat["EntrezID"][[1]] %>%
    clusterProfiler::bitr(
      fromType = "ENTREZID",
      toType = c("SYMBOL"),
      OrgDb = org_db
    ) %>%
    dplyr::rename(EntrezID = ENTREZID, Symbol = SYMBOL) %>%
    dplyr::right_join(dat) %>%
    as_tibble() %>%
    return()
}
