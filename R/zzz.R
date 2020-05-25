.onLoad <- function(libname, pkgname) {
  op <- options()
  op.rmyknife <- list(
    rmyknife.use_memoise = TRUE,
    rmyknife.memoise_path = file.path(
      dirname(tempdir()),
      "rmyknife_memoise"
    ),
    rmyknife.use_biomart_mirror = FALSE,
    rmyknife.biomart_mirror_host = "useast.ensembl.org",
    rmyknife.verbose = FALSE
  )
  toset <- !(names(op.rmyknife) %in% names(op))
  if (any(toset)) options(op.rmyknife[toset])

  if (isTRUE(getOption("rmyknife.use_memoise"))) {
    message("Caching biomaRt results to path ", op.rmyknife$rmyknife.memoise_path)
  }
  
  if (isTRUE(getOption("rmyknife.use_biomart_mirror"))) {
    message(paste0("Using biomart mirror ", getOption("biomart_mirror_host"), "\n\n Note that this only uses the latest Ensembl release version"))
  }
  if (isTRUE(getOption("rmyknife.verbose"))) {
    message("rmyknife verbose is set to TRUE. Use this only for troubleshooting purposes.")
  }

  invisible()
}