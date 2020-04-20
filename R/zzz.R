.onLoad <- function(libname, pkgname) {
  op <- options()
  op.rmyknife <- list(
    rmyknife.use_memoise = TRUE,
    rmyknife.memoise_path = file.path(
      dirname(tempdir()),
      "rmyknife_memoise"
    )
  )
  toset <- !(names(op.rmyknife) %in% names(op))
  if (any(toset)) options(op.rmyknife[toset])

  if (isTRUE(getOption("rmyknife.use_memoise"))) {
    message("Caching biomaRt results to path ", op.rmyknife$rmyknife.memoise_path)
  }

  invisible()
}