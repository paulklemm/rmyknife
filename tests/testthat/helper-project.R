# Fixtures for the reproducible environment tests.
#
# A fake project is a temporary directory carrying the files project_verify()
# reads, plus a small stand-in for the singularity image so checksums and sizes
# are real without moving gigabytes around.

fake_image <- function(dir, name = "mytidyverse-4.6.1-1.simg", content = "not really an image") {
  path <- file.path(dir, name)
  writeLines(content, path)
  normalizePath(path)
}

fake_renv_lock <- function(
  path,
  r_version = "4.6.1",
  snapshot_date = "2026-07-24",
  bioc_version = "3.23",
  packages = list(dplyr = list(Package = "dplyr", Version = "1.2.1", Source = "Repository", Repository = "CRAN"))
) {
  lock <- list(
    R = list(
      Version = r_version,
      Repositories = list(
        list(Name = "CRAN", URL = paste0("https://p3m.dev/cran/", snapshot_date)),
        list(Name = "BioCsoft", URL = paste0("https://bioconductor.org/packages/", bioc_version, "/bioc"))
      )
    ),
    Bioconductor = list(Version = bioc_version),
    Packages = packages
  )
  writeLines(jsonlite::toJSON(lock, auto_unbox = TRUE, pretty = TRUE), file.path(path, "renv.lock"))
}

fake_project <- function(
  r_version = "4.6.1",
  snapshot_date = "2026-07-24",
  bioc_version = "3.23",
  packages = list(dplyr = list(Package = "dplyr", Version = "1.2.1", Source = "Repository", Repository = "CRAN"))
) {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  images <- file.path(root, "images")
  project <- file.path(root, "project")
  dir.create(images)
  dir.create(project)

  image <- fake_image(images)
  dir.create(file.path(project, "renv"))
  writeLines("# stub", file.path(project, "renv", "activate.R"))
  writeLines("{}", file.path(project, "renv", "settings.json"))
  writeLines(
    rprofile_template(snapshot_date, bioc_version, "noble"),
    file.path(project, ".Rprofile")
  )
  fake_renv_lock(project, r_version, snapshot_date, bioc_version, packages)

  env_lock <- list(
    project = "project",
    rmyknife_version = "0.4.0",
    created = "2026-07-24",
    r_version = r_version,
    snapshot_date = snapshot_date,
    bioc_version = bioc_version,
    bind = "/cephfs:/cephfs",
    images = list(describe_image(image, "primary", labels = list(), checksum = TRUE))
  )
  write_env_lock(env_lock, env_lock_path(project))

  list(root = root, project = project, image = image, env_lock = env_lock)
}

status_of <- function(report, check) {
  report$status[report$check == check]
}
