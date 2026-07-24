# Backup and restore of a reproducible project environment.

#' Path of the active renv library, relative to the project root
#'
#' renv keys libraries by platform and R version, so a project can carry trees
#' for several R versions. Only the active one belongs in a backup.
#'
#' @param path Project root
#' @keywords internal
active_library <- function(path) {
  library_path <- tryCatch(renv::paths$library(project = path), error = function(e) NULL)
  if (is.null(library_path) || !dir.exists(library_path)) {
    pattern <- file.path(path, "renv", "library", "*", paste0("R-", getRversion()$major, ".", getRversion()$minor))
    candidates <- Sys.glob(pattern)
    if (length(candidates) == 0) {
      return(NULL)
    }
    library_path <- candidates[1]
  }
  # renv::paths$library() points at the architecture subdirectory; back up the
  # R-version directory above it so the tree restores where renv expects it.
  relative <- sub(paste0("^", normalizePath(path, mustWork = FALSE), "/?"), "", normalizePath(library_path, mustWork = FALSE))
  sub("/x86_64[^/]*$", "", relative)
}

#' Instructions written into every backup
#' @param env_lock Environment lock
#' @param archive Archive file name
#' @param members Files contained in the archive
#' @keywords internal
restore_instructions <- function(env_lock, archive, members) {
  primary <- Filter(function(image) identical(image$role, "primary"), env_lock$images)[[1]]
  # Absent fields round-trip through JSON as NULL rather than NA.
  docker <- primary$docker %||% NA_character_
  sha256 <- primary$sha256 %||% NA_character_
  c(
    paste0("# Restoring ", env_lock$project),
    "",
    paste0("Backed up ", format(Sys.Date()), " by rmyknife ", env_lock$rmyknife_version, "."),
    "",
    "## What this archive pins",
    "",
    paste0("- R ", env_lock$r_version),
    paste0("- CRAN snapshot ", env_lock$snapshot_date, ", Bioconductor ", env_lock$bioc_version),
    paste0("- primary image `", primary$name, "`", if (!is.na(sha256)) paste0(" (sha256 ", substr(sha256, 1, 16), "…)") else ""),
    if (!is.na(docker)) paste0("- rebuildable from docker tag `", docker, "` if the image file is lost") else NULL,
    "",
    "## Contents",
    "",
    paste0("- ", members),
    "",
    "## Restore",
    "",
    "```r",
    paste0('rmyknife::project_restore("', archive, '", destination = "<target-dir>")'),
    "```",
    "",
    "Or by hand:",
    "",
    "```bash",
    paste0("tar --use-compress-program=zstd -xf ", archive),
    "sha256sum -c CHECKSUMS.sha256",
    paste0("cp images/", primary$name, " ", dirname(primary$path), "/"),
    "tar -xf library.tar.zst --use-compress-program=zstd    # if present",
    "```",
    "",
    "The binary library matches this exact image, so no compilation and no network",
    "are needed. `renv::restore()` is only the fallback when the library tarball",
    "was not included.",
    "",
    "## What this archive is not",
    "",
    "This is the compute environment only: the image and the package library. The",
    "project code, its history and its `.Rprofile` are not here, because they live",
    paste0("in git. Check the repository out at commit `", env_lock$git_revision %||% "unknown", "` to pair the two halves back up.")
  )
}

#' Archive the compute environment of a project
#'
#' Produces a single compressed archive holding the singularity images and the
#' built package library, with `renv.lock` to describe it. Together they restore
#' without compilation or network access, so a project stays restorable even when
#' its lockfile has entries that no longer resolve.
#'
#' The project code, its history and its `.Rprofile` are deliberately excluded:
#' they live in git. The archive records the commit the environment served, so
#' the two halves can be paired up again.
#'
#' @param path Project root.
#' @param destination Directory for the archive, relative to `path`.
#' @param include What to archive. Either or both of `"library"` and `"images"`.
#' @return Path to the archive, invisibly.
#' @export
#' @examples
#' \dontrun{
#'   project_backup()
#'   # quick snapshot without the multi-GB images
#'   project_backup(include = "library")
#' }
project_backup <- function(
  path = ".",
  destination = "backup",
  include = c("library", "images")
) {
  include <- match.arg(include, several.ok = TRUE)
  path <- normalizePath(path, mustWork = TRUE)
  env_lock <- read_env_lock(path)

  target <- file.path(path, destination)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)

  revision <- git_run(path, c("rev-parse", "--short", "HEAD"))
  revision <- if (length(revision) == 1) revision else "nogit"
  stem <- paste0(env_lock$project, "_", format(Sys.Date()), "_", revision)
  tarball <- file.path(target, paste0(stem, ".tar"))
  archive <- paste0(tarball, ".zst")
  if (file.exists(archive)) {
    stop("Backup already exists: ", archive)
  }

  staging <- file.path(tempdir(), paste0(stem, "-staging"))
  unlink(staging, recursive = TRUE)
  dir.create(staging, recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  # renv.lock is the only project file carried, because it describes the very
  # library being archived. Everything else the project needs lives in git.
  members <- character()
  if (file.exists(file.path(path, "renv.lock"))) {
    file.copy(file.path(path, "renv.lock"), file.path(staging, "renv.lock"), overwrite = TRUE)
    members <- "renv.lock"
  }

  library_relative <- NULL
  if ("library" %in% include) {
    library_relative <- active_library(path)
    if (is.null(library_relative)) {
      message("⚠️  No renv library found, skipping the library tarball")
    } else {
      members <- c(members, "library.tar.zst")
    }
  }
  if ("images" %in% include) {
    members <- c(members, paste0("images/", vapply(env_lock$images, function(image) image$name, character(1))))
  }
  members <- c(members, "MANIFEST.json", "CHECKSUMS.sha256", "RESTORE.md")

  manifest <- env_lock
  manifest$backup_date <- format(Sys.Date())
  manifest$git_revision <- revision
  manifest$contents <- members
  write_env_lock(manifest, file.path(staging, "MANIFEST.json"))
  writeLines(restore_instructions(manifest, basename(archive), members), file.path(staging, "RESTORE.md"))

  if (!is.null(library_relative)) {
    message("Archiving package library, this takes a while")
    # -h dereferences: renv/library is a tree of symlinks into the shared cache,
    # so without this the archive would contain nothing but dangling links.
    system2("tar", c(
      "-c", "-h", "--use-compress-program=zstd",
      "-f", shQuote(file.path(staging, "library.tar.zst")),
      "-C", shQuote(path), shQuote(library_relative)
    ))
  }

  checksums <- character()
  for (file in list.files(staging, recursive = TRUE)) {
    checksums <- c(checksums, paste0(sha256_file(file.path(staging, file)), "  ", file))
  }
  if ("images" %in% include) {
    for (image in env_lock$images) {
      if (file.exists(image$path)) {
        checksums <- c(checksums, paste0(image$sha256, "  images/", image$name))
      }
    }
  }
  writeLines(checksums, file.path(staging, "CHECKSUMS.sha256"))

  message("Building archive")
  system2("tar", c("-c", "-f", shQuote(tarball), "-C", shQuote(staging), "."))

  if ("images" %in% include) {
    for (image in env_lock$images) {
      if (!file.exists(image$path)) {
        message("⚠️  Image missing, not archived: ", image$path)
        next
      }
      message("Adding image ", image$name, " (", round(image$bytes / 1e9, 2), " GB)")
      system2("tar", c(
        "-r", "-f", shQuote(tarball),
        "--transform", shQuote("s,^,images/,"),
        "-C", shQuote(dirname(image$path)), shQuote(image$name)
      ))
    }
  }

  message("Compressing")
  system2("zstd", c("-T0", "-q", "--rm", shQuote(tarball), "-o", shQuote(archive)))

  message("✅ ", archive, " (", round(file.size(archive) / 1e9, 2), " GB)")
  invisible(archive)
}

#' Restore a compute environment from a backup
#'
#' Unpacks the archive, verifies every member against the recorded checksums and
#' reports what was recovered. Images are unpacked into `images/` for you to move
#' into place; nothing outside `destination` is written.
#'
#' This restores the environment only. Check the project code out of git at the
#' commit named in the manifest to pair the two halves back up.
#'
#' @param archive Path to an archive produced by [project_backup()].
#' @param destination Directory to restore into. Must not already exist.
#' @param verify Whether to check members against `CHECKSUMS.sha256`.
#' @return Path to the restored directory, invisibly.
#' @export
#' @examples
#' \dontrun{
#'   project_restore("backup/myproject_2026-07-24_a1b2c3d.tar.zst", destination = "restored")
#' }
project_restore <- function(archive, destination, verify = TRUE) {
  if (!file.exists(archive)) {
    stop("Archive not found: ", archive)
  }
  if (dir.exists(destination)) {
    stop("Destination already exists: ", destination)
  }
  dir.create(destination, recursive = TRUE)
  destination <- normalizePath(destination, mustWork = TRUE)

  message("Extracting ", basename(archive))
  system2("tar", c("--use-compress-program=zstd", "-x", "-f", shQuote(normalizePath(archive)), "-C", shQuote(destination)))

  manifest_path <- file.path(destination, "MANIFEST.json")
  if (!file.exists(manifest_path)) {
    stop("Archive has no MANIFEST.json, this does not look like a project backup")
  }
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)

  if (verify) {
    message("Verifying checksums")
    checksums <- readLines(file.path(destination, "CHECKSUMS.sha256"), warn = FALSE)
    bad <- character()
    for (line in checksums) {
      expected <- sub("[[:space:]].*$", "", line)
      file <- sub("^[^[:space:]]+[[:space:]]+", "", line)
      full <- file.path(destination, file)
      if (!file.exists(full)) {
        bad <- c(bad, paste0(file, " (missing)"))
      } else if (!identical(sha256_file(full), expected)) {
        bad <- c(bad, paste0(file, " (checksum mismatch)"))
      }
    }
    if (length(bad) > 0) {
      stop("Archive is damaged:\n  ", paste(bad, collapse = "\n  "))
    }
    message("✅ All members verified")
  }

  library_tarball <- file.path(destination, "library.tar.zst")
  if (file.exists(library_tarball)) {
    message("Unpacking package library")
    system2("tar", c("--use-compress-program=zstd", "-x", "-f", shQuote(library_tarball), "-C", shQuote(destination)))
    file.remove(library_tarball)
    message("✅ Library restored, no compilation needed")
  } else {
    message("⚠️  No library tarball; rebuild with renv::restore() after checking project_verify()")
  }

  for (image in manifest$images) {
    restored <- file.path(destination, "images", image$name)
    if (file.exists(restored)) {
      message("✅ Image ", image$name, " restored; move it to ", dirname(image$path))
    } else {
      message("⚠️  Image ", image$name, " not in this archive; it belongs at ", image$path)
      if (!is.null(image$docker) && !is.na(image$docker)) {
        message("   May be rebuildable from docker tag ", image$docker)
      }
    }
  }

  revision <- manifest$git_revision %||% "unknown"
  if (!identical(revision, "nogit")) {
    message("ℹ️  Environment only. Check the project out of git at commit ", revision)
  }

  message("\n✅ Restored to ", destination, ". See RESTORE.md.")
  invisible(destination)
}
