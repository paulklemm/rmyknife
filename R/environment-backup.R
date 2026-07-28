# Backup and restore of a reproducible project environment.

#' Run an archiving command, stopping when it fails
#'
#' A backup that reports success without producing an archive is worse than one
#' that fails, because the failure is only discovered when the archive is needed.
#'
#' @param command Command to run
#' @param args Arguments to the command
#' @keywords internal
run_or_stop <- function(command, args) {
  status <- system2(command, args)
  if (!identical(as.integer(status), 0L)) {
    stop(command, " failed with exit status ", status)
  }
  invisible(status)
}

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
  # Compared as text rather than as a pattern, so a project path holding a regex
  # metacharacter still matches. Tried before resolving symlinks, because
  # `renv/library` is often a link to faster storage: that still archives
  # correctly, since tar is given a project-relative member and dereferences it
  # (-h). Only a library that cannot be named relative to the project at all is
  # a problem, and it is a silent one -- tar strips the leading `/` from an
  # absolute member with a warning and exit status 0, so the library would
  # restore somewhere renv never looks.
  root <- normalizePath(path, mustWork = FALSE)
  under_root <- function(candidate) startsWith(candidate, paste0(root, "/"))
  full <- if (under_root(library_path)) library_path else normalizePath(library_path, mustWork = FALSE)
  if (!under_root(full)) {
    stop(
      "The renv library lives outside the project, so it cannot be archived relative to it:\n",
      "  library: ", full, "\n",
      "  project: ", root, "\n",
      "Unset RENV_PATHS_LIBRARY, or archive the images alone with include = \"images\"."
    )
  }
  relative <- substring(full, nchar(root) + 2L)
  # renv::paths$library() points at the architecture subdirectory; back up the
  # R-version directory above it so the tree restores where renv expects it.
  # Recognised by the parent's name rather than by the architecture's, which is
  # not always x86_64.
  if (grepl("^R-[0-9]", basename(dirname(relative)))) {
    relative <- dirname(relative)
  }
  relative
}

#' Checksums of the images going into a backup
#'
#' Recomputed rather than copied out of `environment.lock`, because the checksum
#' file has to describe what is really in the archive. An image that no longer
#' matches its record stops the backup: the project is then not the environment
#' it claims to be, and archiving it would bake that disagreement in.
#'
#' @param images Image entries that will be archived
#' @keywords internal
image_checksums <- function(images) {
  vapply(images, function(image) {
    actual <- sha256_file(image$path)
    recorded <- lock_field(image$sha256, NA_character_)
    if (!is.na(recorded) && !identical(actual, recorded)) {
      stop(
        "Image has changed since it was recorded: ", image$path, "\n",
        "  recorded: ", recorded, "\n",
        "  on disk:  ", actual, "\n",
        "Diagnose with project_verify(deep = TRUE) before backing up."
      )
    }
    paste0(actual, "  images/", image$name)
  }, character(1))
}

#' Instructions written into every backup
#' @param env_lock Environment lock
#' @param archive Archive file name
#' @param members Files contained in the archive
#' @keywords internal
restore_instructions <- function(env_lock, archive, members) {
  primary <- primary_image(env_lock)
  # Absent fields round-trip through JSON as NULL rather than NA.
  docker <- lock_field(primary$docker, NA_character_)
  sha256 <- lock_field(primary$sha256, NA_character_)
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
    paste0("in git. Check the repository out at commit `", lock_field(env_lock$git_revision, "unknown"), "` to pair the two halves back up.")
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
#' Images are checksummed before they are archived, and a backup of an image
#' that no longer matches `environment.lock` is refused rather than written.
#'
#' @param path Project root.
#' @param destination Directory for the archive, relative to `path`. Staging
#'   happens here too, so it needs room for roughly twice the finished archive.
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

  # Staged beside the archive rather than in tempdir(): the library tarball is
  # gigabytes, and the destination is a directory the caller chose and that has
  # room for the finished archive anyway.
  staging <- file.path(target, paste0(stem, "-staging"))
  # Every leftover of this project, not just this stem's: a hard-killed run
  # cannot clean up after itself, and staging no longer lives in tempdir()
  # where the session exit would have taken it. Matched on the project name as
  # well as the suffix, because `destination` need not be a directory this
  # function owns.
  unlink(Sys.glob(file.path(target, paste0(env_lock$project, "_*-staging"))), recursive = TRUE)
  dir.create(staging, recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  # The uncompressed tarball is an intermediate; a failed run should not leave
  # gigabytes of it behind.
  on.exit(unlink(tarball), add = TRUE)

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
      status_message("warn", "No renv library found, skipping the library tarball")
    } else {
      members <- c(members, "library.tar.zst")
    }
  }
  # Settled before the manifest is written, so that what it lists is what the
  # archive holds rather than what the environment lock hoped for.
  archived_images <- list()
  if ("images" %in% include) {
    for (image in env_lock$images) {
      if (file.exists(image$path)) {
        archived_images[[length(archived_images) + 1L]] <- image
      } else {
        status_message("warn", "Image missing, not archived: ", image$path)
      }
    }
    members <- c(members, paste0("images/", vapply(archived_images, function(image) image$name, character(1))))
  }
  members <- c(members, "MANIFEST.json", "CHECKSUMS.sha256", "RESTORE.md")

  # Up front, so an image that drifted from its record is caught before rather
  # than after the minutes the library tarball costs.
  if (length(archived_images) > 0) {
    message("Checksumming images")
  }
  image_lines <- image_checksums(archived_images)

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
    run_or_stop("tar", c(
      "-c", "-h", "--use-compress-program=zstd",
      "-f", shQuote(file.path(staging, "library.tar.zst")),
      "-C", shQuote(path), shQuote(library_relative)
    ))
  }

  checksums <- character()
  for (file in list.files(staging, recursive = TRUE)) {
    checksums <- c(checksums, paste0(sha256_file(file.path(staging, file)), "  ", file))
  }
  checksums <- c(checksums, image_lines)
  writeLines(checksums, file.path(staging, "CHECKSUMS.sha256"))

  message("Building archive")
  run_or_stop("tar", c("-c", "-f", shQuote(tarball), "-C", shQuote(staging), "."))

  for (image in archived_images) {
    message("Adding image ", image$name, " (", round(image$bytes / 1e9, 2), " GB)")
    run_or_stop("tar", c(
      "-r", "-f", shQuote(tarball),
      "--transform", shQuote("s,^,images/,"),
      "-C", shQuote(dirname(image$path)), shQuote(image$name)
    ))
  }

  message("Compressing")
  run_or_stop("zstd", c("-T0", "-q", "--rm", shQuote(tarball), "-o", shQuote(archive)))

  status_message("ok", archive, " (", round(file.size(archive) / 1e9, 2), " GB)")
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
  # Removed again unless the restore runs to the end, so a failed attempt does
  # not block the retry on "Destination already exists".
  completed <- FALSE
  on.exit(if (!completed) unlink(destination, recursive = TRUE), add = TRUE)
  destination <- normalizePath(destination, mustWork = TRUE)

  message("Extracting ", basename(archive))
  run_or_stop("tar", c("--use-compress-program=zstd", "-x", "-f", shQuote(normalizePath(archive)), "-C", shQuote(destination)))

  manifest_path <- file.path(destination, "MANIFEST.json")
  if (!file.exists(manifest_path)) {
    stop("Archive has no MANIFEST.json, this does not look like a project backup")
  }
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)

  if (verify) {
    checksum_path <- file.path(destination, "CHECKSUMS.sha256")
    if (!file.exists(checksum_path)) {
      stop("Archive has no CHECKSUMS.sha256, this does not look like a project backup")
    }
    message("Verifying checksums")
    checksums <- readLines(checksum_path, warn = FALSE)
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
    status_message("ok", "All members verified")
  }

  library_tarball <- file.path(destination, "library.tar.zst")
  if (file.exists(library_tarball)) {
    message("Unpacking package library")
    run_or_stop("tar", c("--use-compress-program=zstd", "-x", "-f", shQuote(library_tarball), "-C", shQuote(destination)))
    file.remove(library_tarball)
    status_message("ok", "Library restored, no compilation needed")
  } else {
    status_message("warn", "No library tarball; rebuild with renv::restore() after checking project_verify()")
  }

  for (image in manifest$images) {
    restored <- file.path(destination, "images", image$name)
    if (file.exists(restored)) {
      status_message("ok", "Image ", image$name, " restored; move it to ", dirname(image$path))
    } else {
      status_message("warn", "Image ", image$name, " not in this archive; it belongs at ", image$path)
      if (!is.null(image$docker) && !is.na(image$docker)) {
        message("   May be rebuildable from docker tag ", image$docker)
      }
    }
  }

  revision <- lock_field(manifest$git_revision, "unknown")
  if (!identical(revision, "nogit")) {
    status_message("info", "Environment only. Check the project out of git at commit ", revision)
  }

  message("")
  status_message("ok", "Restored to ", destination, ". See RESTORE.md.")
  completed <- TRUE
  invisible(destination)
}
