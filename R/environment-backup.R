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

#' Does the project track files with git-lfs?
#'
#' Matters for backups: `git bundle` carries the whole history, but lfs blobs
#' live outside git, so only their pointers make it into the archive.
#'
#' @param path Project root
#' @keywords internal
uses_git_lfs <- function(path) {
  attributes <- file.path(path, ".gitattributes")
  if (!file.exists(attributes)) {
    return(FALSE)
  }
  any(grepl("filter=lfs", readLines(attributes, warn = FALSE), fixed = TRUE))
}

#' Instructions written into every backup
#' @param env_lock Environment lock
#' @param archive Archive file name
#' @param members Files contained in the archive
#' @param lfs Whether the project tracks files with git-lfs
#' @keywords internal
restore_instructions <- function(env_lock, archive, members, lfs = FALSE) {
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
    "git clone repo.bundle <project>                        # if present",
    "```",
    "",
    "The binary library matches this exact image, so no compilation and no network",
    "are needed. `renv::restore()` is only the fallback when the library tarball",
    "was not included.",
    if (lfs) {
      c(
        "",
        "## git-lfs",
        "",
        "This project tracks files with git-lfs. `repo.bundle` holds the complete",
        "history, but lfs blobs live outside git and are **not** in this archive:",
        "only their pointers are. Clone with `GIT_LFS_SKIP_SMUDGE=1 git clone",
        "repo.bundle <project>`, then fetch the blobs from the lfs remote while it",
        "still exists. Those files are deliverables, not part of the environment."
      )
    } else {
      NULL
    }
  )
}

#' Archive everything needed to restore a project environment
#'
#' Produces a single compressed archive holding the environment lock, the renv
#' configuration, the built package library, the git history and the singularity
#' images. Paired with the archived image, the binary library restores without
#' compilation or network access, so a project stays restorable even when its
#' lockfile has entries that no longer resolve.
#'
#' @param path Project root.
#' @param destination Directory for the archive, relative to `path`.
#' @param include What to archive besides the environment lock and renv
#'   configuration. Any of `"library"`, `"git"` and `"images"`.
#' @return Path to the archive, invisibly.
#' @export
#' @examples
#' \dontrun{
#'   project_backup()
#'   # quick snapshot without the multi-GB images
#'   project_backup(include = c("library", "git"))
#' }
project_backup <- function(
  path = ".",
  destination = "backup",
  include = c("library", "git", "images")
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

  configuration <- c(env_lock_file, "renv.lock", ".Rprofile", "analysis/.Rprofile", "renv/settings.json", "renv/activate.R")
  members <- character()
  for (file in configuration) {
    source <- file.path(path, file)
    if (!file.exists(source)) {
      next
    }
    dir.create(file.path(staging, dirname(file)), recursive = TRUE, showWarnings = FALSE)
    file.copy(source, file.path(staging, file), overwrite = TRUE)
    members <- c(members, file)
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
  if ("git" %in% include && dir.exists(file.path(path, ".git"))) {
    members <- c(members, "repo.bundle")
  }
  if ("images" %in% include) {
    members <- c(members, paste0("images/", vapply(env_lock$images, function(image) image$name, character(1))))
  }
  members <- c(members, "MANIFEST.json", "CHECKSUMS.sha256", "RESTORE.md")

  manifest <- env_lock
  manifest$backup_date <- format(Sys.Date())
  manifest$git_revision <- revision
  manifest$contents <- members
  lfs <- uses_git_lfs(path)
  write_env_lock(manifest, file.path(staging, "MANIFEST.json"))
  writeLines(restore_instructions(env_lock, basename(archive), members, lfs), file.path(staging, "RESTORE.md"))

  if ("git" %in% include && dir.exists(file.path(path, ".git"))) {
    if (lfs) {
      message("⚠️  Project uses git-lfs; the bundle carries pointers, not blobs. See RESTORE.md.")
    }
    message("Bundling git history")
    system2("git", c("-C", shQuote(path), "bundle", "create", shQuote(file.path(staging, "repo.bundle")), "--all"),
      stdout = FALSE, stderr = FALSE)
  }

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

#' Restore a project environment from a backup
#'
#' Unpacks the archive, verifies every member against the recorded checksums and
#' reports what was recovered. Images are unpacked into `images/` for you to move
#' into place; nothing outside `destination` is written.
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

  bundle <- file.path(destination, "repo.bundle")
  if (file.exists(bundle)) {
    message("✅ Git history in repo.bundle; clone it with: git clone repo.bundle <project>")
  }

  message("\n✅ Restored to ", destination, ". See RESTORE.md.")
  invisible(destination)
}
