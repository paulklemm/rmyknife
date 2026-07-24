# Verification of a reproducible project environment.

#' One row of the verification report
#' @param check Name of the check
#' @param status One of "ok", "warn", "fail"
#' @param detail Human readable detail
#' @keywords internal
check_row <- function(check, status, detail) {
  tibble::tibble(check = check, status = status, detail = detail)
}

#' Status of a condition, as "ok" or a chosen failure level
#' @param ok Logical condition
#' @param level Status to use when `ok` is FALSE
#' @keywords internal
status_if <- function(ok, level = "fail") {
  if (isTRUE(ok)) "ok" else level
}

#' Is a package repository serving its index?
#'
#' Probes the package index rather than the repository root: a Posit Package
#' Manager snapshot root answers 400, so only the index tells us whether the
#' snapshot is actually still being served.
#'
#' @param repo Repository URL
#' @keywords internal
repo_ok <- function(repo) {
  tryCatch(
    {
      handle <- curl::new_handle(nobody = TRUE, followlocation = TRUE, timeout = 20)
      status <- curl::curl_fetch_memory(paste0(repo, "/src/contrib/PACKAGES.gz"), handle = handle)$status_code
      status < 400
    },
    error = function(e) FALSE
  )
}

#' Run git in a project and return its output
#' @param path Project root
#' @param args Arguments to git
#' @keywords internal
git_run <- function(path, args) {
  suppressWarnings(tryCatch(
    system2("git", c("-C", shQuote(path), args), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  ))
}

#' Packages in a lockfile that cannot be restored from a repository
#'
#' These are the entries a network restore would fail on: packages renv could
#' not attribute to a repository, packages installed from a local path or the
#' cellar, and git remotes without a pinned commit.
#'
#' @param lockfile Path to renv.lock
#' @keywords internal
unresolvable_packages <- function(lockfile) {
  lock <- jsonlite::fromJSON(lockfile, simplifyVector = FALSE)
  unresolvable <- character()
  for (name in names(lock$Packages)) {
    package <- lock$Packages[[name]]
    source <- package$Source
    if (is.null(source) || tolower(source) %in% c("unknown", "local", "cellar")) {
      unresolvable <- c(unresolvable, name)
    } else if (!is.null(package$RemoteType) && is.null(package$RemoteSha)) {
      unresolvable <- c(unresolvable, name)
    }
  }
  unresolvable
}

#' Check that every recorded image is present and unchanged
#' @param env_lock Environment lock
#' @param deep Whether to verify checksums
#' @keywords internal
verify_images <- function(env_lock, deep) {
  rows <- list()
  for (image in env_lock$images) {
    name <- paste0("image:", image$name)
    if (!file.exists(image$path)) {
      recovery <- if (is.null(image$docker) || is.na(image$docker)) {
        "no docker tag recorded, this container is unrecoverable"
      } else {
        paste0("possibly rebuildable from ", image$docker)
      }
      rows[[length(rows) + 1L]] <- check_row(name, "fail", paste0("missing at ", image$path, "; ", recovery))
      next
    }
    bytes <- as.numeric(file.size(image$path))
    if (!isTRUE(all.equal(bytes, as.numeric(image$bytes)))) {
      rows[[length(rows) + 1L]] <- check_row(
        name, "fail", paste0("size changed: ", bytes, " bytes, expected ", image$bytes)
      )
      next
    }
    if (deep && !is.null(image$sha256) && !is.na(image$sha256)) {
      if (!identical(sha256_file(image$path), image$sha256)) {
        rows[[length(rows) + 1L]] <- check_row(name, "fail", "sha256 mismatch, the image has changed")
        next
      }
      rows[[length(rows) + 1L]] <- check_row(name, "ok", "present, sha256 verified")
      next
    }
    rows[[length(rows) + 1L]] <- check_row(name, "ok", "present, size matches")
  }
  rows
}

#' Check that everything needed to restore this project is in place
#'
#' Reports two independent notions of restorability. *From backup* is offline and
#' exact: the archived containers plus the binary library. *From lockfile* is a
#' network rebuild from scratch and needs every package to resolve, which is
#' commonly partial for a project converted from a pre-renv state. A project can
#' be perfectly restorable from its backup while its lockfile is still messy.
#'
#' @param path Project root.
#' @param deep Whether to verify image checksums. Several GB per image, so slow.
#' @param network Whether to probe the pinned repositories.
#' @param strict Whether to raise an error when any check fails. Used by the
#'   makefile target.
#' @return A tibble of `check`, `status` and `detail`, invisibly.
#' @export
#' @examples
#' \dontrun{
#'   project_verify()
#'   project_verify(deep = TRUE, strict = TRUE)
#' }
project_verify <- function(path = ".", deep = FALSE, network = TRUE, strict = FALSE) {
  path <- normalizePath(path, mustWork = TRUE)
  rows <- list()

  env_lock <- tryCatch(read_env_lock(path), error = function(e) NULL)
  if (is.null(env_lock)) {
    report <- check_row("environment.lock", "fail", "missing or unparseable, run project_init()")
    print_verify(report)
    if (strict) {
      stop("Verification failed")
    }
    return(invisible(report))
  }
  rows[[length(rows) + 1L]] <- check_row("environment.lock", "ok", "present")

  required <- c(".Rprofile", "renv.lock", "renv/activate.R", "renv/settings.json")
  missing <- required[!file.exists(file.path(path, required))]
  rows[[length(rows) + 1L]] <- check_row(
    "project files",
    status_if(length(missing) == 0),
    if (length(missing) == 0) "all present" else paste("missing:", paste(missing, collapse = ", "))
  )

  rows <- c(rows, verify_images(env_lock, deep))

  running <- Sys.getenv("APPTAINER_CONTAINER")
  primary <- Filter(function(image) identical(image$role, "primary"), env_lock$images)[[1]]
  if (nzchar(running)) {
    # Apptainer reports the path it was given, so a session started through a
    # `latest/` symlink names the symlink. Compare resolved paths, or every such
    # session would look like it is running the wrong container.
    resolved <- normalizePath(running, mustWork = FALSE)
    matches <- identical(resolved, primary$path)
    rows[[length(rows) + 1L]] <- check_row(
      "running image",
      status_if(matches, "warn"),
      if (matches && identical(resolved, running)) {
        "matches the recorded primary image"
      } else if (matches) {
        paste0("matches the recorded primary image, via ", running)
      } else {
        paste0("running ", resolved, ", recorded ", primary$path)
      }
    )
  } else if (in_container()) {
    # Inside a container that did not tell us which image file it came from.
    # The build date identifies it just as well.
    running_build <- parse_build_date(label_value(image_labels_self(), "org.label-schema.build-date"))
    recorded_build <- primary$build_date %||% NA_character_
    confirmed <- !is.na(running_build) && identical(running_build, recorded_build)
    rows[[length(rows) + 1L]] <- check_row(
      "running image",
      status_if(confirmed, "warn"),
      if (confirmed) {
        paste0("matches the recorded primary image, by build date ", running_build)
      } else if (is.na(running_build) || is.na(recorded_build)) {
        "APPTAINER_CONTAINER is unset and there is no build date to compare"
      } else {
        paste0("running a container built ", running_build, ", recorded ", recorded_build)
      }
    )
  }

  rows[[length(rows) + 1L]] <- check_row(
    "R version",
    status_if(identical(as.character(getRversion()), env_lock$r_version)),
    paste0("running ", getRversion(), ", recorded ", env_lock$r_version)
  )

  if (!is.null(primary$base) && !is.na(primary$base)) {
    codename <- sub("^.*:", "", primary$base)
    library_root <- file.path(path, "renv", "library")
    trees <- if (dir.exists(library_root)) list.dirs(library_root, recursive = FALSE, full.names = FALSE) else character()
    matches <- length(trees) == 0 || any(grepl(codename, trees, fixed = TRUE))
    rows[[length(rows) + 1L]] <- check_row(
      "library platform",
      status_if(matches, "warn"),
      if (matches) paste0("library matches ", codename) else paste0("library trees ", paste(trees, collapse = ", "), " do not match ", codename)
    )
  }

  synchronized <- tryCatch(
    isTRUE(renv::status(project = path)$synchronized),
    error = function(e) NA
  )
  rows[[length(rows) + 1L]] <- check_row(
    "renv status",
    if (is.na(synchronized)) "warn" else status_if(synchronized),
    if (isTRUE(synchronized)) "library and lockfile agree" else "library and lockfile disagree, see renv::status()"
  )

  lock <- lockfile_info(path)
  if (!is.null(lock)) {
    consistent <- identical(lock$r_version, env_lock$r_version) &&
      identical(lock$snapshot_date, env_lock$snapshot_date) &&
      identical(lock$bioc_version, env_lock$bioc_version)
    rows[[length(rows) + 1L]] <- check_row(
      "lockfile consistency",
      status_if(consistent),
      if (consistent) {
        "renv.lock agrees with environment.lock"
      } else {
        paste0(
          "renv.lock says R ", lock$r_version, " / ", lock$snapshot_date, " / Bioc ", lock$bioc_version,
          "; environment.lock says R ", env_lock$r_version, " / ", env_lock$snapshot_date, " / Bioc ", env_lock$bioc_version
        )
      }
    )

    rprofile_lines <- readLines(file.path(path, ".Rprofile"), warn = FALSE)
    pinned <- any(grepl(env_lock$snapshot_date, rprofile_lines, fixed = TRUE)) &&
      any(grepl(env_lock$bioc_version, rprofile_lines, fixed = TRUE))
    rows[[length(rows) + 1L]] <- check_row(
      ".Rprofile pins",
      status_if(pinned),
      if (pinned) "pins match environment.lock" else "pins missing or disagree with environment.lock"
    )
  }

  images_present <- all(vapply(env_lock$images, function(image) file.exists(image$path), logical(1)))
  backups <- list.files(file.path(path, "backup"), pattern = "\\.tar\\.(zst|gz)$", full.names = TRUE)
  rows[[length(rows) + 1L]] <- check_row(
    "restorable from backup",
    status_if(length(backups) > 0 || images_present, "warn"),
    if (length(backups) > 0) {
      paste0(length(backups), " backup(s) in backup/")
    } else if (images_present) {
      "no backup yet, but all images are present so one can be made now"
    } else {
      "no backup and images are missing"
    }
  )

  if (file.exists(file.path(path, "renv.lock"))) {
    unresolvable <- unresolvable_packages(file.path(path, "renv.lock"))
    rows[[length(rows) + 1L]] <- check_row(
      "restorable from lockfile",
      status_if(length(unresolvable) == 0, "warn"),
      if (length(unresolvable) == 0) {
        "every package resolves to a repository"
      } else {
        paste0(length(unresolvable), " package(s) will not restore over the network: ", paste(unresolvable, collapse = ", "))
      }
    )
  }

  if (network) {
    repos <- pinned_repos(env_lock$snapshot_date, env_lock$bioc_version, sub("^.*:", "", primary$base %||% "noble"))
    unreachable <- names(repos)[!vapply(repos, repo_ok, logical(1))]
    rows[[length(rows) + 1L]] <- check_row(
      "repositories reachable",
      status_if(length(unreachable) == 0, "warn"),
      if (length(unreachable) == 0) "all pinned repositories respond" else paste("unreachable:", paste(unreachable, collapse = ", "))
    )
  }

  makefile <- c(file.path(path, "makefile"), file.path(path, "Makefile"))
  makefile <- makefile[file.exists(makefile)]
  if (length(makefile) > 0) {
    lines <- readLines(makefile[1], warn = FALSE)
    singularity <- grep("^SINGULARITY", lines, value = TRUE)
    referenced <- unlist(regmatches(singularity, gregexpr("[^[:space:]]+\\.(simg|sif)", singularity)))
    # A makefile may name the image through a `latest/` symlink, so resolve
    # before comparing. Paths built from a make or shell variable cannot be
    # resolved here, and are reported as unconfirmable rather than wrong.
    variable <- any(grepl("[$]", referenced))
    matches <- primary$path %in% normalizePath(referenced, mustWork = FALSE)
    rows[[length(rows) + 1L]] <- check_row(
      "makefile image",
      status_if(matches, "warn"),
      if (matches) {
        "points at the recorded primary image"
      } else if (variable) {
        "SINGULARITY line builds the path from a variable, cannot confirm it"
      } else {
        "SINGULARITY line does not reference the recorded image"
      }
    )
  }

  if (dir.exists(file.path(path, ".git"))) {
    tracked <- git_run(path, c("ls-files", "--", ".Rprofile", "renv.lock", env_lock_file))
    expected <- c(".Rprofile", "renv.lock", env_lock_file)
    untracked <- setdiff(expected, tracked)
    rows[[length(rows) + 1L]] <- check_row(
      "git tracked",
      status_if(length(untracked) == 0, "warn"),
      if (length(untracked) == 0) "pins and lock are committed" else paste("not tracked:", paste(untracked, collapse = ", "))
    )
    dirty <- git_run(path, c("status", "--porcelain", "--", ".Rprofile", "renv.lock", env_lock_file))
    rows[[length(rows) + 1L]] <- check_row(
      "git clean",
      status_if(length(dirty) == 0, "warn"),
      if (length(dirty) == 0) "no uncommitted changes to pins" else "pins have uncommitted changes"
    )
  }

  if (length(backups) > 0) {
    newest <- max(file.mtime(backups))
    lock_mtime <- file.mtime(file.path(path, "renv.lock"))
    fresh <- is.na(lock_mtime) || newest >= lock_mtime
    rows[[length(rows) + 1L]] <- check_row(
      "backup freshness",
      status_if(fresh, "warn"),
      if (fresh) "newest backup is up to date" else "renv.lock changed since the newest backup"
    )
  }

  report <- dplyr::bind_rows(rows)
  print_verify(report)
  if (strict && any(report$status == "fail")) {
    stop("Verification failed: ", sum(report$status == "fail"), " check(s)")
  }
  invisible(report)
}

#' Default value for NULL
#' @param x Value
#' @param y Fallback
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

#' Print a verification report
#' @param report Tibble as returned by [project_verify()]
#' @keywords internal
print_verify <- function(report) {
  width <- max(nchar(report$check))
  for (row in seq_len(nrow(report))) {
    message(sprintf(
      "%s %-*s  %s",
      status_tag(report$status[row]),
      width,
      report$check[row],
      report$detail[row]
    ))
  }
  failed <- sum(report$status == "fail")
  warned <- sum(report$status == "warn")
  message(sprintf(
    "\n%d ok, %d warning(s), %d failure(s)",
    sum(report$status == "ok"), warned, failed
  ))
}
