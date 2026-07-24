# Reproducible project environments.
#
# Three layers pin an analysis project: the singularity image (OS, R, system
# libraries), the R package library (renv.lock plus dated repositories), and the
# code (git). renv already covers the middle layer well. These functions add the
# container layer as a first-class, checksummed, in-repo fact, recorded in
# `environment.lock`, and make the whole thing verifiable and archivable.

#' Name of the environment lock file
#' @keywords internal
env_lock_file <- "environment.lock"

#' Path to the environment lock file of a project
#' @param path Project root
#' @keywords internal
env_lock_path <- function(path = ".") {
  file.path(path, env_lock_file)
}

#' Read a label, returning NA when absent or empty
#' @param labels Named list of image labels
#' @param key Label name
#' @keywords internal
label_value <- function(labels, key) {
  value <- labels[[key]]
  if (is.null(value) || !nzchar(value)) {
    return(NA_character_)
  }
  as.character(value)
}

#' sha256 checksum of a file
#' @param path File to checksum
#' @keywords internal
sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = FALSE)
  if (length(out) != 1) {
    stop("Could not compute sha256 for ", path)
  }
  sub("[[:space:]].*$", "", out)
}

#' Convert a singularity `build-date` label to an ISO date
#'
#' Labels look like `Tuesday_21_July_2026_10:1:14_UTC`. `month.name` is always
#' English in base R, so this is locale-independent.
#'
#' @param x Value of the `org.label-schema.build-date` label
#' @keywords internal
parse_build_date <- function(x) {
  if (is.na(x)) {
    return(NA_character_)
  }
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  if (length(parts) < 4) {
    return(NA_character_)
  }
  day <- suppressWarnings(as.integer(parts[2]))
  month <- match(parts[3], month.name)
  year <- suppressWarnings(as.integer(parts[4]))
  if (anyNA(c(day, month, year))) {
    return(NA_character_)
  }
  sprintf("%04d-%02d-%02d", year, month, day)
}

#' Is this session running inside a singularity container?
#'
#' Tests for the metadata directory apptainer places in every container, rather
#' than for `APPTAINER_CONTAINER`. The environment variable answers a different
#' question -- *which* image -- and some launchers drop it while the session is
#' genuinely still inside the container.
#'
#' @param marker Directory that marks a container. Exposed for testing.
#' @keywords internal
in_container <- function(marker = "/.singularity.d") {
  dir.exists(marker)
}

#' Labels of the image the current session runs in
#'
#' Apptainer exposes the image metadata inside the container, so this needs no
#' apptainer binary and no subprocess.
#'
#' @keywords internal
image_labels_self <- function() {
  labels <- "/.singularity.d/labels.json"
  if (!file.exists(labels)) {
    return(list())
  }
  jsonlite::fromJSON(labels, simplifyVector = TRUE)
}

#' Labels of an image we are not running in
#'
#' Falls back to an empty list when apptainer is unavailable or the image cannot
#' be inspected; the checksum remains the authoritative identity either way.
#'
#' @param image Path to the image
#' @keywords internal
image_labels_inspect <- function(image) {
  out <- suppressWarnings(tryCatch(
    system2("apptainer", c("inspect", shQuote(image)), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  ))
  out <- out[grepl(":", out, fixed = TRUE)]
  if (length(out) == 0) {
    return(list())
  }
  keys <- sub(":.*$", "", out)
  values <- trimws(sub("^[^:]*:", "", out))
  stats::setNames(as.list(values), keys)
}

#' Describe a singularity image for the environment lock
#' @param image Path to the image
#' @param role Either "primary" or "aux"
#' @param labels Pre-read labels, or NULL to inspect the image
#' @param checksum Whether to compute the sha256 checksum
#' @keywords internal
describe_image <- function(image, role, labels = NULL, checksum = TRUE) {
  if (!file.exists(image)) {
    stop("Image not found: ", image)
  }
  image <- normalizePath(image, mustWork = TRUE)
  if (is.null(labels)) {
    labels <- image_labels_inspect(image)
  }
  r_version <- label_value(labels, "org.opencontainers.image.version")
  list(
    role = role,
    path = image,
    name = basename(image),
    bytes = as.numeric(file.size(image)),
    sha256 = if (checksum) sha256_file(image) else NA_character_,
    docker = label_value(labels, "org.label-schema.usage.singularity.deffile.from"),
    build_date = parse_build_date(label_value(labels, "org.label-schema.build-date")),
    r_version = sub("^R-", "", r_version),
    base = label_value(labels, "org.opencontainers.image.base.name")
  )
}

#' Ubuntu codename of the image, used for the binary CRAN repository path
#' @param labels Named list of image labels
#' @keywords internal
base_codename <- function(labels) {
  base <- label_value(labels, "org.opencontainers.image.base.name")
  if (is.na(base)) {
    return("noble")
  }
  codename <- sub("^.*:", "", base)
  if (!nzchar(codename)) "noble" else codename
}

#' Pinned repository vector
#' @param snapshot_date CRAN snapshot date
#' @param bioc_version Bioconductor release
#' @param codename Ubuntu codename for binary packages
#' @keywords internal
pinned_repos <- function(snapshot_date, bioc_version, codename = "noble") {
  c(
    CRAN = paste0("https://p3m.dev/cran/__linux__/", codename, "/", snapshot_date),
    BioCsoft = paste0("https://bioconductor.org/packages/", bioc_version, "/bioc"),
    BioCann = paste0("https://bioconductor.org/packages/", bioc_version, "/data/annotation"),
    BioCexp = paste0("https://bioconductor.org/packages/", bioc_version, "/data/experiment")
  )
}

#' The repository pin block written into `.Rprofile`
#' @param snapshot_date CRAN snapshot date
#' @param bioc_version Bioconductor release
#' @param codename Ubuntu codename for binary packages
#' @keywords internal
pin_block <- function(snapshot_date, bioc_version, codename = "noble") {
  c(
    "",
    "# Date-tracked package repositories (managed by rmyknife::project_init).",
    "# CRAN is pinned to a Posit Package Manager snapshot; Bioconductor to a release.",
    "# Bump these two values and re-run renv::snapshot() to move the project forward in time.",
    "local({",
    sprintf('  snapshot_date <- "%s"', snapshot_date),
    sprintf('  bioc_version <- "%s"', bioc_version),
    "  options(repos = c(",
    sprintf('    CRAN     = paste0("https://p3m.dev/cran/__linux__/%s/", snapshot_date),', codename),
    '    BioCsoft = paste0("https://bioconductor.org/packages/", bioc_version, "/bioc"),',
    '    BioCann  = paste0("https://bioconductor.org/packages/", bioc_version, "/data/annotation"),',
    '    BioCexp  = paste0("https://bioconductor.org/packages/", bioc_version, "/data/experiment")',
    "  ))",
    "})",
    ""
  )
}

#' Full `.Rprofile` template for a project that does not have one yet
#' @param snapshot_date CRAN snapshot date
#' @param bioc_version Bioconductor release
#' @param codename Ubuntu codename for binary packages
#' @keywords internal
rprofile_template <- function(snapshot_date, bioc_version, codename = "noble") {
  c(
    "# A project .Rprofile shadows ~/.Rprofile completely, so source it explicitly.",
    "# Must run BEFORE renv activation: ~/.Rprofile calls .libPaths(R_LIBS_USER),",
    "# which would otherwise displace the project library.",
    'if (file.exists("~/.Rprofile")) {',
    '  source("~/.Rprofile")',
    "}",
    pin_block(snapshot_date, bioc_version, codename),
    'source("renv/activate.R")',
    "",
    "# Attach to vscode-R. Interactive only, and must come last so it sees the",
    "# final library paths. Setting TERM_PROGRAM is needed inside tmux sessions.",
    "local({",
    '  if (!interactive() || Sys.getenv("RSTUDIO") != "") {',
    "    return(invisible(NULL))",
    "  }",
    '  Sys.setenv(TERM_PROGRAM = "vscode")',
    '  init <- file.path(Sys.getenv("HOME"), ".vscode-R", "init.R")',
    "  if (file.exists(init)) {",
    "    source(init)",
    "  }",
    "})"
  )
}

#' `.Rprofile` for a subdirectory, delegating to the project root
#' @keywords internal
rprofile_subdir_template <- function() {
  c(
    "# The renv project lives one level up. Without this, starting R directly in",
    "# this directory (e.g. to run the targets pipeline) silently falls back to the",
    "# shared user library instead of the project library.",
    "# Delegate to the root profile rather than duplicating it, so the repo pins and",
    "# ~/.Rprofile are applied here too. Skipped when renv is already active, which",
    "# is the case for the targets subprocess launched by the makefile.",
    'if (Sys.getenv("RENV_PROJECT") == "") {',
    '  owd <- setwd("..")',
    '  source(".Rprofile")',
    "  setwd(owd)",
    "}"
  )
}

#' Does an .Rprofile already pin dated repositories?
#' @param lines Lines of the .Rprofile
#' @keywords internal
has_repo_pins <- function(lines) {
  any(grepl("p3m\\.dev/cran/.*[0-9]{4}-[0-9]{2}-[0-9]{2}", lines)) ||
    any(grepl("snapshot_date *<-", lines))
}

#' Does an .Rprofile set repositories in a way that would fight the pin block?
#' @param lines Lines of the .Rprofile
#' @keywords internal
has_conflicting_repos <- function(lines) {
  !has_repo_pins(lines) && any(grepl("options\\s*\\(\\s*repos", lines))
}

#' Insert the pin block above the renv activation call
#'
#' Ordering matters: the pins must be set before `renv/activate.R` runs, so a
#' project whose profile already activates renv gets the block spliced in above
#' that line rather than appended.
#'
#' @param lines Lines of the existing .Rprofile
#' @param block Pin block to insert
#' @keywords internal
insert_pin_block <- function(lines, block) {
  activate <- grep('source\\("renv/activate\\.R"\\)', lines)
  at <- if (length(activate) > 0) activate[1] - 1L else length(lines)
  append(lines, block, after = at)
}

#' Read version and repository pins out of an existing renv lockfile
#' @param path Project root
#' @keywords internal
lockfile_info <- function(path = ".") {
  lockfile <- file.path(path, "renv.lock")
  if (!file.exists(lockfile)) {
    return(NULL)
  }
  lock <- jsonlite::fromJSON(lockfile, simplifyVector = FALSE)
  cran <- NULL
  for (repo in lock$R$Repositories) {
    if (identical(repo$Name, "CRAN")) {
      cran <- repo$URL
    }
  }
  snapshot_date <- NA_character_
  if (!is.null(cran)) {
    match <- regmatches(cran, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", cran))
    if (length(match) > 0) {
      snapshot_date <- match
    }
  }
  list(
    r_version = lock$R$Version,
    bioc_version = if (is.null(lock$Bioconductor$Version)) NA_character_ else lock$Bioconductor$Version,
    snapshot_date = snapshot_date
  )
}

#' Append missing lines to a file, creating it when absent
#' @param file File to extend
#' @param lines Lines that must be present
#' @keywords internal
ensure_lines <- function(file, lines) {
  existing <- if (file.exists(file)) readLines(file, warn = FALSE) else character()
  missing <- setdiff(lines, existing)
  if (length(missing) == 0) {
    return(FALSE)
  }
  writeLines(c(existing, missing), file)
  TRUE
}

#' Stub makefile for a project that does not have one
#' @param image Path to the primary image
#' @param bind Bind mounts
#' @keywords internal
makefile_template <- function(image, bind) {
  c(
    "SHELL=/bin/bash",
    "current_date := $(shell date +'%Y-%m-%d_%H-%M')",
    sprintf("SINGULARITY=singularity exec --bind %s %s", bind, image),
    "",
    ".PHONY: verify",
    "verify:",
    "\t$(SINGULARITY) Rscript -e 'rmyknife::project_verify(strict = TRUE)'",
    "",
    ".PHONY: backup",
    "backup:",
    "\t$(SINGULARITY) Rscript -e 'rmyknife::project_backup()'"
  )
}

#' Set up a reproducible environment for an analysis project
#'
#' Records the singularity image the project runs in, pins CRAN and Bioconductor
#' to dated snapshots, and initialises renv. Safe to run on an existing project:
#' nothing that is already present and correct is overwritten.
#'
#' Must be run inside the image it records, because [renv::init()] builds the
#' package library with the running R. Recording one image from a session inside
#' a different one would produce a library that matches neither.
#'
#' Auxiliary images (tool containers such as ggsashimi) are recorded and
#' checksummed for archival but never entered.
#'
#' @param path Project root. Defaults to the working directory.
#' @param image Path to the primary singularity image. Defaults to the image the
#'   current session runs in.
#' @param aux_images Character vector of additional images used by the project.
#' @param snapshot_date CRAN snapshot date. Defaults to the date in an existing
#'   `renv.lock`, else the image build date.
#' @param bioc_version Bioconductor release. Defaults to the version in an
#'   existing `renv.lock`, else the running BiocManager version.
#' @param bind Bind mounts for the singularity invocation.
#' @param dirs Directories to create if missing.
#' @param checksum Whether to compute image checksums. Several GB per image over
#'   a network filesystem, so roughly a minute each, once.
#' @param overwrite Whether to rewrite `environment.lock` if it already exists.
#' @return The environment lock, invisibly.
#' @export
#' @examples
#' \dontrun{
#'   # New project, inside the container:
#'   project_init()
#'
#'   # Converting a project that uses a second tool container:
#'   project_init(aux_images = "/cephfs/.../ggsashimi_latest.sif")
#' }
project_init <- function(
  path = ".",
  image = Sys.getenv("APPTAINER_CONTAINER"),
  aux_images = character(),
  snapshot_date = NULL,
  bioc_version = NULL,
  bind = "/cephfs:/cephfs",
  dirs = c("analysis", "docs", "release"),
  checksum = TRUE,
  overwrite = FALSE
) {
  if (!in_container()) {
    stop(
      "project_init() must run inside the singularity image it records, because ",
      "renv builds the library with the running R.\n",
      "Start R with: singularity exec --bind ", bind, " <image> R"
    )
  }

  # Being inside a container and knowing which image file it came from are two
  # separate questions. APPTAINER_CONTAINER answers the second, and some
  # launchers (a tmux server started outside the container, nested shells,
  # --cleanenv) drop it while the session is genuinely still inside.
  running <- Sys.getenv("APPTAINER_CONTAINER")
  if (!nzchar(image)) {
    if (!nzchar(running)) {
      stop(
        "Running inside a container, but APPTAINER_CONTAINER is not set, so the ",
        "image file cannot be detected.\nPass it explicitly, for example:\n",
        '  project_init(image = "', Sys.getenv("SINGULARITY_IMAGES", "<image-dir>"),
        '/latest/mytidyverse.simg")'
      )
    }
    image <- running
  }
  if (nzchar(running)) {
    if (!identical(normalizePath(image, mustWork = FALSE), normalizePath(running, mustWork = FALSE))) {
      stop(
        "`image` is not the image this session runs in:\n",
        "  requested: ", image, "\n",
        "  running:   ", running, "\n",
        "Re-run inside the requested image; a library built by a different R would match neither."
      )
    }
  } else {
    # No path to compare against, so confirm identity by build date instead.
    running_build <- label_value(image_labels_self(), "org.label-schema.build-date")
    given_build <- label_value(image_labels_inspect(image), "org.label-schema.build-date")
    if (!is.na(running_build) && !is.na(given_build) && !identical(running_build, given_build)) {
      stop(
        "`image` is not the container this session runs in:\n",
        "  requested was built ", given_build, "\n",
        "  running was built   ", running_build
      )
    }
    if (is.na(given_build)) {
      message("⚠️  Could not confirm that ", basename(image), " is this session's container.")
    }
  }

  path <- normalizePath(path, mustWork = TRUE)
  labels <- image_labels_self()
  codename <- base_codename(labels)
  lock <- lockfile_info(path)

  if (is.null(snapshot_date)) {
    snapshot_date <- if (!is.null(lock) && !is.na(lock$snapshot_date)) {
      lock$snapshot_date
    } else {
      parse_build_date(label_value(labels, "org.label-schema.build-date"))
    }
  }
  if (is.na(snapshot_date)) {
    stop("Could not determine `snapshot_date` from renv.lock or the image labels; pass it explicitly.")
  }
  if (is.null(bioc_version)) {
    bioc_version <- if (!is.null(lock) && !is.na(lock$bioc_version)) {
      lock$bioc_version
    } else {
      as.character(BiocManager::version())
    }
  }

  message("Project:       ", basename(path))
  message("Image:         ", basename(image))
  message("R version:     ", getRversion())
  message("Snapshot date: ", snapshot_date)
  message("Bioconductor:  ", bioc_version)

  # Pin the running session too, so renv records the pinned repositories in the
  # lockfile rather than whatever was active before.
  old_repos <- getOption("repos")
  on.exit(options(repos = old_repos), add = TRUE)
  options(repos = pinned_repos(snapshot_date, bioc_version, codename))

  rprofile <- file.path(path, ".Rprofile")
  if (!file.exists(rprofile)) {
    writeLines(rprofile_template(snapshot_date, bioc_version, codename), rprofile)
    message("✅ Wrote .Rprofile")
  }

  if (!file.exists(file.path(path, "renv.lock"))) {
    message("Initialising renv, this takes a while")
    renv::init(project = path, bioconductor = bioc_version, restart = FALSE)
    renv::snapshot(project = path, prompt = FALSE)
    message("✅ Initialised renv")
  } else {
    message("✅ renv.lock present, left untouched")
  }

  # renv::init() creates or extends .Rprofile itself, so the pin check happens
  # afterwards. An existing profile keeps everything it had; only the pins are
  # guaranteed. Without them the project would snapshot against a moving CRAN.
  lines <- readLines(rprofile, warn = FALSE)
  if (has_conflicting_repos(lines)) {
    stop(
      ".Rprofile already sets options(repos = ...) to something that is not a dated ",
      "snapshot. Adding the pin block would leave two conflicting settings.\n",
      "Resolve by hand, then re-run project_init()."
    )
  }
  if (!has_repo_pins(lines)) {
    writeLines(insert_pin_block(lines, pin_block(snapshot_date, bioc_version, codename)), rprofile)
    message("✅ Inserted repository pins into existing .Rprofile")
  } else {
    message("✅ .Rprofile already pins dated repositories")
  }

  for (dir in dirs) {
    target <- file.path(path, dir)
    if (!dir.exists(target)) {
      dir.create(target, recursive = TRUE)
      message("✅ Created ", dir, "/")
    }
  }

  analysis_rprofile <- file.path(path, "analysis", ".Rprofile")
  if (dir.exists(dirname(analysis_rprofile)) && !file.exists(analysis_rprofile)) {
    writeLines(rprofile_subdir_template(), analysis_rprofile)
    message("✅ Wrote analysis/.Rprofile")
  }

  images <- list(describe_image(image, "primary", labels = labels, checksum = checksum))
  for (aux in aux_images) {
    images[[length(images) + 1L]] <- describe_image(aux, "aux", checksum = checksum)
  }

  env_lock <- list(
    project = basename(path),
    rmyknife_version = as.character(utils::packageVersion("rmyknife")),
    created = format(Sys.Date()),
    r_version = as.character(getRversion()),
    snapshot_date = snapshot_date,
    bioc_version = bioc_version,
    bind = bind,
    images = images
  )

  lock_path <- env_lock_path(path)
  if (file.exists(lock_path) && !overwrite) {
    message("⚠️  ", env_lock_file, " exists, not overwritten. Use overwrite = TRUE to refresh.")
  } else {
    write_env_lock(env_lock, lock_path)
    message("✅ Wrote ", env_lock_file)
  }

  makefile <- file.path(path, "makefile")
  singularity_line <- sprintf("SINGULARITY=singularity exec --bind %s %s", bind, image)
  if (!file.exists(makefile) && !file.exists(file.path(path, "Makefile"))) {
    writeLines(makefile_template(image, bind), makefile)
    message("✅ Wrote makefile")
  } else {
    message("⚠️  makefile exists, left untouched. It should contain:\n    ", singularity_line)
  }

  ensure_lines(
    file.path(path, ".gitignore"),
    c("_targets", "release", "analysis/*html", "backup")
  )

  invisible(env_lock)
}

#' Write the environment lock
#' @param env_lock Environment lock list
#' @param path Destination file
#' @keywords internal
write_env_lock <- function(env_lock, path) {
  writeLines(jsonlite::toJSON(env_lock, auto_unbox = TRUE, pretty = TRUE, null = "null"), path)
}

#' Read the environment lock of a project
#' @param path Project root
#' @return The environment lock as a list
#' @export
read_env_lock <- function(path = ".") {
  lock_path <- env_lock_path(path)
  if (!file.exists(lock_path)) {
    stop("No ", env_lock_file, " in ", normalizePath(path, mustWork = FALSE), ". Run project_init() first.")
  }
  jsonlite::fromJSON(lock_path, simplifyVector = FALSE)
}
