test_that("a consistent project passes the checks that do not need renv", {
  fixture <- fake_project()
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))

  expect_equal(status_of(report, "environment.lock"), "ok")
  expect_equal(status_of(report, "project files"), "ok")
  expect_equal(status_of(report, "image:analysis-1.2.3.simg"), "ok")
  expect_equal(status_of(report, "lockfile consistency"), "ok")
  expect_equal(status_of(report, ".Rprofile pins"), "ok")
  expect_equal(status_of(report, "restorable from lockfile"), "ok")
})

test_that("a missing image fails and names the docker rebuild route", {
  fixture <- fake_project()
  file.remove(fixture$image)
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))

  row <- report[report$check == "image:analysis-1.2.3.simg", ]
  expect_equal(row$status, "fail")
  expect_match(row$detail, "unrecoverable")
})

test_that("a missing image reports its docker tag when one was recorded", {
  fixture <- fake_project()
  env_lock <- read_env_lock(fixture$project)
  env_lock$images[[1]]$docker <- "example/analysis:1.2.3"
  write_env_lock(env_lock, env_lock_path(fixture$project))
  file.remove(fixture$image)

  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_match(report$detail[report$check == "image:analysis-1.2.3.simg"], "example/analysis")
})

test_that("a tampered image is caught by size, and by checksum when unchanged in size", {
  fixture <- fake_project()
  original <- readLines(fixture$image)

  writeLines(paste0(original, " plus more"), fixture$image)
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "image:analysis-1.2.3.simg"), "fail")

  # Same byte count, different content: only a deep check can see this.
  writeLines(paste0(substr(original, 1, nchar(original) - 1), "X"), fixture$image)
  shallow <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(shallow, "image:analysis-1.2.3.simg"), "ok")

  deep <- suppressMessages(project_verify(fixture$project, network = FALSE, deep = TRUE))
  expect_equal(status_of(deep, "image:analysis-1.2.3.simg"), "fail")
})

test_that("bumping the .Rprofile date without re-snapshotting is caught", {
  fixture <- fake_project()
  writeLines(
    rprofile_template("2026-09-01", "3.23", "noble"),
    file.path(fixture$project, ".Rprofile")
  )
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, ".Rprofile pins"), "fail")
})

test_that("a lockfile disagreeing with the environment lock is caught", {
  fixture <- fake_project()
  fake_renv_lock(fixture$project, r_version = "4.5.1", snapshot_date = "2025-11-05", bioc_version = "3.21")
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))

  expect_equal(status_of(report, "lockfile consistency"), "fail")
  expect_match(report$detail[report$check == "lockfile consistency"], "4\\.5\\.1")
})

test_that("the two restorability tiers are reported independently", {
  fixture <- fake_project(packages = list(
    dplyr = list(Package = "dplyr", Version = "1.2.1", Source = "Repository", Repository = "CRAN"),
    mystery = list(Package = "mystery", Version = "0.1", Source = "unknown")
  ))
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))

  # The image is present, so the project can be backed up right now, even though
  # a package will never restore over the network. This is the normal state of a
  # freshly converted pre-renv project.
  expect_equal(status_of(report, "restorable from backup"), "ok")
  expect_equal(status_of(report, "restorable from lockfile"), "warn")
  expect_match(report$detail[report$check == "restorable from lockfile"], "mystery")
})

test_that("missing project files are reported", {
  fixture <- fake_project()
  file.remove(file.path(fixture$project, "renv", "settings.json"))
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))

  expect_equal(status_of(report, "project files"), "fail")
  expect_match(report$detail[report$check == "project files"], "settings.json")
})

test_that("a symlinked container resolves to the versioned image it points at", {
  fixture <- fake_project()
  symlink <- file.path(dirname(fixture$image), "analysis.simg")
  file.symlink(fixture$image, symlink)

  # init records what the symlink points at, never the moving pointer itself.
  described <- describe_image(symlink, "primary", labels = list(), checksum = TRUE)
  expect_equal(described$path, fixture$image)
  expect_equal(described$name, "analysis-1.2.3.simg")

  # A session started through the symlink is still the recorded container.
  withr::local_envvar(APPTAINER_CONTAINER = symlink)
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "running image"), "ok")
  expect_match(report$detail[report$check == "running image"], "via")
})

test_that("a symlink that has moved on to another image is caught", {
  fixture <- fake_project()
  other <- fake_image(dirname(fixture$image), "analysis-2.0.0.simg", "a newer image")
  symlink <- file.path(dirname(fixture$image), "analysis.simg")
  file.symlink(other, symlink)

  withr::local_envvar(APPTAINER_CONTAINER = symlink)
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "running image"), "warn")
  expect_match(report$detail[report$check == "running image"], "analysis-2\\.0\\.0")
})

test_that("without APPTAINER_CONTAINER the running image is checked by build date", {
  skip_if_not(in_container(), "not running inside a container")
  fixture <- fake_project()
  running_build <- parse_build_date(
    label_value(image_labels_self(), "org.label-schema.build-date")
  )
  skip_if(is.na(running_build), "container has no build-date label")

  env_lock <- read_env_lock(fixture$project)
  env_lock$images[[1]]$build_date <- running_build
  write_env_lock(env_lock, env_lock_path(fixture$project))

  withr::local_envvar(APPTAINER_CONTAINER = "")
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "running image"), "ok")
  expect_match(report$detail[report$check == "running image"], "build date")
})

test_that("a build date that disagrees with the record is caught", {
  skip_if_not(in_container(), "not running inside a container")
  fixture <- fake_project()
  env_lock <- read_env_lock(fixture$project)
  env_lock$images[[1]]$build_date <- "1999-01-01"
  write_env_lock(env_lock, env_lock_path(fixture$project))

  withr::local_envvar(APPTAINER_CONTAINER = "")
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "running image"), "warn")
  expect_match(report$detail[report$check == "running image"], "1999-01-01")
})

test_that("a makefile naming the image through a symlink still matches", {
  fixture <- fake_project()
  symlink <- file.path(dirname(fixture$image), "analysis.simg")
  file.symlink(fixture$image, symlink)
  writeLines(
    paste0("SINGULARITY=singularity exec --bind /data:/data ", symlink),
    file.path(fixture$project, "makefile")
  )

  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "makefile image"), "ok")
})

test_that("a makefile building the image path from a variable says so", {
  fixture <- fake_project()
  writeLines(
    "SINGULARITY=singularity exec --bind /data:/data $(IMAGE_DIR)/latest/analysis.simg",
    file.path(fixture$project, "makefile")
  )

  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "makefile image"), "warn")
  expect_match(report$detail[report$check == "makefile image"], "variable")
})

test_that("a makefile pointing at a different image is flagged", {
  fixture <- fake_project()
  writeLines(
    "SINGULARITY=singularity exec --bind /data:/data /somewhere/else.simg",
    file.path(fixture$project, "makefile")
  )
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "makefile image"), "warn")
})

test_that("strict mode raises exactly when a check fails", {
  fixture <- fake_project()
  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  # The fabricated project has a lockfile but no installed library, so renv
  # rightly reports it out of sync. Everything the container layer owns passes.
  expect_false(any(report$status[report$check != "renv status"] == "fail"))

  file.remove(fixture$image)
  expect_error(
    suppressMessages(project_verify(fixture$project, network = FALSE, strict = TRUE)),
    "Verification failed"
  )
})

test_that("verify fails cleanly on an uninitialised project", {
  path <- withr::local_tempdir()
  report <- suppressMessages(project_verify(path, network = FALSE))
  expect_equal(report$status, "fail")
  expect_equal(nrow(report), 1)
})

test_that("a missing .Rprofile is reported, not raised", {
  fixture <- fake_project()
  file.remove(file.path(fixture$project, ".Rprofile"))

  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "project files"), "fail")
  expect_equal(status_of(report, ".Rprofile pins"), "fail")
  expect_match(report$detail[report$check == ".Rprofile pins"], "no .Rprofile")
})

test_that("a lock without a primary image is reported, not raised", {
  fixture <- fake_project()
  env_lock <- read_env_lock(fixture$project)
  env_lock$images[[1]]$role <- "aux"
  write_env_lock(env_lock, env_lock_path(fixture$project))

  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "environment.lock"), "fail")
  expect_match(report$detail[report$check == "environment.lock"], "primary")
})

test_that("a makefile naming the image through one whole variable says so", {
  # The `.simg` literal is inside the variable here, so there is no path to
  # extract and the variable has to be spotted on the line itself.
  fixture <- fake_project()
  writeLines(
    "SINGULARITY=singularity exec --bind /data:/data $(IMAGE)",
    file.path(fixture$project, "makefile")
  )

  report <- suppressMessages(project_verify(fixture$project, network = FALSE))
  expect_equal(status_of(report, "makefile image"), "warn")
  expect_match(report$detail[report$check == "makefile image"], "variable")
})
