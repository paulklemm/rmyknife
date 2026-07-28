test_that("a backup round-trips through restore with checksums intact", {
  fixture <- fake_project()
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))

  expect_true(file.exists(archive))
  expect_match(basename(archive), "^project_\\d{4}-\\d{2}-\\d{2}_nogit\\.tar\\.zst$")

  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  expect_true(file.exists(file.path(destination, "MANIFEST.json")))
  expect_true(file.exists(file.path(destination, "CHECKSUMS.sha256")))
  expect_true(file.exists(file.path(destination, "RESTORE.md")))
  expect_true(file.exists(file.path(destination, "renv.lock")))
  expect_true(file.exists(file.path(destination, "images", basename(fixture$image))))
})

test_that("the backup carries the environment only, not the repo", {
  fixture <- fake_project()
  dir.create(file.path(fixture$project, "analysis"))
  writeLines("# a repo file", file.path(fixture$project, "analysis", "report.Rmd"))
  writeLines("# delegating profile", file.path(fixture$project, "analysis", ".Rprofile"))
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  # Everything git already holds stays in git.
  expect_false(file.exists(file.path(destination, ".Rprofile")))
  expect_false(file.exists(file.path(destination, "analysis", ".Rprofile")))
  expect_false(file.exists(file.path(destination, "analysis", "report.Rmd")))
  expect_false(file.exists(file.path(destination, "renv", "settings.json")))
  expect_false(file.exists(file.path(destination, "renv", "activate.R")))
  expect_false(file.exists(file.path(destination, "repo.bundle")))
  expect_false(file.exists(file.path(destination, "environment.lock")))

  # renv.lock stays, because it describes the archived library.
  expect_true(file.exists(file.path(destination, "renv.lock")))
})

test_that("RESTORE.md says where the missing half lives", {
  fixture <- fake_project()
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  instructions <- paste(readLines(file.path(destination, "RESTORE.md")), collapse = " ")
  expect_match(instructions, "compute environment only")
  expect_match(instructions, "live in git")
  expect_match(instructions, "Check the repository out at commit")
})

test_that("the restored image is byte-identical to the original", {
  fixture <- fake_project()
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  expect_equal(
    sha256_file(file.path(destination, "images", basename(fixture$image))),
    sha256_file(fixture$image)
  )
})

test_that("restore detects a damaged archive", {
  fixture <- fake_project()
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))

  # Extract, corrupt a member, repack: the checksums no longer describe the tree.
  scratch <- file.path(fixture$root, "scratch")
  dir.create(scratch)
  system2("tar", c("--use-compress-program=zstd", "-x", "-f", shQuote(archive), "-C", shQuote(scratch)))
  writeLines("tampered", file.path(scratch, "renv.lock"))
  damaged <- file.path(fixture$root, "damaged.tar.zst")
  system2("tar", c("-c", "--use-compress-program=zstd", "-f", shQuote(damaged), "-C", shQuote(scratch), "."))

  expect_error(
    suppressMessages(project_restore(damaged, file.path(fixture$root, "restored-damaged"))),
    "Archive is damaged"
  )
})

test_that("the manifest records what went into the backup", {
  fixture <- fake_project()
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  manifest <- jsonlite::fromJSON(file.path(destination, "MANIFEST.json"), simplifyVector = FALSE)
  expect_equal(manifest$project, "project")
  expect_equal(manifest$snapshot_date, "2026-07-24")
  expect_equal(manifest$git_revision, "nogit")
  expect_true("MANIFEST.json" %in% unlist(manifest$contents))
  expect_true(paste0("images/", basename(fixture$image)) %in% unlist(manifest$contents))
})

test_that("RESTORE.md names the docker rebuild route when one is known", {
  fixture <- fake_project()
  env_lock <- read_env_lock(fixture$project)
  env_lock$images[[1]]$docker <- "example/analysis:1.2.3"
  write_env_lock(env_lock, env_lock_path(fixture$project))

  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  instructions <- readLines(file.path(destination, "RESTORE.md"))
  expect_true(any(grepl("example/analysis:1.2.3", instructions, fixed = TRUE)))
})

test_that("backups do not silently overwrite each other", {
  fixture <- fake_project()
  suppressMessages(project_backup(fixture$project, include = "images"))
  expect_error(
    suppressMessages(project_backup(fixture$project, include = "images")),
    "Backup already exists"
  )
})

test_that("restore refuses to write into an existing directory", {
  fixture <- fake_project()
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  dir.create(destination)

  expect_error(suppressMessages(project_restore(archive, destination)), "already exists")
})

test_that("backup requires an initialised project", {
  expect_error(
    suppressMessages(project_backup(withr::local_tempdir())),
    "Run project_init"
  )
})

test_that("a backup made without recorded checksums still restores", {
  # project_init(checksum = FALSE) leaves sha256 absent. The checksum file has
  # to describe the archive regardless, or restore condemns a healthy backup.
  fixture <- fake_project()
  env_lock <- read_env_lock(fixture$project)
  env_lock$images[[1]]$sha256 <- NA_character_
  write_env_lock(env_lock, env_lock_path(fixture$project))

  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  expect_no_error(suppressMessages(project_restore(archive, destination)))

  checksums <- readLines(file.path(destination, "CHECKSUMS.sha256"))
  image_line <- grep("images/", checksums, value = TRUE, fixed = TRUE)
  expect_match(image_line, "^[0-9a-f]{64}  images/")
})

test_that("an image that drifted from its record stops the backup", {
  fixture <- fake_project()
  writeLines("tampered with since project_init", fixture$image)

  expect_error(
    suppressMessages(project_backup(fixture$project, include = "images")),
    "Image has changed since it was recorded"
  )
})

test_that("the manifest lists only the images that were really archived", {
  fixture <- fake_project()
  env_lock <- read_env_lock(fixture$project)
  env_lock$images[[2]] <- list(
    role = "aux", path = "/nowhere/tool.simg", name = "tool.simg",
    bytes = 1, sha256 = NA_character_
  )
  write_env_lock(env_lock, env_lock_path(fixture$project))

  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  contents <- unlist(jsonlite::fromJSON(
    file.path(destination, "MANIFEST.json"), simplifyVector = FALSE
  )$contents)
  expect_true(paste0("images/", basename(fixture$image)) %in% contents)
  expect_false("images/tool.simg" %in% contents)
})

test_that("restore refuses an archive with no checksum file", {
  fixture <- fake_project()
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))

  scratch <- file.path(fixture$root, "scratch")
  dir.create(scratch)
  system2("tar", c("--use-compress-program=zstd", "-x", "-f", shQuote(archive), "-C", shQuote(scratch)))
  file.remove(file.path(scratch, "CHECKSUMS.sha256"))
  stripped <- file.path(fixture$root, "stripped.tar.zst")
  system2("tar", c("-c", "--use-compress-program=zstd", "-f", shQuote(stripped), "-C", shQuote(scratch), "."))

  expect_error(
    suppressMessages(project_restore(stripped, file.path(fixture$root, "restored"))),
    "no CHECKSUMS.sha256"
  )
})

test_that("a failing archive command stops instead of reporting success", {
  expect_error(run_or_stop("tar", c("-c", "-f", shQuote("/proc/nope/x.tar"), "--", "/nonexistent")), "failed with exit status")
  expect_silent(run_or_stop("true", character()))
})

test_that("the library path survives regex metacharacters in the project path", {
  # A project directory holding a `+` used to leave the path absolute, which
  # made the library tarball step archive the wrong tree. Built where renv
  # actually puts it, so the expectation holds on any platform.
  project <- file.path(withr::local_tempdir(), "a+b(c)")
  library_path <- renv::paths$library(project = project)
  dir.create(library_path, recursive = TRUE)

  expected <- dirname(substring(normalizePath(library_path), nchar(normalizePath(project)) + 2L))
  expect_equal(active_library(project), expected)
  expect_false(startsWith(active_library(project), "/"))
})

test_that("the architecture directory is stripped whatever the architecture", {
  # A platform name renv would never compute, so the glob fallback is what runs.
  root <- withr::local_tempdir()
  version <- paste0("R-", getRversion()$major, ".", getRversion()$minor)
  relative <- file.path("renv", "library", "linux-fake-platform", version)
  dir.create(file.path(root, relative, "aarch64-unknown-linux-gnu"), recursive = TRUE)

  expect_equal(active_library(root), relative)
})
