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

test_that("a library outside the project is refused, not silently misrooted", {
  # tar strips a leading `/` with a warning and exit status 0, so an absolute
  # member would restore under destination/<abs path> where renv never looks.
  fixture <- fake_project()
  withr::local_envvar(RENV_PATHS_LIBRARY = file.path(fixture$root, "external-lib"))
  dir.create(renv::paths$library(project = fixture$project), recursive = TRUE)

  expect_error(active_library(fixture$project), "outside the project")
  expect_error(
    suppressMessages(project_backup(fixture$project)),
    "outside the project"
  )
})

test_that("a failed restore does not block the retry", {
  fixture <- fake_project()
  corrupt <- file.path(fixture$root, "corrupt.tar.zst")
  writeLines("not an archive", corrupt)
  destination <- file.path(fixture$root, "restored")

  expect_error(suppressMessages(project_restore(corrupt, destination)), "failed with exit status")
  expect_false(dir.exists(destination))

  # The retry now reaches the archive rather than tripping over the leftovers.
  expect_error(suppressMessages(project_restore(corrupt, destination)), "failed with exit status")
})

test_that("a successful restore keeps its destination", {
  fixture <- fake_project()
  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  expect_true(dir.exists(destination))
  expect_true(file.exists(file.path(destination, "MANIFEST.json")))
})

test_that("a symlinked renv/library is archived, not refused", {
  # Linking renv/library at faster storage is a common habit on a cluster. The
  # member stays project-relative and tar's -h dereferences it, so this is
  # archivable and must not be mistaken for a library outside the project.
  fixture <- fake_project()
  version <- paste0("R-", getRversion()$major, ".", getRversion()$minor)
  elsewhere <- file.path(fixture$root, "fast-storage")
  dir.create(file.path(elsewhere, "linux-ubuntu-noble", version, "x86_64-pc-linux-gnu", "dplyr"), recursive = TRUE)
  writeLines("Package: dplyr", file.path(elsewhere, "linux-ubuntu-noble", version, "x86_64-pc-linux-gnu", "dplyr", "DESCRIPTION"))
  dir.create(file.path(fixture$project, "renv"), showWarnings = FALSE)
  file.symlink(elsewhere, file.path(fixture$project, "renv", "library"))

  expect_equal(
    active_library(fixture$project),
    file.path("renv", "library", "linux-ubuntu-noble", version)
  )

  # And it survives the round trip, dereferenced.
  archive <- suppressMessages(project_backup(fixture$project, include = "library"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))
  expect_true(file.exists(file.path(
    destination, "renv", "library", "linux-ubuntu-noble", version,
    "x86_64-pc-linux-gnu", "dplyr", "DESCRIPTION"
  )))
})

test_that("staging left by a killed run does not accumulate in the project", {
  # Staging lives beside the archive rather than in tempdir(), so a run that
  # died without unwinding leaves gigabytes in the project until swept.
  fixture <- fake_project()
  stale <- file.path(fixture$project, "backup", "project_2020-01-01_deadbee-staging")
  dir.create(stale, recursive = TRUE)
  writeLines("a large library tarball", file.path(stale, "library.tar.zst"))

  # Something the caller owns, which a bare *-staging sweep would have eaten.
  bystander <- file.path(fixture$project, "backup", "my-own-staging")
  dir.create(bystander)

  archive <- suppressMessages(project_backup(fixture$project, include = "images"))

  expect_false(dir.exists(stale))
  expect_true(dir.exists(bystander))
  expect_true(file.exists(archive))
})
