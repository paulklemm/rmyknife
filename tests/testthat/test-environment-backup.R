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
  expect_true(file.exists(file.path(destination, ".Rprofile")))
  expect_true(file.exists(file.path(destination, "renv", "settings.json")))
  expect_true(file.exists(file.path(destination, "images", basename(fixture$image))))
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
  env_lock$images[[1]]$docker <- "paulklemm/mytidyverse:4.6.1-1"
  write_env_lock(env_lock, env_lock_path(fixture$project))

  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  instructions <- readLines(file.path(destination, "RESTORE.md"))
  expect_true(any(grepl("paulklemm/mytidyverse:4.6.1-1", instructions, fixed = TRUE)))
})

test_that("git-lfs projects are told their blobs are not in the bundle", {
  fixture <- fake_project()
  expect_false(uses_git_lfs(fixture$project))

  writeLines(
    "docs/slides.key filter=lfs diff=lfs merge=lfs -text",
    file.path(fixture$project, ".gitattributes")
  )
  expect_true(uses_git_lfs(fixture$project))

  archive <- suppressMessages(project_backup(fixture$project, include = "images"))
  destination <- file.path(fixture$root, "restored")
  suppressMessages(project_restore(archive, destination))

  instructions <- readLines(file.path(destination, "RESTORE.md"))
  expect_true(any(grepl("git-lfs", instructions, fixed = TRUE)))
  expect_true(any(grepl("GIT_LFS_SKIP_SMUDGE", instructions, fixed = TRUE)))
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
