test_that("build-date labels parse regardless of locale", {
  expect_equal(parse_build_date("Tuesday_21_July_2026_10:1:14_UTC"), "2026-07-21")
  expect_equal(parse_build_date("Wednesday_5_November_2025_09:0:00_UTC"), "2025-11-05")
  expect_true(is.na(parse_build_date(NA_character_)))
  expect_true(is.na(parse_build_date("nonsense")))
})

test_that("the ubuntu codename comes from the image base label", {
  expect_equal(base_codename(list("org.opencontainers.image.base.name" = "docker.io/library/ubuntu:noble")), "noble")
  expect_equal(base_codename(list("org.opencontainers.image.base.name" = "ubuntu:jammy")), "jammy")
  expect_equal(base_codename(list()), "noble")
})

test_that("pinned repositories are dated and Bioconductor is versioned", {
  repos <- pinned_repos("2026-07-24", "3.23", "noble")
  expect_equal(unname(repos[["CRAN"]]), "https://p3m.dev/cran/__linux__/noble/2026-07-24")
  expect_match(repos[["BioCsoft"]], "packages/3.23/bioc")
})

test_that("a generated .Rprofile pins repositories before activating renv", {
  lines <- rprofile_template("2026-07-24", "3.23", "noble")
  expect_true(has_repo_pins(lines))
  pins <- grep("snapshot_date <-", lines)
  activate <- grep('source\\("renv/activate\\.R"\\)', lines)
  expect_length(activate, 1)
  expect_lt(pins, activate)
})

test_that("existing repository pins are detected and conflicting ones flagged", {
  pinned <- c('options(repos = c(CRAN = "https://p3m.dev/cran/2026-07-24"))')
  expect_true(has_repo_pins(pinned))
  expect_false(has_conflicting_repos(pinned))

  moving <- c('options(repos = c(CRAN = "https://cloud.r-project.org"))')
  expect_false(has_repo_pins(moving))
  expect_true(has_conflicting_repos(moving))

  silent <- c("# nothing to see here")
  expect_false(has_repo_pins(silent))
  expect_false(has_conflicting_repos(silent))
})

test_that("pins are spliced in above renv activation, preserving the rest", {
  existing <- c(
    'source("~/.Rprofile")',
    "# a comment the user wrote",
    'source("renv/activate.R")',
    "# a trailing comment"
  )
  result <- insert_pin_block(existing, pin_block("2026-07-24", "3.23", "noble"))

  expect_true(all(existing %in% result))
  pins <- grep("snapshot_date <-", result)
  activate <- grep('source\\("renv/activate\\.R"\\)', result)
  expect_lt(pins, activate)
  expect_equal(result[length(result)], "# a trailing comment")
})

test_that("pins are appended when there is no renv activation yet", {
  existing <- c("# just a profile")
  result <- insert_pin_block(existing, pin_block("2026-07-24", "3.23", "noble"))
  expect_true(has_repo_pins(result))
  expect_equal(result[1], "# just a profile")
})

test_that("status tags are fixed width so consecutive lines align", {
  tags <- vapply(c("ok", "warn", "fail", "info"), status_tag, character(1))
  expect_equal(length(unique(nchar(tags))), 1)
  expect_false(any(grepl("[^[:print:]]", tags)))
  expect_equal(status_tag("fail"), "[FAIL]")

  expect_message(status_message("warn", "something"), "\\[WARN\\] something")
})

test_that("bind mounts are site specific and omitted unless asked for", {
  expect_equal(singularity_exec("/img.simg"), "singularity exec /img.simg")
  expect_equal(singularity_exec("/img.simg", NULL), "singularity exec /img.simg")
  expect_equal(singularity_exec("/img.simg", ""), "singularity exec /img.simg")
  expect_equal(
    singularity_exec("/img.simg", "/data:/data"),
    "singularity exec --bind /data:/data /img.simg"
  )

  plain <- makefile_template("/img.simg")
  expect_true(any(grepl("^SINGULARITY=singularity exec /img.simg$", plain)))
  expect_false(any(grepl("--bind", plain)))

  bound <- makefile_template("/img.simg", "/data:/data")
  expect_true(any(grepl("--bind /data:/data", bound, fixed = TRUE)))
})

test_that("lockfile pins are read back out of renv.lock", {
  path <- withr::local_tempdir()
  fake_renv_lock(path, r_version = "4.5.1", snapshot_date = "2025-11-05", bioc_version = "3.21")
  info <- lockfile_info(path)
  expect_equal(info$r_version, "4.5.1")
  expect_equal(info$snapshot_date, "2025-11-05")
  expect_equal(info$bioc_version, "3.21")
})

test_that("lockfile_info returns NULL without a lockfile", {
  expect_null(lockfile_info(withr::local_tempdir()))
})

test_that("ensure_lines appends only what is missing", {
  file <- file.path(withr::local_tempdir(), ".gitignore")
  writeLines(c("release"), file)

  expect_true(ensure_lines(file, c("release", "_targets")))
  expect_equal(readLines(file), c("release", "_targets"))
  expect_false(ensure_lines(file, c("release", "_targets")))
})

test_that("unresolvable packages are exactly those a network restore would miss", {
  path <- withr::local_tempdir()
  fake_renv_lock(path, packages = list(
    dplyr = list(Package = "dplyr", Version = "1.2.1", Source = "Repository", Repository = "CRAN"),
    mystery = list(Package = "mystery", Version = "0.1", Source = "unknown"),
    homemade = list(Package = "homemade", Version = "0.1", Source = "Local"),
    pinned = list(Package = "pinned", Version = "0.1", Source = "GitHub", RemoteType = "github", RemoteSha = "abc123"),
    floating = list(Package = "floating", Version = "0.1", Source = "GitHub", RemoteType = "github")
  ))
  expect_setequal(
    unresolvable_packages(file.path(path, "renv.lock")),
    c("mystery", "homemade", "floating")
  )
})

test_that("being in a container is detected by its metadata, not an env var", {
  # The env var says which image; the marker directory says whether we are in
  # one at all. Conflating them made project_init() refuse to run in sessions
  # that had lost APPTAINER_CONTAINER but were genuinely inside the container.
  expect_false(in_container(file.path(withr::local_tempdir(), "absent")))

  marker <- file.path(withr::local_tempdir(), ".singularity.d")
  dir.create(marker)
  expect_true(in_container(marker))

  withr::local_envvar(APPTAINER_CONTAINER = "")
  expect_true(in_container(marker))
})

test_that("project_init asks for the image when the container did not name it", {
  # This test runs inside a container, so the marker is present but the env var
  # is not: exactly the situation a tmux session started outside produces.
  skip_if_not(in_container(), "not running inside a container")
  withr::local_envvar(APPTAINER_CONTAINER = "")
  expect_error(
    project_init(withr::local_tempdir()),
    "APPTAINER_CONTAINER is not set"
  )
})

test_that("project_init refuses to record an image it is not running in", {
  root <- withr::local_tempdir()
  running <- fake_image(root, "running.simg")
  other <- fake_image(root, "other.simg")
  withr::local_envvar(APPTAINER_CONTAINER = running)
  expect_error(project_init(root, image = other), "not the image this session runs in")
})

test_that("read_env_lock explains itself when the project was never initialised", {
  expect_error(read_env_lock(withr::local_tempdir()), "Run project_init")
})

test_that("project_init rejects a missing image before renv::init() runs", {
  root <- withr::local_tempdir()
  running <- fake_image(root, "running.simg")
  withr::local_envvar(APPTAINER_CONTAINER = running)
  expect_error(
    project_init(root, aux_images = file.path(root, "absent.simg")),
    "Image not found"
  )
  # renv::init() is the expensive step; nothing should have been built yet.
  expect_false(file.exists(file.path(root, "renv.lock")))
})

test_that("a lock field that was NA reads back as absent, not as NULL", {
  # NA round-trips through JSON as null, so absence has two spellings on read.
  expect_equal(lock_field(NULL, "fallback"), "fallback")
  expect_equal(lock_field(NA_character_, "fallback"), "fallback")
  expect_equal(lock_field(character(), "fallback"), "fallback")
  expect_equal(lock_field("value", "fallback"), "value")
})

test_that("a lock without a primary image is refused by name", {
  fixture <- fake_project()
  env_lock <- read_env_lock(fixture$project)
  env_lock$images[[1]]$role <- "aux"
  expect_error(primary_image(env_lock), "No image with role")
})

test_that("image labels are read from the container marker directory", {
  marker <- withr::local_tempdir()
  expect_equal(image_labels_self(marker), list())
  writeLines('{"org.label-schema.build-date": "Tuesday_21_July_2026_10:1:14_UTC"}',
             file.path(marker, "labels.json"))
  expect_equal(
    label_value(image_labels_self(marker), "org.label-schema.build-date"),
    "Tuesday_21_July_2026_10:1:14_UTC"
  )
})
