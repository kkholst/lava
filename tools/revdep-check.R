# Reverse dependency check for the lava package
#
# Checks all CRAN reverse dependencies against the local version of lava.
# Runs R CMD check (with tests) on each revdep in parallel.
#
# Usage:
#   Rscript tools/revdep-check.R
#   Rscript tools/revdep-check.R --no-vignettes
#   REVDEP_PKGS="mets,targeted" Rscript tools/revdep-check.R
#
# Environment variables:
#   REVDEP_PKGS  - comma-separated list of packages to check (default: all)

# -- bootstrap pak -----------------------------------------------------------

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cloud.r-project.org")
}

# -- configuration -----------------------------------------------------------

pkg_name <- "lava"
pkg_dir <- normalizePath(".")
revdep_root <- file.path(pkg_dir, "revdep")
lib_dir <- file.path(revdep_root, "library")
src_dir <- file.path(revdep_root, "src")
checks_dir <- file.path(revdep_root, "checks")
num_workers <- max(1L, parallel::detectCores() - 1L)
cran_repo <- "https://cloud.r-project.org"

args <- commandArgs(trailingOnly = TRUE)
skip_vignettes <- "--no-vignettes" %in% args

# -- helpers -----------------------------------------------------------------

msg <- function(...) {
  message(sprintf(...))
}

timestamp_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
}

parse_check_log <- function(log_path) {
  if (!file.exists(log_path)) {
    return(list(status = "ERROR", warnings = 0L, notes = 0L,
                errors = 1L, details = "check log not found"))
  }
  lines <- readLines(log_path, warn = FALSE)
  status_line <- grep("^Status:", lines, value = TRUE)
  if (length(status_line) == 0L) {
    return(list(status = "ERROR", warnings = 0L, notes = 0L,
                errors = 1L, details = "no status line in log"))
  }
  status_line <- status_line[length(status_line)]
  n_err <- length(grep("ERROR", status_line))
  n_warn <- length(grep("WARNING", status_line))
  n_note <- length(grep("NOTE", status_line))

  if (n_err > 0L) {
    status <- "ERROR"
  } else if (n_warn > 0L) {
    status <- "WARNING"
  } else if (n_note > 0L) {
    status <- "NOTE"
  } else {
    status <- "OK"
  }

  # extract test failure details if any
  test_fail_idx <- grep("tests.*FAILED|test.*failure", lines, ignore.case = TRUE)
  test_details <- if (length(test_fail_idx) > 0L) {
    paste(lines[test_fail_idx], collapse = "; ")
  } else {
    ""
  }

  list(status = status, warnings = n_warn, notes = n_note,
       errors = n_err, details = test_details)
}

# -- setup directories -------------------------------------------------------

dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checks_dir, recursive = TRUE, showWarnings = FALSE)

# -- discover reverse dependencies -------------------------------------------

timestamp_msg("Querying CRAN for reverse dependencies of '%s'...", pkg_name)

available <- utils::available.packages(repos = cran_repo)

revdeps <- tools::package_dependencies(
  pkg_name,
  db = available,
  reverse = TRUE,
  which = c("Depends", "Imports", "LinkingTo", "Suggests")
)[[pkg_name]]

if (is.null(revdeps) || length(revdeps) == 0L) {
  msg("No reverse dependencies found on CRAN. Nothing to do.")
  quit(status = 0L)
}

# filter to user-specified subset if REVDEP_PKGS is set
env_pkgs <- Sys.getenv("REVDEP_PKGS", "")
if (nzchar(env_pkgs)) {
  requested <- trimws(unlist(strsplit(env_pkgs, ",")))
  missing <- setdiff(requested, revdeps)
  if (length(missing) > 0L) {
    msg("WARNING: requested packages not in revdeps: %s", paste(missing, collapse = ", "))
  }
  revdeps <- intersect(revdeps, requested)
}

# only keep packages actually on CRAN
revdeps <- intersect(revdeps, rownames(available))

msg("Found %d reverse dependencies to check.", length(revdeps))
msg("Packages: %s", paste(revdeps, collapse = ", "))
msg("Using %d parallel workers.", num_workers)

# -- install local lava into revdep library ----------------------------------

timestamp_msg("Installing local '%s' and its dependencies into revdep library...", pkg_name)

pak::local_install(
  root = pkg_dir,
  lib = lib_dir,
  upgrade = FALSE,
  ask = FALSE
)

# -- download revdep tarballs ------------------------------------------------

timestamp_msg("Downloading %d reverse dependency tarballs...", length(revdeps))

tarball_paths <- character(length(revdeps))
names(tarball_paths) <- revdeps

for (pkg in revdeps) {
  pkg_row <- available[pkg, ]
  tarball_name <- sprintf("%s_%s.tar.gz", pkg, pkg_row["Version"])
  dest <- file.path(src_dir, tarball_name)
  if (!file.exists(dest)) {
    url <- sprintf("%s/%s", pkg_row["Repository"], tarball_name)
    tryCatch(
      utils::download.file(url, dest, quiet = TRUE, mode = "wb"),
      error = function(e) msg("  Failed to download %s: %s", pkg, conditionMessage(e))
    )
  }
  if (file.exists(dest)) {
    tarball_paths[pkg] <- dest
  }
}

# remove packages that failed to download
failed_dl <- revdeps[!nzchar(tarball_paths)]
if (length(failed_dl) > 0L) {
  msg("WARNING: Failed to download: %s", paste(failed_dl, collapse = ", "))
}
revdeps <- revdeps[nzchar(tarball_paths)]
tarball_paths <- tarball_paths[revdeps]

# -- install revdep dependencies into shared library -------------------------

timestamp_msg("Installing dependencies for all reverse dependencies...")

# collect hard dependencies (Depends, Imports, LinkingTo) recursively
# we avoid Suggests because many revdeps suggest non-CRAN packages
all_revdep_deps <- unique(unlist(
  tools::package_dependencies(
    revdeps, db = available,
    which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  )
))

# filter to CRAN-available, non-base packages
base_pkgs <- rownames(installed.packages(priority = "base"))
to_install <- setdiff(intersect(all_revdep_deps, rownames(available)), base_pkgs)

if (length(to_install) > 0L) {
  timestamp_msg("Installing %d dependency packages (pak will skip up-to-date)...", length(to_install))
  pak::pkg_install(
    to_install,
    lib = lib_dir,
    upgrade = FALSE,
    ask = FALSE
  )
}

# -- run R CMD check on each revdep ------------------------------------------

timestamp_msg("Running R CMD check on %d packages (%d workers)...", length(revdeps), num_workers)

check_one <- function(pkg) {
  pkg_check_dir <- file.path(checks_dir, pkg)
  dir.create(pkg_check_dir, recursive = TRUE, showWarnings = FALSE)

  tarball <- tarball_paths[pkg]
  log_file <- file.path(pkg_check_dir, "check.log")

  # build check arguments
  check_args <- c(
    "--no-manual",
    "--no-multiarch",
    paste0("--library=", lib_dir)
  )
  if (skip_vignettes) {
    check_args <- c(check_args, "--no-vignettes", "--no-build-vignettes")
  }

  # set up environment so the check uses our library
  lib_path <- paste(c(lib_dir, .libPaths()), collapse = .Platform$path.sep)
  env_vars <- c(
    paste0("R_LIBS=", lib_path),
    paste0("R_LIBS_USER=", lib_dir),
    paste0("R_LIBS_SITE=", lib_path),
    "_R_CHECK_FORCE_SUGGESTS_=0"
  )

  cmd <- sprintf(
    "%s CMD check %s %s",
    file.path(R.home("bin"), "R"),
    paste(check_args, collapse = " "),
    shQuote(tarball)
  )

  t0 <- Sys.time()
  result <- tryCatch({
    system2(
      file.path(R.home("bin"), "R"),
      args = c("CMD", "check", check_args, tarball),
      stdout = log_file,
      stderr = log_file,
      env = env_vars,
      timeout = 1800
    )
  }, error = function(e) {
    writeLines(paste("check errored:", conditionMessage(e)), log_file)
    -1L
  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  # find the .Rcheck directory (created in current working dir)
  rcheck_dir <- file.path(getwd(), paste0(pkg, ".Rcheck"))
  check_log <- file.path(rcheck_dir, "00check.log")

  # copy the check log into our checks dir
  if (file.exists(check_log)) {
    file.copy(check_log, file.path(pkg_check_dir, "00check.log"), overwrite = TRUE)
  }

  # copy test output if it exists
  test_dir <- file.path(rcheck_dir, "tests")
  if (dir.exists(test_dir)) {
    file.copy(test_dir, pkg_check_dir, recursive = TRUE, overwrite = TRUE)
  }

  # parse results
  parsed <- if (file.exists(file.path(pkg_check_dir, "00check.log"))) {
    parse_check_log(file.path(pkg_check_dir, "00check.log"))
  } else {
    list(status = "ERROR", warnings = 0L, notes = 0L, errors = 1L,
         details = "no 00check.log produced")
  }

  # clean up .Rcheck directory
  unlink(rcheck_dir, recursive = TRUE)

  list(
    package = pkg,
    status = parsed$status,
    errors = parsed$errors,
    warnings = parsed$warnings,
    notes = parsed$notes,
    details = parsed$details,
    elapsed_min = round(elapsed, 1)
  )
}

# run checks in parallel
# note: mclapply uses forking, works on Linux/macOS; on Windows falls back to serial
results <- parallel::mclapply(revdeps, check_one, mc.cores = num_workers)

# -- build summary -----------------------------------------------------------

timestamp_msg("Building summary...")

summary_df <- data.frame(
  Package = vapply(results, `[[`, character(1), "package"),
  Status = vapply(results, `[[`, character(1), "status"),
  Errors = vapply(results, `[[`, integer(1), "errors"),
  Warnings = vapply(results, `[[`, integer(1), "warnings"),
  Notes = vapply(results, `[[`, integer(1), "notes"),
  Time_min = vapply(results, `[[`, numeric(1), "elapsed_min"),
  Details = vapply(results, `[[`, character(1), "details"),
  stringsAsFactors = FALSE
)

summary_df <- summary_df[order(summary_df$Status, summary_df$Package), ]

# write summary to file
summary_file <- file.path(revdep_root, "summary.md")
lines <- c(
  sprintf("# Reverse dependency check: %s %s", pkg_name,
          as.character(utils::packageVersion(pkg_name, lib.loc = lib_dir))),
  "",
  sprintf("Date: %s", Sys.Date()),
  sprintf("Workers: %d", num_workers),
  sprintf("Packages checked: %d", nrow(summary_df)),
  "",
  "## Results",
  "",
  sprintf("- OK: %d", sum(summary_df$Status == "OK")),
  sprintf("- NOTE: %d", sum(summary_df$Status == "NOTE")),
  sprintf("- WARNING: %d", sum(summary_df$Status == "WARNING")),
  sprintf("- ERROR: %d", sum(summary_df$Status == "ERROR")),
  "",
  "## Details",
  "",
  "| Package | Status | Errors | Warnings | Notes | Time (min) | Details |",
  "|---------|--------|--------|----------|-------|------------|---------|",
  sprintf("| %s | %s | %d | %d | %d | %.1f | %s |",
          summary_df$Package, summary_df$Status, summary_df$Errors,
          summary_df$Warnings, summary_df$Notes, summary_df$Time_min,
          summary_df$Details),
  "",
  "## Check logs",
  "",
  sprintf("Full check logs are in `revdep/checks/<package>/`.")
)

writeLines(lines, summary_file)

# print to stdout
msg("")
msg("=" |> rep(60) |> paste(collapse = ""))
msg("  REVERSE DEPENDENCY CHECK SUMMARY")
msg("=" |> rep(60) |> paste(collapse = ""))
msg("")
msg("  %-30s %s", "Package", "Status")
msg("  %-30s %s", "-------", "------")
for (i in seq_len(nrow(summary_df))) {
  status_marker <- switch(summary_df$Status[i],
    "OK" = "OK",
    "NOTE" = "NOTE",
    "WARNING" = "WARNING",
    "ERROR" = "** ERROR **"
  )
  msg("  %-30s %s (%.1f min)", summary_df$Package[i], status_marker, summary_df$Time_min[i])
}
msg("")
msg("  OK: %d | NOTE: %d | WARNING: %d | ERROR: %d",
    sum(summary_df$Status == "OK"), sum(summary_df$Status == "NOTE"),
    sum(summary_df$Status == "WARNING"), sum(summary_df$Status == "ERROR"))
msg("")
msg("  Full report: %s", summary_file)
msg("  Check logs:  %s", checks_dir)
msg("")

if (length(failed_dl) > 0L) {
  msg("  Skipped (download failed): %s", paste(failed_dl, collapse = ", "))
  msg("")
}

# exit with non-zero if any errors
n_errors <- sum(summary_df$Status == "ERROR")
if (n_errors > 0L) {
  msg("  %d package(s) had ERRORs.", n_errors)
  quit(status = 1L, save = "no", runLast = FALSE)
} else {
  msg("  All reverse dependencies passed.")
  quit(status = 0L, save = "no", runLast = FALSE)
}
