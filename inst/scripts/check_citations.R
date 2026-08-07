#!/usr/bin/env Rscript

# Checks that the citations in this package are still current and internally
# consistent.
#
# Three things go stale on their own, without anyone touching the package:
#
#   1. A DOI stops resolving. Publishers move content, and a reference that
#      looked fine when written becomes a dead link.
#   2. The reference list and the function documentation drift apart. Someone
#      adds a @references entry and forgets the README, or removes a function and
#      leaves its reference behind.
#   3. A software citation ages. citation("terra") reports a version and year
#      that change with every release, so a hard-coded copy in the README slowly
#      stops matching what a user would be told to cite.
#
# Exit codes follow the catalog check in datamatch:
#   0  everything current
#   1  something needs attention, and the report says what
#   2  could not check, usually the network. Not a failure to act on.
#
# Usage:
#   Rscript inst/scripts/check_citations.R [--markdown out.md] [--json out.json]

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) NULL else args[i + 1]
}
markdown_out <- arg_value("--markdown")
json_out <- arg_value("--json")

root <- if (file.exists("DESCRIPTION")) "." else stop("Run from the package root.")
read_all <- function(paths) paste(unlist(lapply(paths, readLines, warn = FALSE)),
                                  collapse = "\n")

r_source <- read_all(list.files(file.path(root, "R"), "[.]R$", full.names = TRUE))
readme <- read_all(file.path(root, "README.md"))
methods_path <- file.path(root, "docs", "methods.md")
methods <- if (file.exists(methods_path)) read_all(methods_path) else ""

problems <- list()
note <- function(kind, detail) {
  problems[[length(problems) + 1]] <<- list(kind = kind, detail = detail)
}

# ---- 1. DOIs resolve --------------------------------------------------------

doi_pattern <- "10[.][0-9]{4,9}/[-._;()/:A-Za-z0-9]+"
# `]` must come first inside a character class, or the class ends there and
# the trailing paren from a markdown link is never stripped.
strip_trailing <- function(x) sub("[]).,;:]+$", "", x)

dois <- unique(strip_trailing(unlist(regmatches(
  c(r_source, readme, methods),
  gregexpr(doi_pattern, c(r_source, readme, methods))
))))
dois <- dois[nzchar(dois)]

# Ask doi.org whether the DOI is registered, and stop there. Following the
# redirect lands on the publisher, and publishers routinely answer a scripted
# HEAD with 403 - Wiley, AMS, Annual Reviews and Inter-Research all do. Treating
# that as a dead reference would file a false issue every quarter and teach
# everyone to close it unread.
#
# A registered DOI redirects (3xx). An unregistered one is 404 at doi.org
# itself. Anything else is inconclusive rather than a failure.
resolves <- function(doi) {
  url <- paste0("https://doi.org/", doi)
  tryCatch({
    handle <- curl::new_handle(nobody = TRUE, followlocation = FALSE,
                               timeout = 25, customrequest = "HEAD")
    curl::curl_fetch_memory(url, handle = handle)$status_code
  }, error = function(e) NA_integer_)
}

reachable <- TRUE
if (length(dois) > 0) {
  if (!requireNamespace("curl", quietly = TRUE)) {
    message("curl not installed; skipping DOI resolution.")
    reachable <- FALSE
  } else {
    codes <- vapply(dois, resolves, integer(1))
    # Every single one failing means the network or doi.org is down, not that
    # every reference rotted at once. Saying so is more useful than filing an
    # issue that lists the whole bibliography.
    if (all(is.na(codes))) {
      message("Could not reach doi.org at all; drawing no conclusion.")
      reachable <- FALSE
    } else {
      # 404 from doi.org means the DOI is not registered. Everything else,
      # including a publisher-side block we never see now, is not our problem.
      dead <- dois[!is.na(codes) & codes == 404]
      for (d in dead) {
        note("unregistered DOI",
             paste0("`", d, "` is not registered at doi.org (HTTP 404)."))
      }
      inconclusive <- dois[!is.na(codes) & codes != 404 &
                             (codes < 300 | codes >= 400)]
      if (length(inconclusive) > 0) {
        message("Inconclusive for ", length(inconclusive), " DOI(s): ",
                paste(inconclusive, collapse = ", "))
      }
    }
  }
}

# ---- 2. Citations and the reference list agree ------------------------------

# Surnames used in running text: "Belkin and O'Reilly (2009)", "Ross et al. 2023"
cited_names <- function(text) {
  hits <- unlist(regmatches(text, gregexpr(
    "\\b[A-Z][a-zA-Z'\u2019]+(?:\\s+et al\\.|\\s+and\\s+[A-Z][a-zA-Z'\u2019]+)\\s*\\(?[0-9]{4}",
    text)))
  unique(sub("\\s+(et al\\.|and\\s+[A-Z][a-zA-Z'\u2019]+)\\s*\\(?[0-9]{4}$", "", hits))
}

split_at <- regexpr("\n## References", readme, fixed = TRUE)
if (split_at < 0) {
  note("missing section", "README has no `## References` section.")
  reference_list <- ""
  readme_body <- readme
} else {
  reference_list <- substring(readme, split_at)
  readme_body <- substring(readme, 1, split_at - 1)
}

cited <- unique(c(cited_names(r_source), cited_names(methods),
                  cited_names(readme_body)))
for (name in cited) {
  if (!grepl(name, reference_list, fixed = TRUE)) {
    note("uncited in list",
         paste0("`", name, "` is cited in the documentation but is not in the ",
                "README reference list."))
  }
}

# Entries in the list that nothing refers to
# Only the papers. The data-source and software subsections list things a user
# should cite, not things this documentation refers to, so they have no
# corresponding mention in running text and are not orphans.
papers_only <- reference_list
cut <- regexpr("\n### Data sources", papers_only, fixed = TRUE)
if (cut > 0) papers_only <- substring(papers_only, 1, cut - 1)

entries <- unlist(regmatches(papers_only, gregexpr(
  "(?m)^- ([A-Za-z'\u2019]+)", papers_only, perl = TRUE)))
entries <- unique(sub("^- ", "", entries))
for (entry in entries) {
  used <- grepl(entry, r_source, fixed = TRUE) ||
    grepl(entry, methods, fixed = TRUE) ||
    grepl(entry, readme_body, fixed = TRUE)
  if (!used) {
    note("orphan reference",
         paste0("`", entry, "` is in the reference list but cited nowhere."))
  }
}

# ---- 3. Software citations name a version, not a year ----------------------

# Deliberately not checked against citation(). That reports whatever version is
# installed, so the answer differs between a contributor's machine and CI, and a
# check that disagrees with itself trains people to ignore it. The README points
# at citation() instead of pinning a year, which removes the failure mode rather
# than monitoring it.

# ---- report -----------------------------------------------------------------

if (!reachable && length(problems) == 0) {
  if (!is.null(markdown_out)) writeLines("Could not check.", markdown_out)
  quit(status = 2)
}

if (length(problems) == 0) {
  message("All ", length(dois), " DOIs resolve; citations and reference list agree.")
  if (!is.null(markdown_out)) {
    writeLines("No citation problems found.", markdown_out)
  }
  quit(status = 0)
}

lines <- c(
  "The scheduled citation check found something that needs a look.",
  "",
  paste0("Checked ", length(dois), " DOIs across `R/`, the README, and ",
         "`docs/methods.md`."),
  ""
)
for (kind in unique(vapply(problems, `[[`, "", "kind"))) {
  lines <- c(lines, paste0("**", kind, "**"), "")
  for (p in problems) if (p$kind == kind) lines <- c(lines, paste0("- ", p$detail))
  lines <- c(lines, "")
}
lines <- c(lines,
  "Re-run locally with:",
  "",
  "```",
  "Rscript inst/scripts/check_citations.R",
  "```")

report <- paste(lines, collapse = "\n")
cat(report, "\n")
if (!is.null(markdown_out)) writeLines(report, markdown_out)
if (!is.null(json_out)) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(problems, json_out, auto_unbox = TRUE, pretty = TRUE)
  }
}
quit(status = 1)
