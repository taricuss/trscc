#!/usr/bin/env Rscript

# Convenience wrapper to knit the manuscript from anywhere in the repo.
# Usage:
#   Rscript manuscript/render_manuscript.R [html|word|all]

suppressPackageStartupMessages({
  library(rmarkdown)
  library(knitr)
})

args <- commandArgs(trailingOnly = TRUE)
output_format <- if (length(args) >= 1L) args[1L] else "all"

repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = FALSE)
if (is.na(repo_root) || !dir.exists(repo_root)) {
  repo_root <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
}

rmd_path <- file.path(repo_root, "manuscript", "TRCSS_manuscript.Rmd")
if (!file.exists(rmd_path)) {
  stop("Manuscript not found at ", rmd_path, call. = FALSE)
}

out_dir <- file.path(repo_root, "manuscript")

cat(sprintf("Repo root : %s\n", repo_root))
cat(sprintf("Rmd path  : %s\n", rmd_path))
cat(sprintf("Output to : %s\n", out_dir))
cat(sprintf("Format    : %s\n", output_format))

render_one <- function(fmt) {
  cat(sprintf("\n=== Rendering %s ===\n", fmt))
  format <- switch(fmt,
    html = "html_document",
    word = "word_document",
    pdf  = "pdf_document",
    stop("Unknown format: ", fmt, call. = FALSE))
  rmarkdown::render(
    input = rmd_path,
    output_format = format,
    output_dir = out_dir,
    clean = TRUE,
    quiet = FALSE,
    envir = new.env(parent = globalenv())
  )
}

if (output_format == "all") {
  render_one("html")
  render_one("word")
} else {
  render_one(output_format)
}

cat("\nDone.\n")
