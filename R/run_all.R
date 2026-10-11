#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage:",
    "Rscript R/run_all.R <GSE206528_counts.tsv.gz> <GSE206528_schema.tsv> <GSE144803_non-normalized.txt.gz> <output_dir>"
  ))
}
counts_file <- normalizePath(args[[1]], mustWork=TRUE)
schema_file <- normalizePath(args[[2]], mustWork=TRUE)
tnf_file <- normalizePath(args[[3]], mustWork=TRUE)
outdir <- normalizePath(args[[4]], winslash="/", mustWork=FALSE)

# Resolve script directory independently of working directory.
cmd <- commandArgs(trailingOnly=FALSE)
file_arg <- sub("^--file=","",cmd[grep("^--file=",cmd)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1]],mustWork=TRUE)) else getwd()

source(file.path(script_dir,"00_helpers.R"))
source(file.path(script_dir,"01_GSE206528_primary_DE_and_robustness.R"))
source(file.path(script_dir,"02_GSE206528_upstream_triad.R"))
source(file.path(script_dir,"03_GSE144803_TNF_CEMIP.R"))
source(file.path(script_dir,"04_audit_against_historical.R"))
source(file.path(script_dir,"05_session_info.R"))

if (!requireNamespace("edgeR",quietly=TRUE)) {
  stop("Required Bioconductor package 'edgeR' is not installed. No fallback/substitution is permitted.")
}

safe_dir(outdir)
results_dir <- file.path(outdir,"RESULTS")
audit_dir <- file.path(outdir,"AUDIT")
runtime_dir <- file.path(outdir,"RUNTIME")
safe_dir(results_dir); safe_dir(audit_dir); safe_dir(runtime_dir)

run_gse206528_primary(counts_file, schema_file, results_dir)
run_gse206528_upstream_triad(counts_file, schema_file, results_dir)
run_gse144803_tnf_cemip(tnf_file, results_dir)
aud <- run_final_audit(results_dir, audit_dir)
write_session_info(runtime_dir)

status <- if (isTRUE(aud$material)) "MATERIAL_DISCREPANCY_STOP" else "READY_FOR_MANUSCRIPT_DRAFTING"
writeLines(status, file.path(outdir,"FINAL_STATUS.txt"))

if (isTRUE(aud$material)) {
  stop("Final audit detected a MATERIAL DISCREPANCY. Manuscript integration must stop; see AUDIT/FINAL_AUDIT_REPORT.md")
}

cat("\n# READY FOR MANUSCRIPT DRAFTING\n")
