#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

job_id <- if (length(args) >= 1L) args[[1]] else ""
app_root_path <- if (length(args) >= 2L) args[[2]] else ""

if (!nzchar(job_id) || !nzchar(app_root_path)) {
  stop("Usage: structure_job_worker.R <job_id> <app_root>")
}

setwd(app_root_path)
source(file.path(app_root_path, "global.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "structure", "structure.server.R"), local = TRUE, encoding = "UTF-8")

structure_job_worker_bootstrap(app_root_path, job_id)
