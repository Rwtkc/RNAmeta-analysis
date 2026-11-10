motif_job_input_filename <- function(file_name) {
  cleaned <- gsub("[^A-Za-z0-9._-]+", "_", basename(as.character(file_name %||% "input.bed")[[1]]))
  if (!nzchar(cleaned)) {
    cleaned <- "input.bed"
  }

  cleaned
}

motif_job_worker_upload_context <- function(upload_context, selected_file, job_dir, controls = list()) {
  input_dir <- file.path(job_dir, "inputs")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  files <- if (identical(controls$analysis_mode, "matched_full")) upload_run_context_saved_files(upload_context) else list(selected_file)
  copied_files <- lapply(seq_along(files), function(i) {
    file <- files[[i]]
    target_dir <- file.path(input_dir, as.character(i))
    dir.create(target_dir, showWarnings = FALSE)
    target_path <- file.path(target_dir, motif_job_input_filename(file$file_name))
    if (!file.copy(file$file_path, target_path, overwrite = FALSE)) stop("Failed to snapshot a saved BED file into the motif job workspace.")
    file$file_path <- normalizePath(target_path, winslash = "/", mustWork = TRUE)
    file
  })
  selected_index <- which(vapply(files, function(file) identical(file$file_path, selected_file$file_path), logical(1)))[[1]]
  copied_file <- copied_files[[selected_index]]

  worker_context <- upload_context
  worker_context$file_count <- length(copied_files)
  worker_context$file_name <- copied_file$file_name
  worker_context$file_path <- copied_file$file_path
  worker_context$file_names <- vapply(copied_files, function(file) file$file_name, character(1))
  worker_context$saved_files <- copied_files
  worker_context$preview <- NULL
  worker_context
}
