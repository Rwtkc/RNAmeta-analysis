`%||%` <- function(value, fallback) {
  if (is.null(value)) {
    return(fallback)
  }

  value
}

MOTIF_JOB_ACTIVE_LIMIT <- 3L
MOTIF_JOB_ACTIVE_STATUSES <- c("queued", "running")

motif_job_cache_root <- function() {
  cache_root <- app_path("cache", "motif_runtime_jobs")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  cache_root
}

motif_job_submission_lock_path <- function() {
  file.path(motif_job_cache_root(), ".submit-lock")
}

motif_job_with_submission_lock <- function(expr, timeout_seconds = 5) {
  lock_path <- motif_job_submission_lock_path()
  deadline <- as.numeric(Sys.time()) + timeout_seconds
  acquired <- FALSE

  while (as.numeric(Sys.time()) < deadline) {
    if (dir.create(lock_path, recursive = FALSE, showWarnings = FALSE)) {
      acquired <- TRUE
      break
    }

    Sys.sleep(0.1)
  }

  if (!acquired) {
    stop("The motif server is busy handling another submission. Please retry in a few seconds.")
  }

  on.exit(unlink(lock_path, recursive = TRUE, force = TRUE), add = TRUE)
  force(expr)
}

motif_job_normalize_id <- function(job_id) {
  job_value <- trimws(as.character(job_id %||% ""))
  if (!nzchar(job_value)) {
    return("")
  }

  gsub("[^A-Za-z0-9._-]+", "-", job_value)
}

motif_job_dir <- function(job_id) {
  normalized_id <- motif_job_normalize_id(job_id)
  if (!nzchar(normalized_id)) {
    stop("A valid motif job id is required.")
  }

  file.path(motif_job_cache_root(), normalized_id)
}

motif_job_state_path <- function(job_id) {
  file.path(motif_job_dir(job_id), "job.json")
}

motif_job_spec_path <- function(job_id) {
  file.path(motif_job_dir(job_id), "spec.rds")
}

motif_job_result_path <- function(job_id) {
  file.path(motif_job_dir(job_id), "result.rds")
}

motif_job_worker_stdout_path <- function(job_id) {
  file.path(motif_job_dir(job_id), "worker.stdout.log")
}

motif_job_worker_stderr_path <- function(job_id) {
  file.path(motif_job_dir(job_id), "worker.stderr.log")
}

motif_job_timestamp <- function(time = Sys.time()) {
  format(time, "%Y-%m-%d %H:%M:%S %z")
}

motif_job_random_token <- function(length = 6L) {
  paste(sample(c(letters, LETTERS, 0:9), size = length, replace = TRUE), collapse = "")
}

motif_job_new_id <- function() {
  sprintf(
    "motif-%s-%s",
    format(Sys.time(), "%Y%m%d-%H%M%S"),
    motif_job_random_token()
  )
}

motif_job_atomic_write_json <- function(path, data) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- sprintf("%s.tmp-%s", path, motif_job_random_token(8L))
  jsonlite::write_json(data, path = temp_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (!file.rename(temp_path, path)) {
    file.copy(temp_path, path, overwrite = TRUE)
    unlink(temp_path, force = TRUE)
  }
}

motif_job_atomic_save_rds <- function(path, data) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- sprintf("%s.tmp-%s", path, motif_job_random_token(8L))
  saveRDS(data, file = temp_path)
  if (!file.rename(temp_path, path)) {
    file.copy(temp_path, path, overwrite = TRUE)
    unlink(temp_path, force = TRUE)
  }
}

motif_job_load_state <- function(job_id) {
  state_path <- motif_job_state_path(job_id)
  if (!file.exists(state_path)) {
    return(NULL)
  }

  jsonlite::fromJSON(state_path, simplifyVector = FALSE)
}

motif_job_list_states <- function() {
  job_root <- motif_job_cache_root()
  job_dirs <- list.dirs(job_root, recursive = FALSE, full.names = TRUE)
  if (!length(job_dirs)) {
    return(list())
  }

  states <- lapply(job_dirs, function(job_dir) {
    state_path <- file.path(job_dir, "job.json")
    if (!file.exists(state_path)) {
      return(NULL)
    }

    tryCatch(
      jsonlite::fromJSON(state_path, simplifyVector = FALSE),
      error = function(...) NULL
    )
  })

  states <- Filter(Negate(is.null), states)
  if (exists("motif_job_reap_timed_out_states", mode = "function")) {
    states <- motif_job_reap_timed_out_states(states)
  }

  states
}

motif_job_active_states <- function() {
  Filter(function(state) {
    status_value <- tolower(as.character(state$status %||% ""))
    status_value %in% MOTIF_JOB_ACTIVE_STATUSES
  }, motif_job_list_states())
}

motif_job_assert_capacity <- function(limit = MOTIF_JOB_ACTIVE_LIMIT) {
  active_states <- motif_job_active_states()
  active_count <- length(active_states)

  if (active_count >= limit) {
    stop(sprintf(
      "The motif server is at capacity. %d motif jobs are already queued or running (global limit: %d). Please try again later.",
      active_count,
      limit
    ))
  }

  invisible(active_states)
}

motif_job_load_result <- function(job_id) {
  result_path <- motif_job_result_path(job_id)
  if (!file.exists(result_path)) {
    return(NULL)
  }

  readRDS(result_path)
}

motif_job_save_state <- function(job_id, state) {
  motif_job_atomic_write_json(motif_job_state_path(job_id), state)
  state
}

motif_job_update_state <- function(job_id, patch = list()) {
  current_state <- motif_job_load_state(job_id)
  if (is.null(current_state)) {
    stop(sprintf("Motif job '%s' does not exist.", job_id))
  }

  next_state <- utils::modifyList(current_state, patch)
  motif_job_save_state(job_id, next_state)
}

motif_job_submit_state <- function(upload_context, controls, selected_file, job_id) {
  engine_key <- tolower(as.character(controls$engine)[[1]])
  list(
    jobId = job_id,
    status = "queued",
    engine = engine_key,
    engineLabel = motif_engine_label(engine_key),
    sampleName = basename(selected_file$file_name),
    sourceFile = selected_file$file_name,
    species = upload_context$species,
    subsetLabel = motif_subset_label(controls$sites_subset),
    motifWidth = motif_effective_width(engine_key, controls$motif_width),
    maxMotifs = max(1L, motif_safe_integer(controls$max_motifs, 6L)),
    timeoutSeconds = motif_engine_timeout_seconds(engine_key),
    timeoutMinutes = motif_engine_timeout_seconds(engine_key) / 60L,
    submittedAt = motif_job_timestamp(),
    startedAt = NULL,
    finishedAt = NULL,
    message = sprintf("Queued %s motif analysis.", motif_engine_label(engine_key)),
    workspaceDir = normalizePath(motif_job_dir(job_id), winslash = "/", mustWork = FALSE),
    resultDir = NULL,
    errorMessage = NULL
  )
}

motif_job_rscript_path <- function() {
  configured_path <- "/public/liuqi/miniconda3/envs/R4.3.0/bin/Rscript"
  if (file.exists(configured_path)) {
    return(configured_path)
  }

  fallback_path <- file.path(R.home("bin"), "Rscript")
  if (file.exists(fallback_path)) {
    return(fallback_path)
  }

  stop("Unable to locate Rscript for motif background jobs.")
}

motif_job_launch <- function(job_id) {
  rscript_path <- motif_job_rscript_path()
  worker_script <- app_path("services", "motif_job_worker.R")
  worker_stdout <- motif_job_worker_stdout_path(job_id)
  worker_stderr <- motif_job_worker_stderr_path(job_id)
  worker_args <- c(worker_script, job_id, app_root)

  dir.create(motif_job_dir(job_id), recursive = TRUE, showWarnings = FALSE)
  unlink(worker_stdout, force = TRUE)
  unlink(worker_stderr, force = TRUE)

  if (.Platform$OS.type == "windows") {
    argument_list <- paste(
      vapply(
        c(worker_script, job_id, app_root),
        function(value) shQuote(value),
        character(1)
      ),
      collapse = ", "
    )
    command <- sprintf(
      "Start-Process -FilePath %s -ArgumentList @(%s) -WorkingDirectory %s -RedirectStandardOutput %s -RedirectStandardError %s -WindowStyle Hidden",
      shQuote(rscript_path),
      argument_list,
      shQuote(app_root),
      shQuote(worker_stdout),
      shQuote(worker_stderr)
    )
    invisible(system2("powershell", c("-NoProfile", "-Command", command), wait = FALSE))
    return(invisible(TRUE))
  }

  process_id <- suppressWarnings(
    system2(
      command = rscript_path,
      args = worker_args,
      stdout = worker_stdout,
      stderr = worker_stderr,
      wait = FALSE
    )
  )

  is.numeric(process_id) || is.integer(process_id)
}

motif_submit_background_job <- function(upload_context, controls) {
  motif_job_with_submission_lock({
    controls <- motif_validate_analysis_controls(controls)
    saved_files <- upload_run_context_saved_files(upload_context)
    if (!length(saved_files)) {
      stop("Save a BED file in Upload / Run before running motif analysis.")
    }

    motif_job_assert_capacity()

    selected_file <- motif_selected_saved_file(upload_context, controls$sample_name)
    if (is.null(selected_file) || !nzchar(selected_file$file_path) || !file.exists(selected_file$file_path)) {
      stop("No saved BED file is available for motif analysis.")
    }

    job_id <- motif_job_new_id()
    job_dir <- motif_job_dir(job_id)
    dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
    worker_upload_context <- motif_job_worker_upload_context(upload_context, selected_file, job_dir, controls)

    state <- motif_job_submit_state(upload_context, controls, selected_file, job_id)
    spec <- list(
      job_id = job_id,
      upload_context = worker_upload_context,
      controls = controls
    )

    motif_job_atomic_save_rds(motif_job_spec_path(job_id), spec)
    motif_job_save_state(job_id, state)
    motif_job_launch(job_id)

    list(
      job_id = job_id,
      state = state
    )
  })
}

motif_job_query <- function(job_id) {
  normalized_id <- motif_job_normalize_id(job_id)
  if (!nzchar(normalized_id)) {
    return(NULL)
  }

  list(
    state = motif_job_load_state(normalized_id),
    result = motif_job_load_result(normalized_id)
  )
}

motif_job_load_spec <- function(job_id) {
  spec_path <- motif_job_spec_path(job_id)
  if (!file.exists(spec_path)) {
    stop(sprintf("Motif job spec for '%s' is missing.", job_id))
  }

  readRDS(spec_path)
}

motif_job_result_state_patch <- function(result) {
  list(
    status = "completed",
    startedAt = if (is.null(result$started_at)) NULL else result$started_at,
    finishedAt = motif_job_timestamp(),
    message = sprintf(
      "%s completed with %d parsed motif%s.",
      result$engine_label,
      length(result$motifs %||% list()),
      if (length(result$motifs %||% list()) == 1L) "" else "s"
    ),
    resultDir = result$result_output_dir %||% result$job_package_dir,
    errorMessage = NULL
  )
}

motif_job_failure_state_patch <- function(message, started_at = NULL) {
  list(
    status = "failed",
    startedAt = started_at,
    finishedAt = motif_job_timestamp(),
    message = "Motif analysis failed.",
    errorMessage = as.character(message)[[1]]
  )
}

motif_background_worker_run <- function(job_id) {
  spec <- motif_job_load_spec(job_id)
  started_at <- motif_job_timestamp()
  state <- motif_job_update_state(job_id, list(
    status = "running",
    startedAt = started_at,
    message = sprintf("Running %s on the server workspace.", motif_engine_label(spec$controls$engine)),
    errorMessage = NULL
  ))

  result <- tryCatch(
    run_motif_analysis(
      upload_context = spec$upload_context,
      controls = spec$controls,
      progress = NULL,
      job_id = job_id
    ),
    error = function(err) {
      motif_job_update_state(job_id, motif_job_failure_state_patch(conditionMessage(err), started_at = started_at))
      NULL
    }
  )

  if (is.null(result)) {
    return(invisible(NULL))
  }

  motif_job_atomic_save_rds(motif_job_result_path(job_id), result)
  motif_job_update_state(job_id, motif_job_result_state_patch(result))
  invisible(TRUE)
}
