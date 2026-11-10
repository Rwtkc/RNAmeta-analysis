structure_job_cache_root <- function() {
  cache_root <- app_path("cache", "structure_runtime_jobs")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  cache_root
}

structure_job_submission_lock_path <- function() {
  file.path(structure_job_cache_root(), ".submit-lock")
}

structure_job_with_submission_lock <- function(expr, timeout_seconds = 5) {
  lock_path <- structure_job_submission_lock_path()
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
    stop("The Structure server is busy handling another submission. Please retry in a few seconds.")
  }

  on.exit(unlink(lock_path, recursive = TRUE, force = TRUE), add = TRUE)
  force(expr)
}

structure_job_normalize_id <- function(job_id) {
  cleaned <- trimws(as.character(job_id %||% ""))
  cleaned <- gsub("[^A-Za-z0-9_-]", "", cleaned)

  if (!startsWith(cleaned, "structure_")) {
    return("")
  }

  cleaned
}

structure_job_dir <- function(job_id) {
  normalized_id <- structure_job_normalize_id(job_id)
  if (!nzchar(normalized_id)) {
    stop("A valid Structure job id is required.")
  }

  file.path(structure_job_cache_root(), normalized_id)
}

structure_job_timestamp <- function(time = Sys.time()) {
  format(time, "%Y-%m-%d %H:%M:%S %z")
}

structure_job_expires_at <- function() {
  structure_job_timestamp(Sys.time() + as.difftime(7, units = "days"))
}

structure_job_random_token <- function(length = 6L) {
  paste(sample(c(letters, LETTERS, 0:9), size = length, replace = TRUE), collapse = "")
}

structure_job_new_id <- function() {
  sprintf(
    "structure_%s_%s",
    format(Sys.time(), "%Y%m%d%H%M%S"),
    structure_job_random_token()
  )
}

structure_job_state_path <- function(job_id) {
  file.path(structure_job_dir(job_id), "job.json")
}

structure_job_spec_path <- function(job_id) {
  file.path(structure_job_dir(job_id), "spec.rds")
}

structure_job_result_path <- function(job_id) {
  file.path(structure_job_dir(job_id), "result.rds")
}

structure_job_worker_stdout_path <- function(job_id) {
  file.path(structure_job_dir(job_id), "worker.stdout.log")
}

structure_job_worker_stderr_path <- function(job_id) {
  file.path(structure_job_dir(job_id), "worker.stderr.log")
}

structure_job_atomic_write_json <- function(path, data) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- tempfile(pattern = "structure_job_", tmpdir = dirname(path), fileext = ".json")
  jsonlite::write_json(data, temp_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (!file.rename(temp_path, path)) {
    file.copy(temp_path, path, overwrite = TRUE)
    unlink(temp_path, force = TRUE)
  }

  invisible(path)
}

structure_job_atomic_save_rds <- function(path, data) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- tempfile(pattern = "structure_result_", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(data, temp_path)
  if (!file.rename(temp_path, path)) {
    file.copy(temp_path, path, overwrite = TRUE)
    unlink(temp_path, force = TRUE)
  }

  invisible(path)
}

structure_job_save <- function(job_id, state, result = NULL) {
  structure_job_atomic_write_json(structure_job_state_path(job_id), state)
  if (!is.null(result)) {
    structure_job_atomic_save_rds(structure_job_result_path(job_id), result)
  }

  invisible(list(state = state, result = result))
}

structure_job_load_state <- function(job_id) {
  state_path <- structure_job_state_path(job_id)
  if (!file.exists(state_path)) {
    return(NULL)
  }

  tryCatch(
    jsonlite::fromJSON(state_path, simplifyVector = FALSE),
    error = function(...) NULL
  )
}

structure_job_load_result <- function(job_id) {
  result_path <- structure_job_result_path(job_id)
  if (!file.exists(result_path)) {
    return(NULL)
  }

  tryCatch(readRDS(result_path), error = function(...) NULL)
}

structure_job_save_state <- function(job_id, state) {
  structure_job_atomic_write_json(structure_job_state_path(job_id), state)
  state
}

structure_job_update_state <- function(job_id, patch = list()) {
  current_state <- structure_job_load_state(job_id)
  if (is.null(current_state)) {
    stop(sprintf("Structure job '%s' does not exist.", job_id))
  }

  next_state <- utils::modifyList(current_state, patch)
  structure_job_save_state(job_id, next_state)
}

structure_job_list_states <- function() {
  job_root <- structure_job_cache_root()
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

  Filter(Negate(is.null), states)
}

structure_job_active_states <- function() {
  Filter(function(state) {
    tolower(as.character(state$status %||% "")) %in% STRUCTURE_JOB_ACTIVE_STATUSES
  }, structure_job_list_states())
}

structure_job_capacity_snapshot <- function(limit = STRUCTURE_JOB_ACTIVE_LIMIT) {
  states <- structure_job_list_states()
  statuses <- vapply(
    states,
    function(state) tolower(as.character(state$status %||% "")),
    character(1)
  )
  queued_count <- sum(statuses == "queued")
  running_count <- sum(statuses == "running")
  active_count <- queued_count + running_count
  limit_value <- as.integer(limit)

  list(
    queued = as.integer(queued_count),
    running = as.integer(running_count),
    active = as.integer(active_count),
    limit = limit_value,
    available = as.integer(max(0L, limit_value - active_count))
  )
}

structure_job_assert_capacity <- function(limit = STRUCTURE_JOB_ACTIVE_LIMIT) {
  active_states <- structure_job_active_states()
  active_count <- length(active_states)
  if (active_count >= limit) {
    stop(sprintf(
      "The Structure server is at capacity. %d Structure jobs are already queued or running (global limit: %d). Please try again later.",
      active_count,
      limit
    ))
  }

  invisible(active_states)
}

structure_job_query <- function(job_id) {
  if (exists("rnameta_cache_cleanup_due", mode = "function")) {
    rnameta_cache_cleanup_due()
  }

  normalized_id <- structure_job_normalize_id(job_id)
  if (!nzchar(normalized_id)) {
    return(NULL)
  }

  list(
    state = structure_job_load_state(normalized_id),
    result = structure_job_load_result(normalized_id)
  )
}

structure_result_payload <- function(result) {
  transcript <- result$extracted$transcript %||% list()

  list(
    sequence = result$extracted$sequence,
    structure = result$folded$structure,
    fileName = result$selection$file$file_name %||% "",
    coordinate = result$selection$coordinate$raw %||% "",
    locus = result$extracted$locus,
    bedLocus = result$extracted$bed_locus %||% "",
    strand = result$extracted$strand,
    sequenceLength = nchar(result$extracted$sequence),
    method = result$folded$method,
    windowSize = result$extracted$window_size,
    intervalLength = result$extracted$interval_length,
    usedBedInterval = isTRUE(result$extracted$used_bed_interval),
    matchedBedRow = isTRUE(result$extracted$matched_bed_row),
    transcriptId = transcript$transcript_id %||% "",
    transcriptName = transcript$transcript_name %||% "",
    geneId = transcript$gene_id %||% "",
    biotype = transcript$biotype %||% "",
    transcriptLength = transcript$transcript_length %||% nchar(result$extracted$sequence),
    highlightStart = transcript$highlight_start %||% NULL,
    highlightEnd = transcript$highlight_end %||% NULL,
    highlightRanges = result$extracted$highlight_ranges %||% list(),
    overlapBases = transcript$overlap_bases %||% 0L,
    candidateCount = transcript$candidate_count %||% 0L
  )
}
