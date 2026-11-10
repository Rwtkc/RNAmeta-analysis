MOTIF_ENGINE_TIMEOUT_SECONDS <- 40L * 60L
MOTIF_ENGINE_TIMEOUT_MINUTES <- MOTIF_ENGINE_TIMEOUT_SECONDS / 60L

motif_engine_timeout_seconds <- function(engine = NULL) {
  as.integer(MOTIF_ENGINE_TIMEOUT_SECONDS)
}

motif_job_timeout_field <- function(state, name, fallback = "") {
  value <- state[[name]]
  if (is.null(value) || !length(value) || is.na(value[[1]])) {
    return(fallback)
  }

  as.character(value[[1]])
}

motif_job_timeout_timestamp <- function(time = Sys.time()) {
  format(time, "%Y-%m-%d %H:%M:%S %z")
}

motif_job_parse_timestamp <- function(value) {
  if (is.null(value) || !length(value) || is.na(value[[1]])) {
    return(as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC"))
  }

  timestamp <- trimws(as.character(value[[1]]))
  if (!nzchar(timestamp)) {
    return(as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC"))
  }

  parsed <- suppressWarnings(as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S %z", tz = "UTC"))
  if (!is.na(parsed)) {
    return(parsed)
  }

  suppressWarnings(as.POSIXct(timestamp, tz = "UTC"))
}

motif_job_timeout_reference_time <- function(state) {
  started_at <- motif_job_parse_timestamp(state$startedAt)
  if (!is.na(started_at)) {
    return(started_at)
  }

  motif_job_parse_timestamp(state$submittedAt)
}

motif_job_is_timed_out <- function(state, now = Sys.time(), timeout_seconds = MOTIF_ENGINE_TIMEOUT_SECONDS) {
  status <- tolower(motif_job_timeout_field(state, "status"))
  if (!status %in% c("queued", "running")) {
    return(FALSE)
  }

  reference_time <- motif_job_timeout_reference_time(state)
  if (is.na(reference_time)) {
    return(FALSE)
  }

  elapsed_seconds <- as.numeric(difftime(now, reference_time, units = "secs"))
  is.finite(elapsed_seconds) && elapsed_seconds >= as.numeric(timeout_seconds)
}

motif_job_timeout_state_patch <- function(state, now = Sys.time()) {
  engine_label <- motif_job_timeout_field(state, "engineLabel", "Motif job")

  list(
    status = "failed",
    finishedAt = motif_job_timeout_timestamp(now),
    message = "Motif analysis timed out.",
    errorMessage = sprintf(
      "%s exceeded the %d minute runtime limit and was marked failed automatically.",
      engine_label,
      MOTIF_ENGINE_TIMEOUT_MINUTES
    )
  )
}

motif_job_reap_timed_out_states <- function(states, now = Sys.time()) {
  lapply(states, function(state) {
    if (!motif_job_is_timed_out(state, now = now)) {
      return(state)
    }

    job_id <- motif_job_timeout_field(state, "jobId")
    if (!nzchar(job_id)) {
      return(state)
    }

    next_state <- utils::modifyList(state, motif_job_timeout_state_patch(state, now = now))
    tryCatch(
      motif_job_save_state(job_id, next_state),
      error = function(...) next_state
    )
  })
}
