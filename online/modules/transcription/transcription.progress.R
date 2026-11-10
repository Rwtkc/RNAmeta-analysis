transcription_progress_start <- function(session) {
  progress <- shiny::Progress$new(session = session, min = 0, max = 100)
  attr(progress, "rnameta_last_value") <- 0L
  progress$set(
    value = 0,
    message = "Running Transcription Analysis",
    detail = "0% | Initializing analysis workspace"
  )
  progress
}

transcription_progress_update <- function(progress, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  last_value <- attr(progress, "rnameta_last_value", exact = TRUE)
  if (!is.numeric(last_value) || length(last_value) != 1L || is.na(last_value)) {
    last_value <- 0
  }

  normalized_value <- as.integer(round(max(0, min(100, value))))
  display_value <- max(last_value, normalized_value)
  attr(progress, "rnameta_last_value") <- display_value

  progress$set(
    value = display_value,
    message = "Running Transcription Analysis",
    detail = sprintf("%s%% | %s", display_value, detail)
  )

  invisible(NULL)
}

transcription_progress_close <- function(progress) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  Sys.sleep(2)
  progress$close()
  invisible(NULL)
}
