motif_progress_start <- function(session) {
  progress <- shiny::Progress$new(session = session, min = 0, max = 100)
  progress$set(
    value = 0,
    message = "Running Motif Analysis",
    detail = "0% | Initializing motif workspace"
  )
  progress
}

motif_progress_update <- function(progress, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  safe_value <- as.integer(round(max(0, min(100, value))))
  progress$set(
    value = safe_value,
    message = "Running Motif Analysis",
    detail = sprintf("%s%% | %s", safe_value, detail)
  )

  invisible(NULL)
}

motif_progress_close <- function(progress) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  Sys.sleep(1)
  progress$close()
  invisible(NULL)
}
