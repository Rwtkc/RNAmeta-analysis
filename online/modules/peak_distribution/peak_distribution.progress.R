peak_distribution_progress_start <- function(session) {
  progress <- shiny::Progress$new(session = session, min = 0, max = 100)
  progress$set(
    value = 0,
    message = "Running Peak Distribution Analysis",
    detail = "0% | Initializing analysis workspace"
  )
  progress
}

peak_distribution_progress_update <- function(progress, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  progress$set(
    value = value,
    message = "Running Peak Distribution Analysis",
    detail = sprintf("%s%% | %s", value, detail)
  )

  invisible(NULL)
}

peak_distribution_progress_close <- function(progress) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  Sys.sleep(2)
  progress$close()
  invisible(NULL)
}
