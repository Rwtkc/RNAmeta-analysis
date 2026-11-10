splicesite_progress_start <- function(session) {
  progress <- shiny::Progress$new(session = session, min = 0, max = 100)
  progress$set(
    value = 0,
    message = "Running Splicesite Analysis",
    detail = "0% | Initializing analysis workspace"
  )
  progress
}

splicesite_progress_update <- function(progress, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  progress$set(
    value = value,
    message = "Running Splicesite Analysis",
    detail = sprintf("%s%% | %s", value, detail)
  )

  invisible(NULL)
}

splicesite_progress_close <- function(progress) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  Sys.sleep(2)
  progress$close()
  invisible(NULL)
}
