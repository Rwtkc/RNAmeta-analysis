gene_type_progress_start <- function(session) {
  progress <- shiny::Progress$new(session = session, min = 0, max = 100)
  progress$set(
    value = 0,
    message = "Running Gene Type Analysis",
    detail = "0% | Initializing analysis workspace"
  )
  progress
}

gene_type_progress_update <- function(progress, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  progress$set(
    value = value,
    message = "Running Gene Type Analysis",
    detail = sprintf("%s%% | %s", value, detail)
  )

  invisible(NULL)
}

gene_type_progress_close <- function(progress) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  Sys.sleep(2)
  progress$close()
  invisible(NULL)
}
