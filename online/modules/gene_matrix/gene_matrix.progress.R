gene_matrix_progress_start <- function(session) {
  shiny::Progress$new(
    session = session,
    min = 0,
    max = 100,
    style = "notification"
  )
}

gene_matrix_progress_update <- function(progress, value, detail = NULL) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  progress$set(
    value = value,
    message = "Running Gene Matrix Analysis",
    detail = detail
  )
}

gene_matrix_progress_close <- function(progress) {
  if (!is.null(progress)) {
    progress$close()
  }

  invisible(NULL)
}
