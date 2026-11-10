motif_placeholder_chart_payload <- function(job_state = NULL) {
  if (is.null(job_state)) {
    return(NULL)
  }

  list(
    status = job_state$status %||% "idle",
    jobId = job_state$jobId %||% "",
    jobStatus = job_state$status %||% "idle",
    jobMessage = job_state$message %||% "",
    engineLabel = job_state$engineLabel %||% "",
    startedAt = job_state$startedAt %||% NULL,
    finishedAt = job_state$finishedAt %||% NULL,
    motifs = list()
  )
}

motif_job_notice_card <- function(job_state = NULL, error_message = NULL) {
  if (!is.null(error_message) && nzchar(error_message)) {
    return(
      div(
        class = "site-profile-result-card site-profile-result-card--error",
        tags$p(class = "site-profile-result-card__copy", error_message)
      )
    )
  }

  if (is.null(job_state)) {
    return(NULL)
  }

  status_label <- toupper(job_state$status %||% "idle")
  job_copy <- c(
    sprintf("Job ID: %s", job_state$jobId %||% "pending"),
    if (!is.null(job_state$engineLabel) && nzchar(job_state$engineLabel)) {
      sprintf("Engine: %s", job_state$engineLabel)
    },
    if (!is.null(job_state$message) && nzchar(job_state$message)) {
      job_state$message
    },
    if (!is.null(job_state$errorMessage) && nzchar(job_state$errorMessage)) {
      sprintf("Error: %s", job_state$errorMessage)
    }
  )
  job_copy <- Filter(Negate(is.null), job_copy)

  div(
    class = "site-profile-result-card",
    tags$h4(class = "site-profile-result-card__title", status_label),
    tags$p(class = "site-profile-result-card__copy", paste(job_copy, collapse = " | "))
  )
}
