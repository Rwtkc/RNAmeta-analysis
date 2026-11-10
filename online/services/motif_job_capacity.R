motif_job_capacity_snapshot <- function(limit = MOTIF_JOB_ACTIVE_LIMIT) {
  states <- motif_job_list_states()
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
