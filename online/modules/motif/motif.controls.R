motif_controls_config <- function(id, upload_context = NULL, current_values = list(), active_job_state = NULL, queue_status = NULL) {
  ns <- NS(id)
  defaults <- utils::modifyList(motif_default_controls(), current_values)
  sample_choices <- motif_build_sample_choices(upload_context)
  sample_values <- unname(sample_choices)
  engine_values <- unname(motif_engine_choices())
  subset_values <- as.character(motif_sites_subset_choices())
  sample_name <- defaults$sample_name
  if (is.null(sample_name)) {
    sample_name <- ""
  }
  sample_name <- as.character(sample_name)[[1]]
  if (!length(sample_values)) {
    sample_name <- ""
  } else if (!nzchar(sample_name) || !sample_name %in% sample_values) {
    sample_name <- sample_values[[1]]
  }

  engine <- defaults$engine
  if (is.null(engine)) {
    engine <- motif_default_controls()$engine
  }
  engine <- as.character(engine)[[1]]
  if (!engine %in% engine_values) {
    engine <- motif_default_controls()$engine
  }

  sites_subset <- defaults$sites_subset
  if (is.null(sites_subset)) {
    sites_subset <- motif_default_controls()$sites_subset
  }
  sites_subset <- motif_normalize_subset_value(sites_subset)

  motif_width <- motif_effective_width(engine, defaults$motif_width)
  max_motifs <- max(1L, motif_safe_integer(defaults$max_motifs, motif_default_controls()$max_motifs))
  job_lookup <- as.character(defaults$job_lookup %||% "")[[1]]
  queue_snapshot <- queue_status %||% list()
  queue_limit <- max(1L, motif_safe_integer(queue_snapshot$limit, MOTIF_JOB_ACTIVE_LIMIT))
  queue_active <- max(0L, motif_safe_integer(queue_snapshot$active, 0L))
  queue_queued <- max(0L, motif_safe_integer(queue_snapshot$queued, 0L))
  queue_running <- max(0L, motif_safe_integer(queue_snapshot$running, 0L))

  list(
    ids = list(
      sample_name = ns("sample_name"),
      engine = ns("engine"),
      analysis_mode = ns("analysis_mode"),
      sites_subset = ns("sites_subset"),
      motif_width = ns("motif_width"),
      sequence_window = ns("sequence_window"),
      max_motifs = ns("max_motifs"),
      result_bundle = ns("result_bundle"),
      job_lookup = ns("job_lookup"),
      query_job = ns("query_job")
    ),
    defaults = list(
      sample_name = sample_name,
      engine = engine,
      analysis_mode = motif_validate_analysis_mode(defaults$analysis_mode),
      sites_subset = sites_subset,
      motif_width = motif_width,
      sequence_window = defaults$sequence_window %||% motif_default_controls()$sequence_window,
      max_motifs = max_motifs,
      job_lookup = job_lookup
    ),
    choices = list(
      samples = motif_choice_options(sample_choices),
      engines = motif_choice_options(motif_engine_choices()),
      analysis_modes = motif_choice_options(motif_analysis_mode_choices()),
      sites_subset = motif_choice_options(motif_sites_subset_choices())
    ),
    constraints = list(
      defaultMinWidth = motif_min_width("homer"),
      homerMinWidth = motif_min_width("homer"),
      memeMinWidth = motif_min_width("meme"),
      stremeMinWidth = motif_min_width("streme")
    ),
    job = list(
      activeJobId = if (is.null(active_job_state$jobId)) "" else active_job_state$jobId,
      status = if (is.null(active_job_state$status)) "idle" else active_job_state$status,
      engineLabel = if (is.null(active_job_state$engineLabel)) "" else active_job_state$engineLabel,
      message = if (is.null(active_job_state$message)) "" else active_job_state$message
    ),
    queue = list(
      active = queue_active,
      queued = queue_queued,
      running = queue_running,
      limit = queue_limit,
      available = max(0L, queue_limit - queue_active)
    )
  )
}
