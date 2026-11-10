structure_export_defaults <- function() {
  list(
    format = "png",
    width = 3000,
    height = 1800,
    dpi = 300
  )
}

structure_export_config <- function(id) {
  ns <- NS(id)
  defaults <- structure_export_defaults()

  list(
    ids = list(
      format = ns("export_format"),
      width = ns("export_width"),
      height = ns("export_height"),
      dpi = ns("export_dpi"),
      trigger = ns("export_figure")
    ),
    defaults = defaults,
    sectioned = TRUE,
    choices = list(
      format = list(
        list(label = "PNG", value = "png"),
        list(label = "PDF", value = "pdf")
      )
    )
  )
}

structure_controls_defaults <- function() {
  list(
    sample = "1",
    coordinate = "",
    window_size = "fullsize",
    job_lookup = ""
  )
}

structure_has_saved_upload <- function(upload_context) {
  saved_files <- if (is.null(upload_context)) list() else upload_run_context_saved_files(upload_context)

  !is.null(upload_context) &&
    length(saved_files) > 0 &&
    all(vapply(saved_files, function(file) {
      !is.null(file$file_path) && nzchar(file$file_path) && file.exists(file$file_path)
    }, logical(1)))
}

structure_bed_table <- function(file_entry, limit = 250L) {
  path <- file_entry$file_path
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(data.table::data.table())
  }

  limit <- suppressWarnings(as.integer(limit))
  use_limit <- is.finite(limit) && limit > 0L
  bed <- if (use_limit) {
    data.table::fread(
      path,
      header = FALSE,
      sep = "\t",
      fill = TRUE,
      quote = "",
      data.table = TRUE,
      nrows = limit
    )
  } else {
    browser_read_bed_table(path)
  }
  if (!nrow(bed) || ncol(bed) < 3) {
    return(data.table::data.table())
  }

  if (use_limit) {
    bed <- bed[seq_len(min(nrow(bed), limit))]
  }
  data.table::data.table(
    row_index = seq_len(nrow(bed)),
    chrom = browser_canonical_chr_vector(bed[[1]]),
    start = suppressWarnings(as.integer(bed[[2]])),
    end = suppressWarnings(as.integer(bed[[3]])),
    name = if (ncol(bed) >= 4) as.character(bed[[4]]) else "",
    strand = if (ncol(bed) >= 6) as.character(bed[[6]]) else "."
  )[is.finite(start) & is.finite(end) & end > start & nzchar(chrom)]
}

structure_interval_label <- function(row) {
  label <- sprintf(
    "%s:%s-%s",
    row$chrom,
    format(row$start + 1L, scientific = FALSE, trim = TRUE),
    format(row$end, scientific = FALSE, trim = TRUE)
  )
  if (!is.null(row$name) && nzchar(row$name) && !identical(row$name, ".")) {
    label <- sprintf("%s | %s", row$name, label)
  }
  if (!is.null(row$strand) && row$strand %in% c("+", "-")) {
    label <- sprintf("%s (%s)", label, row$strand)
  }
  label
}

structure_sample_choices <- function(saved_files) {
  if (!length(saved_files)) {
    return(list())
  }

  lapply(seq_along(saved_files), function(index) {
    list(
      label = saved_files[[index]]$file_name %||% sprintf("BED file %d", index),
      value = as.character(index)
    )
  })
}

structure_interval_choices <- function(file_entry) {
  bed <- structure_bed_table(file_entry)
  if (!nrow(bed)) {
    return(list())
  }

  lapply(seq_len(nrow(bed)), function(index) {
    row <- bed[index]
    list(
      label = structure_interval_label(row),
      value = as.character(row$row_index)
    )
  })
}

structure_coordinate_placeholder <- function(file_entry, limit = 50L) {
  path <- file_entry$file_path
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return("chrom:start-end")
  }

  bed <- tryCatch(
    data.table::fread(
      path,
      header = FALSE,
      sep = "\t",
      fill = TRUE,
      quote = "",
      data.table = TRUE,
      nrows = limit
    ),
    error = function(...) data.table::data.table()
  )
  if (!nrow(bed) || ncol(bed) < 3) {
    return("chrom:start-end")
  }

  preview <- data.table::data.table(
    chrom = trimws(as.character(bed[[1]])),
    start = suppressWarnings(as.integer(bed[[2]])),
    end = suppressWarnings(as.integer(bed[[3]]))
  )
  preview <- preview[nzchar(chrom) & is.finite(start) & is.finite(end) & end > start]
  if (!nrow(preview)) {
    return("chrom:start-end")
  }

  sprintf("%s:%s-%s", preview$chrom[[1]], preview$start[[1]], preview$end[[1]])
}

structure_controls_config <- function(
  id,
  upload_context,
  selected_sample = 1L,
  coordinate = NULL,
  window_size = NULL,
  job_lookup = NULL,
  loaded = FALSE,
  active_job_state = NULL,
  queue_status = NULL
) {
  ns <- NS(id)
  defaults <- structure_controls_defaults()
  saved_files <- if (is.null(upload_context)) list() else upload_run_context_saved_files(upload_context)
  sample_index <- suppressWarnings(as.integer(selected_sample))
  if (!is.finite(sample_index)) {
    sample_index <- 1L
  }
  sample_index <- min(max(sample_index, 1L), max(1L, length(saved_files)))
  selected_file <- if (length(saved_files)) saved_files[[sample_index]] else NULL
  coordinate_placeholder <- if (is.null(selected_file)) "chrom:start-end" else structure_coordinate_placeholder(selected_file)
  queue_snapshot <- queue_status %||% list()
  queue_limit <- max(1L, suppressWarnings(as.integer(queue_snapshot$limit %||% STRUCTURE_JOB_ACTIVE_LIMIT)))
  queue_active <- max(0L, suppressWarnings(as.integer(queue_snapshot$active %||% 0L)))
  queue_queued <- max(0L, suppressWarnings(as.integer(queue_snapshot$queued %||% 0L)))
  queue_running <- max(0L, suppressWarnings(as.integer(queue_snapshot$running %||% 0L)))

  list(
    ids = list(
      sample = ns("structure_sample"),
      coordinate = ns("structure_coordinate"),
      window_size = ns("structure_window_size"),
      job_lookup = ns("job_lookup"),
      query_job = ns("query_job")
    ),
    defaults = utils::modifyList(defaults, list(
      sample = as.character(sample_index),
      coordinate = as.character(coordinate %||% defaults$coordinate),
      coordinate_placeholder = coordinate_placeholder,
      window_size = as.character(window_size %||% defaults$window_size),
      job_lookup = as.character(job_lookup %||% defaults$job_lookup)
    )),
    enabled = structure_has_saved_upload(upload_context),
    loaded = isTRUE(loaded),
    choices = list(
      samples = structure_sample_choices(saved_files),
      window_sizes = list(list(label = "Full size", value = "fullsize"))
    ),
    job = list(
      activeJobId = active_job_state$jobId %||% "",
      status = active_job_state$status %||% "idle",
      engineLabel = active_job_state$engineLabel %||% "",
      message = active_job_state$message %||% ""
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
