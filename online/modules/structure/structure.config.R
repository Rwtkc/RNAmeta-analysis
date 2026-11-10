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
  job_lookup = NULL,
  loaded = FALSE,
  active_job_state = NULL,
  queue_status = NULL,
  mapping_status = "idle",
  mapping_message = "",
  analysis_locked = FALSE
) {
  ns <- NS(id)
  defaults <- structure_controls_defaults()
  queue_snapshot <- queue_status %||% list()
  queue_limit <- max(1L, suppressWarnings(as.integer(queue_snapshot$limit %||% STRUCTURE_JOB_ACTIVE_LIMIT)))
  queue_active <- max(0L, suppressWarnings(as.integer(queue_snapshot$active %||% 0L)))
  queue_queued <- max(0L, suppressWarnings(as.integer(queue_snapshot$queued %||% 0L)))
  queue_running <- max(0L, suppressWarnings(as.integer(queue_snapshot$running %||% 0L)))

  list(
    ids = list(
      job_lookup = ns("job_lookup"),
      query_job = ns("query_job"),
      map_structure = ns("map_structure")
    ),
    defaults = utils::modifyList(defaults, list(
      job_lookup = as.character(job_lookup %||% defaults$job_lookup)
    )),
    enabled = structure_has_saved_upload(upload_context),
    analysisLocked = isTRUE(analysis_locked),
    loaded = isTRUE(loaded),
    mapping = list(
      status = as.character(mapping_status %||% "idle"),
      message = as.character(mapping_message %||% "")
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

structure_react_mapping_table_payload <- function(table, page = 1L, page_size = 10L, query = "") {
  empty_payload <- list(columns = list(), rows = list(), total = 0L)
  if (is.null(table) || !nrow(table)) {
    return(empty_payload)
  }

  bed_cols <- grep("^bed_col_[0-9]+$", names(table), value = TRUE)
  bed_cols <- bed_cols[order(as.integer(sub("^bed_col_", "", bed_cols)))]
  first_bed_cols <- bed_cols[seq_len(min(3L, length(bed_cols)))]
  remaining_bed_cols <- if (length(bed_cols) > 3L) bed_cols[4:length(bed_cols)] else character()
  display_cols <- c(
    first_bed_cols,
    "gene_id",
    "transcript_id",
    "transcript_interval",
    "submit_rnafold",
    remaining_bed_cols
  )
  display_cols <- display_cols[display_cols %in% c(names(table), "submit_rnafold")]
  query <- trimws(as.character(query %||% ""))[[1]]
  searchable_cols <- intersect(setdiff(display_cols, "submit_rnafold"), names(table))
  if (nzchar(query) && length(searchable_cols)) {
    matches <- Reduce(
      `|`,
      lapply(searchable_cols, function(column) {
        grepl(query, as.character(table[[column]]), fixed = TRUE, ignore.case = TRUE)
      }),
      init = rep(FALSE, nrow(table))
    )
    table <- table[matches]
  }

  total_rows <- nrow(table)
  page_size <- max(1L, suppressWarnings(as.integer(page_size)))
  page <- max(1L, suppressWarnings(as.integer(page)))
  total_pages <- max(1L, ceiling(total_rows / page_size))
  page <- min(page, total_pages)
  page_start <- if (total_rows) ((page - 1L) * page_size) + 1L else 1L
  page_end <- if (total_rows) min(total_rows, page_start + page_size - 1L) else 0L
  page_table <- if (total_rows) table[page_start:page_end] else table[0]
  bed_headers <- browser_bed_column_names(length(bed_cols))
  headers <- c(
    bed_headers[seq_along(first_bed_cols)],
    "Gene ID",
    "Transcript ID",
    "Transcript interval",
    "Submit to RNAfold",
    if (length(remaining_bed_cols)) bed_headers[seq(4L, length(bed_cols))] else character()
  )
  headers <- headers[seq_along(display_cols)]

  columns <- lapply(seq_along(display_cols), function(index) {
    list(id = display_cols[[index]], label = headers[[index]])
  })
  rows <- lapply(seq_len(nrow(page_table)), function(index) {
    row <- page_table[index]
    values <- lapply(display_cols, function(column) {
      if (identical(column, "submit_rnafold")) {
        return("Submit to RNAfold")
      }
      value <- row[[column]][[1]] %||% ""
      if (length(value) != 1L || is.na(value)) "" else as.character(value)
    })
    names(values) <- display_cols
    values$row_key <- as.character(row$row_key[[1]] %||% "")
    values$can_submit <- nzchar(as.character(row$transcript_id[[1]] %||% ""))
    values
  })

  list(
    columns = columns,
    rows = rows,
    total = as.integer(total_rows),
    page = as.integer(page),
    pageSize = as.integer(page_size),
    query = query
  )
}

structure_mapping_config <- function(
  id,
  upload_context,
  table,
  status = "idle",
  message = "",
  revision = 0L,
  submission_revision = 0L,
  active_job_state = NULL,
  result = NULL,
  page = 1L,
  query = ""
) {
  ns <- NS(id)
  list(
    hasSavedUpload = structure_has_saved_upload(upload_context),
    workspaceStatus = as.character(status %||% "idle"),
    message = as.character(message %||% ""),
    table = structure_react_mapping_table_payload(table, page = page, query = query),
    page = as.integer(page),
    query = as.character(query %||% ""),
    tableRevision = as.integer(revision),
    submissionRevision = as.integer(submission_revision),
    renderedInput = ns("mapping_rendered"),
    renderToken = as.character(revision),
    submitInput = ns("structure_submit_mapping"),
    pageInput = ns("mapping_page"),
    queryInput = ns("mapping_query"),
    job = list(
      activeJobId = active_job_state$jobId %||% "",
      status = active_job_state$status %||% "",
      engineLabel = active_job_state$engineLabel %||% "RNAfold",
      message = active_job_state$message %||% "",
      errorMessage = active_job_state$errorMessage %||% "",
      method = result$method %||% "",
      locus = result$locus %||% ""
    )
  )
}
