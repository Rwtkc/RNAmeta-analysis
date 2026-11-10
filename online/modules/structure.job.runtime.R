structure_send_render <- function(session, result_payload) {
  if (is.null(result_payload) ||
      is.null(result_payload$sequence) ||
      is.null(result_payload$structure)) {
    return(invisible(FALSE))
  }

  session$sendCustomMessage(
    "rnameta-structure-render",
    list(
      frameId = session$ns("forna_frame"),
      sequence = result_payload$sequence,
      structure = result_payload$structure,
      highlights = result_payload$highlightRanges %||% list(),
      highlightStart = result_payload$highlightStart %||% NULL,
      highlightEnd = result_payload$highlightEnd %||% NULL
    )
  )
}

structure_job_notice_card <- function(job_state = NULL, error_message = NULL) {
  if (!is.null(error_message) && nzchar(error_message)) {
    return(
      div(
        class = "site-profile-result-card site-profile-result-card--error structure-result-card",
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
    class = "site-profile-result-card structure-result-card",
    tags$h4(class = "site-profile-result-card__title", status_label),
    tags$p(class = "site-profile-result-card__copy", paste(job_copy, collapse = " | "))
  )
}

structure_highlight_sequence <- function(result) {
  sequence <- as.character(result$sequence %||% "")
  start <- suppressWarnings(as.integer(result$highlightStart %||% NA_integer_))
  end <- suppressWarnings(as.integer(result$highlightEnd %||% NA_integer_))

  if (!nzchar(sequence) || !is.finite(start) || !is.finite(end) || end < start) {
    return("")
  }

  start <- max(1L, start)
  end <- min(nchar(sequence), end)
  if (end < start) {
    return("")
  }

  substr(sequence, start, end)
}

structure_basic_info_rows <- function(result) {
  fragment_sequence <- structure_highlight_sequence(result)
  fragment_length <- nchar(fragment_sequence)
  if (!nzchar(fragment_sequence)) {
    fragment_length <- suppressWarnings(as.integer(result$overlapBases %||% 0L))
  }

  list(
    list(label = "Fragment Type", value = "BED-derived RNA fragment"),
    list(label = "Source Transcript", value = result$transcriptId %||% ""),
    list(label = "Source Gene", value = result$geneId %||% ""),
    list(label = "Transcript Length", value = sprintf("%s nt", result$transcriptLength %||% result$sequenceLength %||% "")),
    list(label = "Fragment Length", value = sprintf("%s nt", fragment_length)),
    list(label = "Locus on Transcript", value = sprintf("%s-%s", result$highlightStart %||% "", result$highlightEnd %||% "")),
    list(label = "Sequence", value = fragment_sequence, type = "sequence")
  )
}

structure_sequence_preview <- function(sequence, width = 100L) {
  sequence <- as.character(sequence %||% "")
  width <- suppressWarnings(as.integer(width))
  if (!is.finite(width) || width <= 0L) {
    width <- 100L
  }

  if (nchar(sequence) <= width) {
    return(sequence)
  }

  paste0(substr(sequence, 1L, width), "...")
}

structure_basic_info_value_tag <- function(row) {
  value <- as.character(row$value %||% "")
  if (identical(row$type %||% "", "sequence") && nchar(value) > 100L) {
    return(
      tags$details(
        class = "structure-sequence-toggle",
        tags$summary(
          tags$span(class = "structure-sequence-toggle__preview", structure_sequence_preview(value, 100L)),
          tags$span(class = "structure-sequence-toggle__more", "Show full sequence")
        ),
        tags$p(class = "structure-sequence-toggle__full", value)
      )
    )
  }

  value
}

structure_basic_info_card <- function(result) {
  rows <- structure_basic_info_rows(result)
  rows <- rows[vapply(rows, function(row) nzchar(as.character(row$value %||% "")), logical(1))]

  div(
    class = "structure-basic-info",
    tags$h3(class = "structure-basic-info__title", "Basic information of RNA fragment"),
    tags$dl(
      class = "structure-basic-info__grid",
      lapply(rows, function(row) {
        tagList(
          tags$dt(row$label),
          tags$dd(structure_basic_info_value_tag(row))
        )
      })
    )
  )
}

structure_selected_saved_file <- function(upload_context, sample_index) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (!length(saved_files)) {
    return(NULL)
  }

  sample_index <- suppressWarnings(as.integer(sample_index))
  if (!is.finite(sample_index)) {
    sample_index <- 1L
  }
  sample_index <- min(max(sample_index, 1L), length(saved_files))
  saved_files[[sample_index]]
}

structure_job_input_filename <- function(file_name) {
  cleaned <- gsub("[^A-Za-z0-9._-]+", "_", basename(as.character(file_name %||% "input.bed")[[1]]))
  if (!nzchar(cleaned)) {
    cleaned <- "input.bed"
  }

  cleaned
}

structure_job_worker_upload_context <- function(upload_context, selected_file, job_dir) {
  input_dir <- file.path(job_dir, "inputs")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

  target_path <- file.path(input_dir, structure_job_input_filename(selected_file$file_name))
  if (!file.copy(selected_file$file_path, target_path, overwrite = TRUE)) {
    stop("Failed to copy selected BED file into the Structure job workspace.")
  }

  copied_file <- selected_file
  copied_file$file_path <- normalizePath(target_path, winslash = "/", mustWork = TRUE)

  worker_context <- upload_context
  worker_context$file_count <- 1L
  worker_context$file_name <- copied_file$file_name
  worker_context$file_path <- copied_file$file_path
  worker_context$file_names <- copied_file$file_name
  worker_context$saved_files <- list(copied_file)
  worker_context$preview <- NULL
  worker_context
}

structure_job_rscript_path <- function() {
  configured_path <- "/public/liuqi/miniconda3/envs/R4.3.0/bin/Rscript"
  if (file.exists(configured_path)) {
    return(configured_path)
  }

  fallback_path <- file.path(R.home("bin"), "Rscript")
  if (file.exists(fallback_path)) {
    return(fallback_path)
  }

  stop("Unable to locate Rscript for Structure background jobs.")
}

structure_job_launch <- function(job_id) {
  rscript_path <- structure_job_rscript_path()
  worker_script <- app_path("services", "structure_job_worker.R")
  worker_stdout <- structure_job_worker_stdout_path(job_id)
  worker_stderr <- structure_job_worker_stderr_path(job_id)
  worker_args <- c(worker_script, job_id, app_root)

  dir.create(structure_job_dir(job_id), recursive = TRUE, showWarnings = FALSE)
  unlink(worker_stdout, force = TRUE)
  unlink(worker_stderr, force = TRUE)

  if (.Platform$OS.type == "windows") {
    argument_list <- paste(
      vapply(c(worker_script, job_id, app_root), shQuote, character(1)),
      collapse = ", "
    )
    command <- sprintf(
      "Start-Process -FilePath %s -ArgumentList @(%s) -WorkingDirectory %s -RedirectStandardOutput %s -RedirectStandardError %s -WindowStyle Hidden",
      shQuote(rscript_path),
      argument_list,
      shQuote(app_root),
      shQuote(worker_stdout),
      shQuote(worker_stderr)
    )
    invisible(system2("powershell", c("-NoProfile", "-Command", command), wait = FALSE))
    return(invisible(TRUE))
  }

  process_id <- suppressWarnings(
    system2(
      command = rscript_path,
      args = worker_args,
      stdout = worker_stdout,
      stderr = worker_stderr,
      wait = FALSE
    )
  )

  is.numeric(process_id) || is.integer(process_id)
}

structure_job_submit_state <- function(upload_context, controls, selected_file, job_id) {
  list(
    jobId = job_id,
    status = "queued",
    engine = "rnafold",
    engineLabel = "RNAfold",
    sampleName = basename(selected_file$file_name),
    sourceFile = selected_file$file_name,
    species = upload_context$species,
    coordinate = trimws(as.character(controls$coordinate %||% "")),
    windowSize = as.character(controls$window_size %||% structure_controls_defaults()$window_size),
    submittedAt = structure_job_timestamp(),
    startedAt = NULL,
    finishedAt = NULL,
    expiresAt = structure_job_expires_at(),
    message = "Queued RNAfold Structure analysis.",
    workspaceDir = normalizePath(structure_job_dir(job_id), winslash = "/", mustWork = FALSE),
    resultDir = NULL,
    errorMessage = NULL
  )
}

structure_submit_background_job <- function(upload_context, controls) {
  structure_job_with_submission_lock({
    controls <- utils::modifyList(structure_controls_defaults(), controls)
    if (!structure_has_saved_upload(upload_context)) {
      stop("Save a BED file in Upload / Run before running Structure.")
    }
    if (!nzchar(trimws(as.character(controls$coordinate %||% "")))) {
      stop("Enter Coordinate before running Structure.")
    }

    structure_job_assert_capacity()

    selected_file <- structure_selected_saved_file(upload_context, controls$sample)
    if (is.null(selected_file) || !nzchar(selected_file$file_path) || !file.exists(selected_file$file_path)) {
      stop("No saved BED file is available for Structure analysis.")
    }
    structure_resolve_coordinate_selection(
      upload_context = upload_context,
      sample_index = controls$sample,
      coordinate = controls$coordinate
    )

    job_id <- structure_job_new_id()
    job_dir <- structure_job_dir(job_id)
    dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
    worker_upload_context <- structure_job_worker_upload_context(upload_context, selected_file, job_dir)

    state <- structure_job_submit_state(upload_context, controls, selected_file, job_id)
    spec <- list(
      job_id = job_id,
      upload_context = worker_upload_context,
      controls = controls
    )

    structure_job_atomic_save_rds(structure_job_spec_path(job_id), spec)
    structure_job_save_state(job_id, state)
    structure_job_launch(job_id)

    list(job_id = job_id, state = state)
  })
}

structure_job_load_spec <- function(job_id) {
  spec_path <- structure_job_spec_path(job_id)
  if (!file.exists(spec_path)) {
    stop(sprintf("Structure job spec for '%s' is missing.", job_id))
  }

  readRDS(spec_path)
}

run_structure_analysis <- function(upload_context, controls) {
  selection <- structure_resolve_coordinate_selection(
    upload_context = upload_context,
    sample_index = controls$sample,
    coordinate = controls$coordinate
  )
  extracted <- structure_extract_sequence(
    upload_context = upload_context,
    selection = selection,
    window_size = controls$window_size
  )
  folded <- structure_fold_sequence(extracted$sequence)

  structure_result_payload(list(selection = selection, extracted = extracted, folded = folded))
}

structure_job_result_state_patch <- function(result) {
  list(
    status = "completed",
    finishedAt = structure_job_timestamp(),
    expiresAt = structure_job_expires_at(),
    message = "Structure preview completed.",
    resultDir = result$workspaceDir %||% NULL,
    errorMessage = NULL,
    fileName = result$fileName %||% "",
    coordinate = result$coordinate %||% "",
    locus = result$locus %||% "",
    strand = result$strand %||% "",
    sequenceLength = result$sequenceLength %||% 0L,
    method = result$method %||% "",
    windowSize = result$windowSize %||% "",
    intervalLength = result$intervalLength %||% 0L,
    usedBedInterval = isTRUE(result$usedBedInterval),
    matchedBedRow = isTRUE(result$matchedBedRow),
    transcriptId = result$transcriptId %||% "",
    geneId = result$geneId %||% "",
    bedLocus = result$bedLocus %||% "",
    highlightStart = result$highlightStart %||% NULL,
    highlightEnd = result$highlightEnd %||% NULL
  )
}

structure_job_failure_state_patch <- function(message, started_at = NULL) {
  list(
    status = "failed",
    startedAt = started_at,
    finishedAt = structure_job_timestamp(),
    expiresAt = structure_job_expires_at(),
    message = "Structure analysis failed.",
    errorMessage = as.character(message)[[1]]
  )
}

structure_background_worker_run <- function(job_id) {
  spec <- structure_job_load_spec(job_id)
  started_at <- structure_job_timestamp()
  structure_job_update_state(job_id, list(
    status = "running",
    startedAt = started_at,
    message = "Running RNAfold on the server workspace.",
    errorMessage = NULL
  ))

  result <- tryCatch(
    run_structure_analysis(
      upload_context = spec$upload_context,
      controls = spec$controls
    ),
    error = function(err) {
      structure_job_update_state(job_id, structure_job_failure_state_patch(conditionMessage(err), started_at = started_at))
      NULL
    }
  )

  if (is.null(result)) {
    return(invisible(NULL))
  }

  result$jobId <- job_id
  result$startedAt <- started_at
  result$finishedAt <- structure_job_timestamp()
  result$workspaceDir <- normalizePath(structure_job_dir(job_id), winslash = "/", mustWork = FALSE)
  structure_job_atomic_save_rds(structure_job_result_path(job_id), result)
  structure_job_update_state(job_id, structure_job_result_state_patch(result))
  invisible(TRUE)
}
