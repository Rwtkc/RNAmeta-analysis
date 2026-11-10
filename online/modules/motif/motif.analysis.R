motif_initialize_result <- function(upload_context, selected_file, controls, motif_peaks, job_id = NULL) {
  engine_key <- tolower(as.character(controls$engine)[[1]])
  engine_label <- motif_engine_label(engine_key)
  subset_label <- motif_subset_label(controls$sites_subset)
  motif_width <- motif_effective_width(engine_key, controls$motif_width)
  analysis_mode <- motif_validate_analysis_mode(controls$analysis_mode)
  sequence_window <- if (identical(analysis_mode, "matched_full")) NULL else motif_validate_sequence_window(controls$sequence_window, motif_width)
  max_motifs <- max(1L, motif_safe_integer(controls$max_motifs, 6L))
  command_preview <- motif_engine_command_preview(engine_key, basename(selected_file$file_name), list(
    motif_width = motif_width,
    max_motifs = max_motifs,
    analysis_mode = analysis_mode
  ))

  list(
    status = "prepared",
    job_id = job_id,
    job_status = "prepared",
    job_message = sprintf("%s workspace prepared.", engine_label),
    started_at = motif_job_timestamp(),
    finished_at = NULL,
    engine = engine_key,
    analysis_mode = analysis_mode,
    engine_label = engine_label,
    sample_name = basename(selected_file$file_name),
    source_file = selected_file$file_name,
    source_path = selected_file$file_path,
    species = upload_context$species,
    sites_subset = subset_label,
    motif_width = motif_width,
    max_motifs = max_motifs,
    peak_count = length(motif_peaks),
    sequence_window = sequence_window,
    query_count = length(motif_peaks),
    control_count = 0L,
    bundle_name = "",
    bundle_status = "server-run",
    bundle_message = sprintf("%s results will be parsed directly from the server workspace.", engine_label),
    command_preview = command_preview,
    expected_outputs = unique(c(motif_engine_expected_outputs(engine_key),
      if (identical(analysis_mode, "matched_full")) "control.fa")),
    motifs = list(),
    notes = c(
      sprintf("Selected engine: %s", engine_label),
      sprintf("Selected sample: %s", basename(selected_file$file_name)),
      sprintf("Peak subset rule: %s", subset_label),
      if (identical(analysis_mode, "matched_full")) "Complete BED intervals; RNA strand; matched non-peak background" else sprintf("Peak-centred sequence window: %d nt", sequence_window),
      sprintf("Peak rows entering motif stage: %d", length(motif_peaks))
    )
  )
}

motif_execute_analysis <- function(upload_context, controls = motif_default_controls(), progress = NULL, job_id = NULL) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (!length(saved_files)) {
    stop("Save a BED file in Upload / Run before running motif analysis.")
  }

  controls <- motif_validate_analysis_controls(controls)
  selected_file <- motif_selected_saved_file(upload_context, controls$sample_name)
  if (is.null(selected_file)) {
    stop("No saved BED file is available for motif analysis.")
  }

  motif_progress_update(progress, 16, sprintf("Preparing motif input for %s", basename(selected_file$file_name)))
  annotation_bundle <- peak_distribution_load_annotation_bundle(upload_context$species)
  raw_peaks <- peak_distribution_read_peaks(
    source_file = selected_file$file_path,
    txdb = annotation_bundle$txdb,
    sample_name = basename(selected_file$file_name),
    subset_for_motif = "All"
  )

  motif_progress_update(progress, 42, "Applying motif subset selection and engine parameters")
  motif_peaks <- motif_subset_peak_ranges(raw_peaks, controls$sites_subset)
  result <- motif_initialize_result(
    upload_context = upload_context,
    selected_file = selected_file,
    controls = controls,
    motif_peaks = motif_peaks,
    job_id = job_id
  )

  motif_progress_update(progress, 72, sprintf("Preparing %s workspace and command", result$engine_label))
  motif_progress_update(progress, 82, "Generating server-side motif input package")
  exclusion <- raw_peaks
  if (identical(result$analysis_mode, "matched_full")) {
    all_peaks <- lapply(saved_files, function(file) {
      if (identical(file$file_path, selected_file$file_path)) return(raw_peaks)
      peak_distribution_read_peaks(file$file_path, annotation_bundle$txdb, file$file_name, "All")
    })
    exclusion <- do.call(c, unname(all_peaks))
  }
  result$source_input_md5 <- as.list(tools::md5sum(vapply(
    if (identical(result$analysis_mode, "matched_full")) saved_files else list(selected_file),
    function(file) file$file_path, character(1))))
  job_package <- motif_prepare_job_package(result, motif_peaks, annotation_bundle, exclusion, progress)
  result$job_package_dir <- job_package$package_dir
  result$job_package_files <- job_package$package_files
  result$job_script_name <- job_package$script_name
  result$job_manifest_name <- job_package$manifest_name
  result$control_count <- job_package$control_count
  result$query_count <- job_package$query_count
  result$command_preview <- job_package$command_preview
  result$matching_summary <- job_package$manifest$matching
  result$analysis_metadata <- job_package$manifest
  result$notes <- c(
    result$notes,
    sprintf("Server motif workspace: %s", result$job_package_dir),
    sprintf("Generated workspace files: %s", paste(result$job_package_files, collapse = ", ")),
    sprintf("Engine command: %s", result$command_preview)
  )

  engine_run <- motif_run_engine(
    engine = result$engine,
    package_dir = result$job_package_dir,
    controls = list(
      motif_width = result$motif_width,
      max_motifs = result$max_motifs,
      analysis_mode = result$analysis_mode
    ),
    progress = progress
  )

  motif_progress_update(progress, 96, sprintf("Parsing %s result files", result$engine_label))
  parsed_result <- motif_detect_local_result(
    engine = result$engine,
    job_package_dir = result$job_package_dir,
    max_motifs = result$max_motifs
  )

  if (is.null(parsed_result) || !length(parsed_result$motifs)) {
    result$notes <- c(
      result$notes,
      sprintf("No parseable %s motifs were detected in output directory: %s", result$engine_label, engine_run$output_dir)
    )
    stop(sprintf(
      "%s",
      motif_result_parse_failure_message(
        engine_label = result$engine_label,
        engine = result$engine,
        motif_width = result$motif_width
      )
    ))
  }

  result$status <- "completed"
  result$job_status <- "completed"
  result$job_message <- sprintf("%s completed successfully.", result$engine_label)
  result$bundle_name <- parsed_result$result_name
  result$bundle_status <- parsed_result$result_status
  result$bundle_message <- parsed_result$result_message
  result$motifs <- parsed_result$motifs
  result$result_output_dir <- parsed_result$output_dir
  result$finished_at <- motif_job_timestamp()
  result$notes <- c(
    result$notes,
    sprintf("Parsed %s output directory: %s", result$engine_label, parsed_result$output_dir)
  )

  if (!is.null(engine_run$stdout_log) && !is.null(engine_run$stderr_log)) {
    result$notes <- c(
      result$notes,
      sprintf("Engine stdout log: %s", engine_run$stdout_log),
      sprintf("Engine stderr log: %s", engine_run$stderr_log)
    )
  }

  motif_progress_update(progress, 100, "Finalizing motif payload")
  result$chart_payload <- motif_build_chart_payload(result)
  result
}

run_motif_analysis <- function(upload_context, controls = motif_default_controls(), result_bundle = NULL, progress = NULL, job_id = NULL) {
  motif_execute_analysis(
    upload_context = upload_context,
    controls = controls,
    progress = progress,
    job_id = job_id
  )
}
