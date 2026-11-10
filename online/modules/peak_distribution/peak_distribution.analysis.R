peak_distribution_sample_label <- function(file_info) {
  if (!is.null(file_info$file_name) && nzchar(file_info$file_name)) {
    return(tools::file_path_sans_ext(basename(file_info$file_name)))
  }

  "Group1"
}

peak_distribution_build_series_metadata <- function(saved_files) {
  original_names <- vapply(saved_files, peak_distribution_sample_label, character(1))
  display_names <- sprintf("Sample %d", seq_along(original_names))

  lapply(seq_along(saved_files), function(index) {
    list(
      display_name = display_names[[index]],
      original_name = original_names[[index]]
    )
  })
}

run_peak_distribution_analysis <- function(upload_context, controls = peak_distribution_default_controls(), progress = NULL) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (is.null(upload_context) || !length(saved_files)) {
    stop("Save a BED file in Upload / Run before running Peak Distribution.")
  }

  if (is.null(upload_context$species) || !nzchar(upload_context$species)) {
    stop("A species must be saved in Upload / Run before running Peak Distribution.")
  }

  controls <- utils::modifyList(peak_distribution_default_controls(), controls)
  peak_distribution_progress_update(progress, 10, "Loading saved BED files and species annotation bundle")
  annotation_bundle <- peak_distribution_cached_annotation_bundle(upload_context$species)
  peak_distribution_progress_update(progress, 30, "Constructing feature-level transcript annotation structure")
  txdb_features <- peak_distribution_cached_txdb_features(upload_context$species, minimal_component_length = 100)
  feature_counts_list <- vector("list", length(saved_files))
  annotation_results <- vector("list", length(saved_files))
  detail_files <- character()
  series_metadata <- peak_distribution_build_series_metadata(saved_files)
  sample_labels <- vapply(series_metadata, function(entry) entry$display_name, character(1))
  interval_count <- 0L

  for (index in seq_along(saved_files)) {
    file_info <- saved_files[[index]]
    sample_name <- sprintf("Group%d", index)
    peak_distribution_progress_update(
      progress,
      45 + floor(index / length(saved_files) * 30),
      sprintf("Annotating %s (%d/%d)", basename(file_info$file_name), index, length(saved_files))
    )
    raw_peaks <- peak_distribution_read_peaks(file_info$file_path, annotation_bundle$txdb, sample_name = sample_name)
    analysis_peaks <- GenomicRanges::resize(raw_peaks, width = 1, fix = "center")
    annotation_result <- peak_distribution_annotate_peaks(
      peaks = analysis_peaks,
      txdb_features = txdb_features,
      source_file = file_info$file_path,
      detail_output_dir = controls$detail_output_dir,
      sample_index = index
    )
    feature_counts <- peak_distribution_compute_feature_counts(annotation_result)
    feature_counts[, Sample := sample_name]
    feature_counts_list[[index]] <- feature_counts
    annotation_results[[index]] <- annotation_result
    detail_files <- c(detail_files, annotation_result$detail_file)
    interval_count <- interval_count + length(raw_peaks)
  }

  peak_distribution_progress_update(progress, 88, "Summarizing feature counts and preparing grouped bar output")
  feature_counts <- data.table::rbindlist(feature_counts_list, use.names = TRUE, fill = TRUE)
  plot_data <- peak_distribution_prepare_plot_data(
    feature_counts = feature_counts,
    selected_features = controls$selected_features,
    sample_labels = sample_labels
  )

  peak_distribution_progress_update(progress, 100, "Building chart payload and finalizing results")
  list(
    plot = peak_distribution_build_plot(
      plot_data = plot_data
    ),
    chart_payload = peak_distribution_build_chart_payload(
      plot_data = plot_data,
      series_metadata = series_metadata
    ),
    feature_counts = plot_data,
    annotation_result = annotation_results,
    detail_file = detail_files,
    species = upload_context$species,
    sample_name = vapply(series_metadata, function(entry) entry$original_name, character(1)),
    sample_display_names = sample_labels,
    sample_count = length(saved_files),
    interval_count = interval_count
  )
}
