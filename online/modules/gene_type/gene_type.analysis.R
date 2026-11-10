gene_type_compute_counts <- function(peaks, feature_ranges, sample_name) {
  transcript_type <- gene_type_find_overlapping_type_values(peaks, feature_ranges)

  if (!length(transcript_type)) {
    return(data.table::data.table(
      GeneType = character(),
      Frequency = numeric(),
      Sample = character()
    ))
  }

  counts <- data.table::as.data.table(table(transcript_type))
  data.table::setnames(counts, c("GeneType", "Frequency"))
  counts[, Frequency := as.numeric(Frequency)]
  counts[, Sample := sample_name]
  counts[order(-Frequency, GeneType)]
}

run_gene_type_analysis <- function(upload_context, controls = gene_type_default_controls(), progress = NULL) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (is.null(upload_context) || !length(saved_files)) {
    stop("Save a BED file in Upload / Run before running Gene Type.")
  }

  if (is.null(upload_context$species) || !nzchar(upload_context$species)) {
    stop("A species must be saved in Upload / Run before running Gene Type.")
  }

  controls <- utils::modifyList(gene_type_default_controls(), controls)
  if (exists("gene_type_cached_analysis_result", mode = "function")) {
    cached_result <- gene_type_cached_analysis_result(upload_context)
    if (!is.null(cached_result)) {
      gene_type_progress_update(progress, 100, "Reusing cached Gene Type result for the saved BED context")
      return(cached_result)
    }
  }

  gene_type_progress_update(progress, 12, "Loading saved BED files and species annotation bundle")
  annotation_bundle <- gene_type_load_annotation_bundle(upload_context$species)
  gene_type_progress_update(progress, 28, "Preparing cached transcript biotype feature ranges")
  feature_ranges <- gene_type_feature_ranges(upload_context$species, annotation_bundle)
  gene_type_counts_list <- vector("list", length(saved_files))
  sample_metadata <- rnameta_build_sample_metadata(saved_files)
  interval_count <- 0L

  for (index in seq_along(saved_files)) {
    file_info <- saved_files[[index]]
    sample_name <- sample_metadata[[index]]$display_name
    gene_type_progress_update(
      progress,
      40 + floor(index / length(saved_files) * 40),
      sprintf("Counting transcript biotypes for %s (%d/%d)", basename(file_info$file_name), index, length(saved_files))
    )
    analysis_peaks <- gene_type_read_peaks(file_info$file_path, annotation_bundle$txdb, sample_name = sample_name)
    gene_type_counts_list[[index]] <- gene_type_compute_counts(analysis_peaks, feature_ranges, sample_name)
    interval_count <- interval_count + length(analysis_peaks)
  }
  gene_type_counts <- data.table::rbindlist(gene_type_counts_list, use.names = TRUE, fill = TRUE)
  gene_type_progress_update(progress, 100, "Building chart payload and finalizing results")

  result <- list(
    chart_payload = gene_type_build_chart_payload(gene_type_counts, sample_metadata = sample_metadata),
    gene_type_counts = gene_type_counts,
    species = upload_context$species,
    sample_name = rnameta_sample_original_names(sample_metadata),
    sample_display_names = rnameta_sample_display_names(sample_metadata),
    sample_count = length(saved_files),
    interval_count = interval_count,
    overlap_count = sum(gene_type_counts$Frequency, na.rm = TRUE)
  )

  if (exists("gene_type_store_analysis_result", mode = "function")) {
    gene_type_store_analysis_result(upload_context, result)
  }

  result
}
