peak_exon_sample_label <- function(file_info) {
  rnameta_sample_original_name(file_info)
}

peak_exon_group_palette <- function() {
  c("#de7a5b", "#859b7a", "#6676a0", "#c9886d", "#b9c78f")
}

peak_exon_location_levels <- function() {
  c("First", "Middle", "Last")
}

peak_exon_progress_set <- function(progress, module_label, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  progress$set(
    value = value,
    message = sprintf("Running %s Analysis", module_label),
    detail = sprintf("%s%% | %s", value, detail)
  )

  invisible(NULL)
}

peak_exon_require_saved_upload <- function(upload_context, module_label) {
  saved_files <- upload_run_context_saved_files(upload_context)

  if (is.null(upload_context) || !length(saved_files)) {
    stop(sprintf("Save a BED file in Upload / Run before running %s.", module_label))
  }

  if (is.null(upload_context$species) || !nzchar(upload_context$species)) {
    stop(sprintf("A species must be saved in Upload / Run before running %s.", module_label))
  }

  saved_files
}

peak_exon_load_annotation_bundle <- function(species) {
  if (exists("peak_distribution_cached_annotation_bundle", mode = "function")) {
    return(peak_distribution_cached_annotation_bundle(species))
  }

  peak_distribution_load_annotation_bundle(species)
}

peak_exon_load_txdb_features <- function(species, minimal_component_length = 100) {
  if (exists("peak_distribution_cached_txdb_features", mode = "function")) {
    return(peak_distribution_cached_txdb_features(species, minimal_component_length = minimal_component_length))
  }

  annotation_bundle <- peak_exon_load_annotation_bundle(species)

  peak_distribution_build_txdb_features(
    txdb = annotation_bundle$txdb,
    txlens_summary = annotation_bundle$txlens_summary,
    minimal_component_length = minimal_component_length
  )
}

peak_exon_compute_transcript_stats <- function(peaks, transcript_ranges, transcript_lookup, sample_name) {
  subject_index <- peak_exon_first_hit_subject_index(peaks, transcript_ranges)

  if (!length(subject_index)) {
    return(data.table::data.table(
      tx_name = character(),
      tx_len = numeric(),
      nexon = numeric(),
      Sample = character()
    ))
  }

  transcript_ids <- unique(as.character(transcript_ranges$tx_name[subject_index]))
  transcript_stats <- transcript_lookup[tx_name %in% transcript_ids, .(tx_name, tx_len, nexon)]
  transcript_stats <- transcript_stats[!is.na(tx_name)]
  transcript_stats[, Sample := sample_name]
  transcript_stats[order(Sample, tx_name)]
}

peak_exon_compute_exon_size_data <- function(peaks, exon_ranges, sample_name) {
  subject_index <- peak_exon_first_hit_subject_index(peaks, exon_ranges)

  if (!length(subject_index)) {
    return(data.table::data.table(
      exonID = character(),
      tx_name = character(),
      rank = numeric(),
      length = numeric(),
      exonLoc = character(),
      Sample = character()
    ))
  }

  exon_info <- data.table::as.data.table(as.data.frame(S4Vectors::mcols(exon_ranges[subject_index])))
  output_columns <- intersect(
    c("exonID", "exon_id", "exon_name", "exon_rank", "tx_name", "length", "rank", "exonLoc"),
    names(exon_info)
  )
  exon_info <- unique(exon_info[, ..output_columns], by = "exonID")
  exon_info[, Sample := sample_name]
  exon_info[, exonLoc := factor(exonLoc, levels = peak_exon_location_levels())]
  exon_info[order(Sample, exonLoc, rank)]
}

peak_exon_collect_analysis_data <- function(
  upload_context,
  module_label,
  progress = NULL,
  include_exon_size = FALSE,
  include_transcript_stats = FALSE
) {
  saved_files <- peak_exon_require_saved_upload(upload_context, module_label)
  annotation_bundle <- peak_exon_load_annotation_bundle(upload_context$species)
  txdb_features <- peak_exon_load_txdb_features(upload_context$species, minimal_component_length = 100)
  feature_bundle <- peak_exon_feature_bundle(
    species = upload_context$species,
    txdb_features = txdb_features,
    annotation_bundle = annotation_bundle
  )
  sample_metadata <- rnameta_build_sample_metadata(saved_files)
  sample_names <- rnameta_sample_display_names(sample_metadata)
  interval_count <- 0L

  if (isTRUE(include_exon_size)) {
    exon_size_data_list <- vector("list", length(saved_files))
  } else {
    exon_size_data_list <- NULL
  }

  if (isTRUE(include_transcript_stats)) {
    transcript_stats_list <- vector("list", length(saved_files))
  } else {
    transcript_stats_list <- NULL
  }

  peak_exon_progress_set(progress, module_label, 16, "Loading saved BED files and species exon annotation bundle")
  peak_exon_progress_set(progress, module_label, 32, "Preparing transcript and exon feature structures from the legacy TxDb path")

  for (index in seq_along(saved_files)) {
    file_info <- saved_files[[index]]
    sample_name <- sample_names[[index]]
    progress_value <- 42 + floor(index / length(saved_files) * 44)

    peak_exon_progress_set(
      progress,
      module_label,
      progress_value,
      sprintf("Annotating exon-aware transcript hits for %s (%d/%d)", basename(file_info$file_name), index, length(saved_files))
    )

    raw_peaks <- peak_distribution_read_peaks(file_info$file_path, annotation_bundle$txdb, sample_name = sample_name)
    analysis_peaks <- GenomicRanges::resize(raw_peaks, width = 1, fix = "center")

    if (isTRUE(include_exon_size)) {
      exon_size_data_list[[index]] <- peak_exon_compute_exon_size_data(
        peaks = analysis_peaks,
        exon_ranges = feature_bundle$exon_ranges,
        sample_name = sample_name
      )
    }

    if (isTRUE(include_transcript_stats)) {
      transcript_stats_list[[index]] <- peak_exon_compute_transcript_stats(
        peaks = analysis_peaks,
        transcript_ranges = feature_bundle$transcript_ranges,
        transcript_lookup = feature_bundle$transcript_lookup,
        sample_name = sample_name
      )
    }

    interval_count <- interval_count + length(raw_peaks)
  }

  list(
    species = upload_context$species,
    sample_metadata = sample_metadata,
    sample_names = sample_names,
    sample_display_names = sample_names,
    sample_original_names = rnameta_sample_original_names(sample_metadata),
    sample_count = length(saved_files),
    interval_count = interval_count,
    exon_size_data = if (isTRUE(include_exon_size)) data.table::rbindlist(exon_size_data_list, use.names = TRUE, fill = TRUE) else NULL,
    transcript_stats = if (isTRUE(include_transcript_stats)) data.table::rbindlist(transcript_stats_list, use.names = TRUE, fill = TRUE) else NULL
  )
}
