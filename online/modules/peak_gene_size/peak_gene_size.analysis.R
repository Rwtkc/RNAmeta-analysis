peak_gene_size_compute_lengths <- function(peaks, transcript_ranges, txlens, sample_name) {
  transcript_ranges_unstranded <- transcript_ranges
  GenomicRanges::strand(transcript_ranges_unstranded) <- "*"
  overlaps <- GenomicRanges::findOverlaps(peaks, transcript_ranges_unstranded)

  if (!length(overlaps)) {
    return(data.table::data.table(
      tx_name = character(),
      tx_len = numeric(),
      Sample = character()
    ))
  }

  hit_index <- which(!duplicated(S4Vectors::queryHits(overlaps)))
  subject_index <- S4Vectors::subjectHits(overlaps[hit_index])
  transcript_hits <- unique(as.character(transcript_ranges$tx_name[subject_index]))

  txlens_dt <- data.table::as.data.table(txlens)
  lengths_dt <- txlens_dt[tx_name %in% transcript_hits, .(tx_name, tx_len)]
  lengths_dt <- unique(lengths_dt, by = "tx_name")
  lengths_dt <- lengths_dt[!is.na(tx_len) & tx_len > 0]
  lengths_dt[, Sample := sample_name]
  lengths_dt[order(Sample, tx_len)]
}

run_peak_gene_size_analysis <- function(upload_context, controls = peak_gene_size_default_controls(), progress = NULL) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (is.null(upload_context) || !length(saved_files)) {
    stop("Save a BED file in Upload / Run before running Peak Gene Size.")
  }

  if (is.null(upload_context$species) || !nzchar(upload_context$species)) {
    stop("A species must be saved in Upload / Run before running Peak Gene Size.")
  }

  controls <- utils::modifyList(peak_gene_size_default_controls(), controls)
  peak_gene_size_progress_update(progress, 14, "Loading saved BED files and species transcript annotation bundle")
  annotation_bundle <- peak_gene_size_load_annotation_bundle(upload_context$species)
  peak_gene_size_progress_update(progress, 36, "Preparing transcript ranges from the legacy TxDb path")
  transcript_ranges <- peak_gene_size_prepare_transcript_ranges(annotation_bundle$txdb)
  peak_gene_size_data_list <- vector("list", length(saved_files))
  sample_metadata <- rnameta_build_sample_metadata(saved_files)
  interval_count <- 0L

  for (index in seq_along(saved_files)) {
    file_info <- saved_files[[index]]
    sample_name <- sample_metadata[[index]]$display_name
    peak_gene_size_progress_update(
      progress,
      48 + floor(index / length(saved_files) * 34),
      sprintf("Collecting transcript lengths for %s (%d/%d)", basename(file_info$file_name), index, length(saved_files))
    )
    analysis_peaks <- peak_gene_size_read_peaks(file_info$file_path, annotation_bundle$txdb, sample_name = sample_name)
    peak_gene_size_data_list[[index]] <- peak_gene_size_compute_lengths(analysis_peaks, transcript_ranges, annotation_bundle$txlens, sample_name)
    interval_count <- interval_count + length(analysis_peaks)
  }
  peak_gene_size_data <- data.table::rbindlist(peak_gene_size_data_list, use.names = TRUE, fill = TRUE)
  peak_gene_size_progress_update(progress, 100, "Building chart payload and finalizing results")

  list(
    chart_payload = peak_gene_size_build_chart_payload(peak_gene_size_data, sample_metadata = sample_metadata),
    peak_gene_size_data = peak_gene_size_data,
    species = upload_context$species,
    sample_name = rnameta_sample_original_names(sample_metadata),
    sample_display_names = rnameta_sample_display_names(sample_metadata),
    sample_count = length(saved_files),
    interval_count = interval_count,
    transcript_count = nrow(peak_gene_size_data)
  )
}
