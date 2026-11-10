run_gene_matrix_analysis <- function(upload_context, progress = NULL) {
  saved_files <- gene_matrix_saved_files(upload_context)

  if (!length(saved_files)) {
    stop("Save BED files in Upload / Run before running Gene Matrix.")
  }

  if (!gene_matrix_sample_count_ok(length(saved_files))) {
    stop(gene_matrix_required_sample_message())
  }

  if (is.null(upload_context$species) || !nzchar(upload_context$species)) {
    stop("A species must be saved in Upload / Run before running Gene Matrix.")
  }

  gene_matrix_progress_update(progress, 10, "Loading species annotation bundle for multi-sample overlap analysis")
  annotation_bundle <- gene_matrix_load_annotation_bundle(upload_context$species)
  gene_matrix_progress_update(progress, 28, "Constructing transcript-level feature annotation structure")
  txdb_features <- gene_matrix_txdb_features(upload_context$species, minimal_component_length = 100)
  prepared_features <- gene_matrix_prepared_features(
    species = upload_context$species,
    txdb_features = txdb_features,
    minimal_component_length = 100
  )
  transcript_meta <- gene_matrix_transcript_meta(upload_context$species, annotation_bundle = annotation_bundle)
  gene_name_lookup <- gene_matrix_build_gene_name_lookup(transcript_meta)

  sample_metadata <- rnameta_build_sample_metadata(saved_files)
  sample_labels <- rnameta_sample_display_names(sample_metadata)
  sample_results <- vector("list", length(saved_files))

  for (index in seq_along(saved_files)) {
    file_info <- saved_files[[index]]
    progress_value <- 28 + floor(index / length(saved_files) * 54)

    gene_matrix_progress_update(
      progress,
      progress_value,
      sprintf("Collecting overlap gene set for %s (%d/%d)", basename(file_info$file_name), index, length(saved_files))
    )

    raw_peaks <- peak_distribution_read_peaks(file_info$file_path, annotation_bundle$txdb, sample_name = sample_labels[[index]])
    analysis_peaks <- GenomicRanges::resize(raw_peaks, width = 1, fix = "center")
    gene_set <- gene_matrix_compute_gene_set_fast(
      peaks = analysis_peaks,
      prepared_features = prepared_features,
      gene_name_lookup = gene_name_lookup,
      transcript_meta = transcript_meta
    )

    sample_results[[index]] <- list(
      sample_name = sample_labels[[index]],
      interval_count = length(raw_peaks),
      gene_set = gene_set
    )
  }

  gene_sets <- setNames(lapply(sample_results, `[[`, "gene_set"), sample_labels)
  gene_matrix_progress_update(progress, 92, "Building venn diagram payload and finalizing shared-gene summary")
  union_genes <- sort(unique(unlist(gene_sets, use.names = FALSE)))
  chart_payload <- gene_matrix_build_chart_payload(gene_sets, sample_labels, sample_metadata = sample_metadata)

  list(
    chart_payload = chart_payload,
    species = upload_context$species,
    sample_count = length(saved_files),
    sample_labels = sample_labels,
    sample_name = rnameta_sample_original_names(sample_metadata),
    sample_display_names = sample_labels,
    union_gene_count = length(union_genes),
    sample_results = sample_results,
    gene_sets = gene_sets
  )
}
