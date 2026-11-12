site_profile_density_points <- function(hit_table, x_min, x_max, feature_order = NULL, point_count = 512) {
  if (is.null(hit_table) || !nrow(hit_table)) {
    return(data.table::data.table(
      Panel = character(),
      Sample = character(),
      Position = numeric(),
      Density = numeric()
    ))
  }

  density_rows <- lapply(split(hit_table, by = c("feature", "Sample"), keep.by = TRUE), function(sample_hits) {
    density_values <- sample_hits$rel_pos[is.finite(sample_hits$rel_pos)]
    if (!length(density_values)) {
      density_values <- c(0, 0)
    } else if (length(density_values) < 2) {
      density_values <- rep(density_values[[1]], 2)
    }

    density_curve <- stats::density(
      density_values,
      from = x_min,
      to = x_max,
      n = point_count,
      bw = if (length(unique(density_values)) < 2) 1 else "nrd0",
      na.rm = TRUE
    )

    data.table::data.table(
      Panel = unique(as.character(sample_hits$feature))[[1]],
      Sample = unique(as.character(sample_hits$Sample))[[1]],
      Position = as.numeric(density_curve$x),
      Density = as.numeric(density_curve$y)
    )
  })

  density_data <- data.table::rbindlist(density_rows, use.names = TRUE, fill = TRUE)

  if (!is.null(feature_order) && length(feature_order)) {
    density_data[, Panel := factor(Panel, levels = feature_order)]
    data.table::setorder(density_data, Panel, Sample, Position)
    density_data[, Panel := as.character(Panel)]
  }

  density_data
}

site_profile_load_annotation_bundle <- function(species) {
  if (exists("peak_distribution_cached_annotation_bundle", mode = "function")) {
    return(peak_distribution_cached_annotation_bundle(species))
  }

  peak_distribution_load_annotation_bundle(species)
}

site_profile_load_txdb_features <- function(species, minimal_component_length = 100) {
  if (exists("peak_distribution_cached_txdb_features", mode = "function")) {
    return(peak_distribution_cached_txdb_features(species, minimal_component_length = minimal_component_length))
  }

  annotation_bundle <- site_profile_load_annotation_bundle(species)
  peak_distribution_build_txdb_features(
    txdb = annotation_bundle$txdb,
    txlens_summary = annotation_bundle$txlens_summary,
    minimal_component_length = minimal_component_length
  )
}

site_profile_center_peaks <- function(peaks) {
  GenomicRanges::resize(peaks, width = 1, fix = "center")
}

site_profile_heatmap_peaks <- function(peaks) {
  GenomicRanges::resize(peaks, width = 15, fix = "center")
}

site_profile_query_table <- function(query_ranges, sample_name) {
  peak_table <- data.table::as.data.table(query_ranges)
  peak_table[, index := .I]
  peak_table[, Sample := sample_name]
  peak_table
}

site_profile_empty_hits <- function(sample_name, feature_label, x_min, x_max, flank_size) {
  data.table::data.table(
    seqnames = character(),
    start = integer(),
    end = integer(),
    width = integer(),
    strand = character(),
    Sample = sample_name,
    feature = feature_label,
    rel_pos = numeric(),
    flankSize = flank_size,
    xmin = x_min,
    xmax = x_max
  )
}

site_profile_collect_boundary_hits <- function(query_ranges, windows, sample_name, feature_label, rel_offset, x_min, x_max, flank_size) {
  hits <- GenomicFeatures::mapToTranscripts(query_ranges, windows, ignore.strand = FALSE)

  if (!length(hits)) {
    return(site_profile_empty_hits(sample_name, feature_label, x_min, x_max, flank_size))
  }

  peak_table <- site_profile_query_table(query_ranges, sample_name)
  hits_table <- data.table::as.data.table(hits)
  data.table::setnames(hits_table, c("seqnames", "xHits"), c("transcriptID", "index"))

  peak_hits <- peak_table[hits_table, on = "index"]
  peak_hits[, c("index", "transcriptsHits") := NULL]
  peak_hits[, feature := feature_label]
  peak_hits[, rel_pos := i.start - rel_offset]
  peak_hits[, flankSize := flank_size]
  peak_hits[, xmin := x_min]
  peak_hits[, xmax := x_max]
  peak_hits
}

site_profile_map_to_transcript_space <- function(source_file, species, peaks, exons, cache_tag = "default") {
  site_profile_cached_transcript_space_peaks(source_file, species, cache_tag = cache_tag, function() {
    transcript_peaks <- GenomicFeatures::mapToTranscripts(peaks, exons)
    if (length(transcript_peaks)) {
      GenomicRanges::strand(transcript_peaks) <- "*"
    }
    transcript_peaks
  })
}

site_profile_evenly_sample_windows <- function(windows, sample_n = 10000) {
  if (is.null(windows) || !length(windows) || length(windows) <= sample_n) {
    return(windows)
  }

  sampled_index <- unique(round(seq(1, length(windows), length.out = sample_n)))
  windows[sampled_index]
}

site_profile_score_matrix_rna <- function(target, windows) {
  widths <- unique(IRanges::width(windows))
  if (length(widths) > 1L) stop("Heatmap windows must have equal widths.")
  ncols <- if (length(widths)) widths[[1]] else 0L
  if (!length(windows)) return(matrix(numeric(), 0L, ncols))
  chromosomes <- as.character(GenomeInfoDb::seqnames(windows))
  limits <- GenomeInfoDb::seqlengths(windows)[chromosomes]
  valid <- IRanges::start(windows) >= 1L & (is.na(limits) | IRanges::end(windows) <= limits)
  windows <- windows[valid]
  chromosomes <- chromosomes[valid]
  result <- matrix(0, nrow = length(windows), ncol = ncols)
  if (!length(windows) || !length(target)) return(result)
  # Pad coverage to the full annotation window; never shrink its coordinate axis.
  coverage_widths <- tapply(IRanges::end(windows), chromosomes, max)
  target_levels <- GenomeInfoDb::seqlevels(target)
  coverage_widths <- coverage_widths[intersect(names(coverage_widths), target_levels)]
  target <- GenomeInfoDb::keepSeqlevels(target, names(coverage_widths), pruning.mode = "coarse")
  target_rle <- GenomicRanges::coverage(target)
  for (chromosome in names(coverage_widths)) {
    indices <- which(chromosomes == chromosome)
    signal <- target_rle[[chromosome]]
    required <- coverage_widths[[chromosome]]
    if (length(signal) < required) signal <- c(signal, S4Vectors::Rle(0, required - length(signal)))
    views <- IRanges::Views(signal, start = IRanges::start(windows)[indices], end = IRanges::end(windows)[indices])
    result[indices, ] <- t(IRanges::viewApply(views, as.vector))
  }
  minus <- as.character(GenomicRanges::strand(windows)) == "-"
  if (any(minus)) result[minus, ] <- result[minus, ncols:1L, drop = FALSE]
  result
}

site_profile_prepare_heatmap_matrix <- function(tag_matrix) {
  if (is.null(tag_matrix)) return(matrix(numeric(), nrow = 0, ncol = 0))
  if (!nrow(tag_matrix)) return(tag_matrix)

  row_order <- order(rowSums(tag_matrix))
  tag_matrix <- tag_matrix[row_order, , drop = FALSE]
  non_zero_rows <- rowSums(tag_matrix) > 0

  tag_matrix <- tag_matrix[non_zero_rows, , drop = FALSE]

  tag_matrix
}

site_profile_runtime_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$boundary_windows <- new.env(hash = TRUE, parent = emptyenv())
  env$transcript_space_peaks <- new.env(hash = TRUE, parent = emptyenv())
  env$heatmap_panels <- new.env(hash = TRUE, parent = emptyenv())
  env
})

site_profile_runtime_cache_key <- function(...) {
  paste(vapply(list(...), as.character, character(1)), collapse = "::")
}

site_profile_file_fingerprint <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return("missing")
  }

  info <- file.info(path)
  paste(
    normalizePath(path, winslash = "/"),
    unname(info$size[[1]]),
    as.numeric(unname(info$mtime[[1]])),
    sep = "|"
  )
}

site_profile_cached_boundary_windows <- function(kind, species, flank_size, builder) {
  key <- site_profile_runtime_cache_key(kind, species, flank_size)

  if (exists(key, envir = site_profile_runtime_cache_env$boundary_windows, inherits = FALSE)) {
    return(get(key, envir = site_profile_runtime_cache_env$boundary_windows, inherits = FALSE))
  }

  windows <- builder()
  assign(key, windows, envir = site_profile_runtime_cache_env$boundary_windows)
  windows
}

site_profile_cached_transcript_space_peaks <- function(source_file, species, cache_tag = "default", builder) {
  key <- site_profile_runtime_cache_key(cache_tag, species, site_profile_file_fingerprint(source_file))

  if (exists(key, envir = site_profile_runtime_cache_env$transcript_space_peaks, inherits = FALSE)) {
    return(get(key, envir = site_profile_runtime_cache_env$transcript_space_peaks, inherits = FALSE))
  }

  mapped_peaks <- builder()
  assign(key, mapped_peaks, envir = site_profile_runtime_cache_env$transcript_space_peaks)
  mapped_peaks
}

site_profile_cached_heatmap_panels <- function(kind, source_file, species, flank_size, builder, cache_token = NULL) {
  cache_token <- if (is.null(cache_token) || !nzchar(cache_token)) "default" else cache_token
  key <- site_profile_runtime_cache_key(kind, species, flank_size, cache_token, site_profile_file_fingerprint(source_file))

  if (exists(key, envir = site_profile_runtime_cache_env$heatmap_panels, inherits = FALSE)) {
    return(get(key, envir = site_profile_runtime_cache_env$heatmap_panels, inherits = FALSE))
  }

  panels <- builder()
  assign(key, panels, envir = site_profile_runtime_cache_env$heatmap_panels)
  panels
}

site_profile_runtime_cache_status <- function() {
  list(
    boundary_window_keys = ls(envir = site_profile_runtime_cache_env$boundary_windows),
    transcript_space_peak_keys = ls(envir = site_profile_runtime_cache_env$transcript_space_peaks),
    heatmap_panel_keys = ls(envir = site_profile_runtime_cache_env$heatmap_panels)
  )
}

site_profile_runtime_cache_reset <- function() {
  rm(
    list = ls(envir = site_profile_runtime_cache_env$boundary_windows),
    envir = site_profile_runtime_cache_env$boundary_windows
  )
  rm(
    list = ls(envir = site_profile_runtime_cache_env$transcript_space_peaks),
    envir = site_profile_runtime_cache_env$transcript_space_peaks
  )
  rm(
    list = ls(envir = site_profile_runtime_cache_env$heatmap_panels),
    envir = site_profile_runtime_cache_env$heatmap_panels
  )
  invisible(TRUE)
}
