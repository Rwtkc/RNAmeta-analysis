peak_exon_runtime_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$feature_bundle <- new.env(hash = TRUE, parent = emptyenv())
  env
})

peak_exon_runtime_cache_key <- function(...) {
  paste(vapply(list(...), as.character, character(1)), collapse = "::")
}

peak_exon_cached_feature_bundle <- function(species, builder) {
  key <- peak_exon_runtime_cache_key(species)

  if (exists(key, envir = peak_exon_runtime_cache_env$feature_bundle, inherits = FALSE)) {
    return(get(key, envir = peak_exon_runtime_cache_env$feature_bundle, inherits = FALSE))
  }

  feature_bundle <- builder()
  assign(key, feature_bundle, envir = peak_exon_runtime_cache_env$feature_bundle)
  feature_bundle
}

peak_exon_prepare_transcript_lookup <- function(txlens, exons = NULL) {
  txlens_dt <- data.table::as.data.table(txlens)

  if (!"tx_name" %in% names(txlens_dt)) {
    row_names <- rownames(as.data.frame(txlens))
    if (length(row_names) == nrow(txlens_dt) && any(nzchar(row_names))) {
      txlens_dt[, tx_name := row_names]
    }
  }

  if (!"nexon" %in% names(txlens_dt) && !is.null(exons)) {
    exon_counts <- data.table::data.table(
      tx_name = names(exons),
      nexon = S4Vectors::elementNROWS(exons)
    )
    txlens_dt <- merge(txlens_dt, exon_counts, by = "tx_name", all.x = TRUE)
  }

  txlens_dt[, tx_len := suppressWarnings(as.numeric(tx_len))]
  txlens_dt[, nexon := suppressWarnings(as.numeric(nexon))]
  unique(txlens_dt[, .(tx_name, tx_len, nexon)], by = "tx_name")
}

peak_exon_prepare_exon_ranges <- function(exons) {
  exon_ranges <- BiocGenerics::unlist(exons)
  exon_ranges$tx_name <- names(exon_ranges)
  exon_ranges$length <- IRanges::width(exon_ranges)

  exon_counts <- data.table::data.table(name = exon_ranges$tx_name)[, .N, by = name]
  exon_count_values <- exon_counts$N
  names(exon_count_values) <- exon_counts$name

  # exonsBy(..., by = "tx") supplies exon_rank in transcript order on both
  # strands. Reversing negative-strand ranks would swap First and Last.
  exon_ranges$rank <- as.integer(exon_ranges$exon_rank)

  exon_total <- unname(exon_count_values[as.character(exon_ranges$tx_name)])
  exon_loc <- rep("Middle", length(exon_ranges))
  exon_loc[exon_ranges$rank == 1] <- "First"
  exon_loc[!is.na(exon_total) & exon_ranges$rank == exon_total] <- "Last"

  exon_ranges$exonLoc <- factor(exon_loc, levels = peak_exon_location_levels())
  exon_ranges$exonID <- paste0(exon_ranges$tx_name, ":", exon_ranges$rank, recycle0 = TRUE)
  exon_ranges
}

peak_exon_prepare_transcript_ranges <- function(transcripts) {
  transcript_ranges <- transcripts

  if (!"tx_name" %in% names(S4Vectors::mcols(transcript_ranges))) {
    if ("tx_name" %in% names(S4Vectors::mcols(transcript_ranges))) {
      transcript_ranges$tx_name <- as.character(transcript_ranges$tx_name)
    } else if (!is.null(names(transcript_ranges)) && any(nzchar(names(transcript_ranges)))) {
      transcript_ranges$tx_name <- names(transcript_ranges)
    } else if ("tx_id" %in% names(S4Vectors::mcols(transcript_ranges))) {
      transcript_ranges$tx_name <- as.character(transcript_ranges$tx_id)
    } else {
      transcript_ranges$tx_name <- NA_character_
    }
  }

  transcript_ranges[!is.na(transcript_ranges$tx_name) & nzchar(transcript_ranges$tx_name)]
}

peak_exon_feature_bundle <- function(species, txdb_features, annotation_bundle) {
  builder <- function() {
    list(
      exon_ranges = peak_exon_prepare_exon_ranges(txdb_features$exons),
      transcript_ranges = peak_exon_prepare_transcript_ranges(txdb_features$transcripts),
      transcript_lookup = peak_exon_prepare_transcript_lookup(annotation_bundle$txlens, exons = txdb_features$exons)
    )
  }

  peak_exon_cached_feature_bundle(species, builder)
}

peak_exon_first_hit_subject_index <- function(peaks, feature_ranges) {
  overlaps <- GenomicRanges::findOverlaps(peaks, feature_ranges, ignore.strand = TRUE)

  if (!length(overlaps)) {
    return(integer())
  }

  subject_hits <- S4Vectors::subjectHits(overlaps[peak_distribution_first_hit_index(S4Vectors::queryHits(overlaps))])
  unique(subject_hits)
}

peak_exon_runtime_cache_status <- function() {
  list(
    feature_bundle_keys = ls(envir = peak_exon_runtime_cache_env$feature_bundle)
  )
}

peak_exon_runtime_cache_reset <- function() {
  rm(list = ls(envir = peak_exon_runtime_cache_env$feature_bundle), envir = peak_exon_runtime_cache_env$feature_bundle)
  invisible(TRUE)
}
