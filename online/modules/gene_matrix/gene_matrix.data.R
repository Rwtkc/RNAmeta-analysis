gene_matrix_load_annotation_bundle <- function(species) {
  if (exists("peak_distribution_cached_annotation_bundle", mode = "function")) {
    return(peak_distribution_cached_annotation_bundle(species))
  }

  bundle <- peak_distribution_load_annotation_bundle(species)
  bundle$txlens_summary <- data.table::as.data.table(bundle$txlens)[, lapply(.SD, max), by = gene_id]
  bundle
}

gene_matrix_txdb_features <- function(species, minimal_component_length = 100) {
  if (exists("peak_distribution_cached_txdb_features", mode = "function")) {
    return(peak_distribution_cached_txdb_features(species, minimal_component_length = minimal_component_length))
  }

  bundle <- gene_matrix_load_annotation_bundle(species)
  peak_distribution_build_txdb_features(
    txdb = bundle$txdb,
    txlens_summary = bundle$txlens_summary,
    minimal_component_length = minimal_component_length
  )
}

gene_matrix_build_transcript_meta <- function(gff_data) {
  transcript_meta <- data.table::as.data.table(as.data.frame(gff_data))
  transcript_meta <- transcript_meta[type == "transcript", .(gene_id, transcript_id, gene_type, gene_name)]
  transcript_meta[!is.na(transcript_id)]
}

gene_matrix_transcript_meta <- function(species, annotation_bundle = NULL) {
  builder <- function() {
    bundle <- annotation_bundle
    if (is.null(bundle)) {
      bundle <- gene_matrix_load_annotation_bundle(species)
    }
    gene_matrix_build_transcript_meta(bundle$gff)
  }

  if (exists("gene_matrix_cached_transcript_meta", mode = "function")) {
    return(gene_matrix_cached_transcript_meta(species, builder))
  }

  builder()
}

gene_matrix_build_gene_name_lookup <- function(transcript_meta) {
  transcript_table <- data.table::as.data.table(transcript_meta)
  transcript_table <- transcript_table[!is.na(transcript_id) & nzchar(transcript_id)]
  transcript_table <- transcript_table[!is.na(gene_name) & nzchar(gene_name)]

  if (!nrow(transcript_table)) {
    return(stats::setNames(character(), character()))
  }

  transcript_table <- transcript_table[, .(gene_name = gene_name[[1L]]), by = transcript_id]
  stats::setNames(transcript_table$gene_name, transcript_table$transcript_id)
}

gene_matrix_feature_tx_names <- function(genomic_region) {
  if (inherits(genomic_region, "GRangesList") || inherits(genomic_region, "CompressedGRangesList")) {
    feature_ranges <- BiocGenerics::unlist(genomic_region)
    feature_ranges$tx_name <- names(feature_ranges)
  } else {
    feature_ranges <- genomic_region

    if (!"tx_name" %in% names(S4Vectors::mcols(feature_ranges))) {
      if (!is.null(names(feature_ranges)) && any(nzchar(names(feature_ranges)))) {
        feature_ranges$tx_name <- names(feature_ranges)
      } else if ("tx_id" %in% names(S4Vectors::mcols(feature_ranges))) {
        feature_ranges$tx_name <- as.character(feature_ranges$tx_id)
      } else {
        feature_ranges$tx_name <- NA_character_
      }
    }
  }

  feature_ranges <- feature_ranges[!is.na(feature_ranges$tx_name) & nzchar(feature_ranges$tx_name)]

  if (length(feature_ranges)) {
    feature_ranges$tx_name <- as.character(feature_ranges$tx_name)
  }

  feature_ranges
}

gene_matrix_build_prepared_features <- function(txdb_features) {
  feature_sets <- list(
    promoters = txdb_features$promoters,
    fiveUTRs = txdb_features$fiveUTRs,
    cds = txdb_features$cds,
    threeUTRs = txdb_features$threeUTRs,
    introns = txdb_features$introns,
    stopCodon = txdb_features$stopCodon,
    transcripts = txdb_features$transcripts
  )

  lapply(feature_sets, gene_matrix_feature_tx_names)
}

gene_matrix_prepared_features <- function(species, txdb_features = NULL, minimal_component_length = 100) {
  builder <- function() {
    features <- txdb_features
    if (is.null(features)) {
      features <- gene_matrix_txdb_features(species, minimal_component_length = minimal_component_length)
    }
    gene_matrix_build_prepared_features(features)
  }

  if (exists("gene_matrix_cached_prepared_features", mode = "function")) {
    return(gene_matrix_cached_prepared_features(species, minimal_component_length, builder))
  }

  builder()
}

gene_matrix_feature_transcript_ids <- function(peaks, feature_ranges) {
  if (!length(feature_ranges)) {
    return(character())
  }

  overlaps <- GenomicRanges::findOverlaps(peaks, feature_ranges, ignore.strand = TRUE)

  if (!length(overlaps)) {
    return(character())
  }

  overlaps <- overlaps[peak_distribution_first_hit_index(S4Vectors::queryHits(overlaps))]
  sort(unique(as.character(feature_ranges$tx_name[S4Vectors::subjectHits(overlaps)])))
}

gene_matrix_collect_transcript_ids <- function(peaks, prepared_features) {
  transcript_ids <- unique(unlist(
    lapply(prepared_features, function(feature_set) gene_matrix_feature_transcript_ids(peaks, feature_set)),
    use.names = FALSE
  ))

  sort(transcript_ids[nzchar(transcript_ids)])
}

gene_matrix_compute_gene_set_fast <- function(peaks, prepared_features, gene_name_lookup = NULL, transcript_meta = NULL) {
  transcript_ids <- gene_matrix_collect_transcript_ids(peaks, prepared_features)

  if (!length(transcript_ids)) {
    return(character())
  }

  if (is.null(gene_name_lookup)) {
    gene_name_lookup <- gene_matrix_build_gene_name_lookup(transcript_meta)
  }

  gene_names <- unname(gene_name_lookup[transcript_ids])
  gene_names <- stats::na.omit(gene_names)
  sort(unique(gene_names[nzchar(gene_names)]))
}

gene_matrix_extract_annotation_table <- function(annotation_vector) {
  if (is.null(annotation_vector) || !length(annotation_vector)) {
    return(data.table::data.table(Index = integer(), feature = character(), transcript_id = character()))
  }

  annotation_table <- data.table::data.table(
    Index = seq_along(annotation_vector),
    annotation = as.character(annotation_vector)
  )
  annotation_table <- annotation_table[!is.na(annotation) & nzchar(annotation)]

  if (!nrow(annotation_table)) {
    return(data.table::data.table(Index = integer(), feature = character(), transcript_id = character()))
  }

  annotation_table[, feature := sub(" \\(.+$", "", annotation)]
  annotation_table[, transcript_id := sub("^.*\\((.+)\\)$", "\\1", annotation)]
  annotation_table[, transcript_id := sub(",.*$", "", transcript_id)]
  annotation_table[!is.na(feature) & !is.na(transcript_id), .(Index, feature, transcript_id)]
}

gene_matrix_peak_table <- function(peaks) {
  data.table::data.table(
    seqnames = as.character(GenomeInfoDb::seqnames(peaks)),
    start = BiocGenerics::start(peaks),
    end = BiocGenerics::end(peaks),
    peakname = as.character(peaks$peakname)
  )[, `:=`(
    Index = .I,
    peak = paste0(seqnames, ":", start, "-", end, "|", peakname)
  )]
}

gene_matrix_compute_feature_table <- function(annotation_result, peaks, gff_data) {
  main_annotations <- gene_matrix_extract_annotation_table(annotation_result$annotation)
  stop_annotations <- gene_matrix_extract_annotation_table(annotation_result$annotation_stopcodon)
  transcript_annotations <- gene_matrix_extract_annotation_table(annotation_result$annotation_transcript)
  peak_table <- gene_matrix_peak_table(peaks)

  peak_detail <- main_annotations[peak_table, on = "Index", nomatch = 0]
  peak_detail <- peak_detail[!is.na(feature), .(peakDetail = paste(peak, collapse = ",")), by = transcript_id]

  feature_counts <- data.table::rbindlist(
    list(
      main_annotations[, .(transcript_id, feature)],
      stop_annotations[, .(transcript_id, feature)],
      transcript_annotations[, .(transcript_id, feature)]
    ),
    use.names = TRUE,
    fill = TRUE
  )

  if (!nrow(feature_counts)) {
    return(data.table::data.table(
      gene_id = character(),
      transcript_id = character(),
      gene_name = character(),
      gene_type = character(),
      Promoter = integer(),
      UTR5 = integer(),
      CDS = integer(),
      UTR3 = integer(),
      `Stop Codon` = integer(),
      Intron = integer(),
      Transcript = integer(),
      peakDetail = character()
    ))
  }

  feature_counts <- feature_counts[, .N, by = .(transcript_id, feature)]
  feature_table <- data.table::dcast(feature_counts, transcript_id ~ feature, value.var = "N")
  feature_table[is.na(feature_table)] <- 0

  transcript_parts <- intersect(c("UTR5", "CDS", "UTR3", "Intron"), names(feature_table))
  feature_table[, Transcript := if (length(transcript_parts)) rowSums(.SD) else 0, .SDcols = transcript_parts]

  transcript_meta <- data.table::as.data.table(as.data.frame(gff_data))
  transcript_meta <- transcript_meta[type == "transcript", .(gene_id, transcript_id, gene_type, gene_name)]

  results <- merge(feature_table, transcript_meta, by = "transcript_id", all.x = TRUE)
  results <- merge(results, peak_detail, by = "transcript_id", all.x = TRUE)

  output_order <- c("gene_id", "transcript_id", "gene_name", "gene_type", "Promoter", "UTR5", "CDS", "UTR3", "Stop Codon", "Intron", "Transcript", "peakDetail")

  for (column_name in output_order) {
    if (!column_name %in% names(results)) {
      results[, (column_name) := if (column_name %in% c("gene_id", "transcript_id", "gene_name", "gene_type", "peakDetail")) NA_character_ else 0]
    }
  }

  data.table::setcolorder(results, output_order)
  results
}

gene_matrix_gene_set <- function(feature_table) {
  sort(unique(stats::na.omit(feature_table$gene_name)))
}
