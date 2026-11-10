peak_distribution_build_txdb_features <- function(txdb, txlens_summary, minimal_component_length = 50) {
  transcripts <- GenomicFeatures::transcripts(txdb)
  exons <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
  introns <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
  promoters <- GenomicFeatures::promoters(txdb, upstream = 1000, downstream = 0)
  five_utrs <- GenomicFeatures::fiveUTRsByTranscript(x = txdb, use.names = TRUE)
  upstream <- GenomicRanges::flank(five_utrs, 500, start = TRUE, both = FALSE)
  three_utrs <- GenomicFeatures::threeUTRsByTranscript(x = txdb, use.names = TRUE)
  downstream <- GenomicRanges::flank(three_utrs, 500, start = FALSE, both = FALSE)
  cds <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = TRUE)
  cds_list <- range(cds)
  stop_codon <- GenomicRanges::flank(cds_list, 200, start = FALSE, both = TRUE)
  start_codon <- GenomicRanges::flank(cds_list, 200, start = TRUE, both = TRUE)

  flag_utr5 <- sum(IRanges::width(five_utrs)) > minimal_component_length
  flag_utr3 <- sum(IRanges::width(three_utrs)) > minimal_component_length
  flag_cds <- sum(IRanges::width(cds)) > minimal_component_length

  name_utr5 <- names(five_utrs)[flag_utr5]
  name_utr3 <- names(three_utrs)[flag_utr3]
  name_cds <- names(cds)[flag_cds]
  name_mrna <- Reduce(intersect, list(name_utr5, name_utr3, name_cds))
  name_filtered_mrna <- intersect(name_mrna, names(exons))
  length_filtered_mrna <- txlens_summary$tx_name

  list(
    transcripts = transcripts,
    exons = exons,
    promoters = promoters,
    fiveUTRs = five_utrs,
    upstream = upstream,
    introns = introns,
    cds = cds,
    threeUTRs = three_utrs,
    stopCodon = stop_codon,
    startCodon = start_codon,
    downstream = downstream,
    filtered_mrna = intersect(name_filtered_mrna, length_filtered_mrna)
  )
}

peak_distribution_first_hit_index <- function(index_vector) {
  which(!duplicated(index_vector))
}

peak_distribution_collect_annotation_hits <- function(peaks, genomic_region, type, same_strand = FALSE) {
  if (inherits(genomic_region, "GRangesList") || inherits(genomic_region, "CompressedGRangesList")) {
    genomic_hits <- BiocGenerics::unlist(genomic_region)
    genomic_hits$tx_name <- names(genomic_hits)
    genomic_hits$length <- IRanges::width(genomic_hits)
  } else {
    genomic_hits <- genomic_region
    genomic_hits$length <- IRanges::width(genomic_hits)
  }

  if (type %in% c("Intron", "CDS")) {
    unique_regions <- genomic_hits[!duplicated(genomic_hits$tx_name)]
    temp <- data.table::data.table(name = genomic_hits$tx_name)[, .N, by = name]
    region_lengths <- temp$N
    names(region_lengths) <- temp$name
    region_strand <- as.character(GenomicRanges::strand(unique_regions))

    genomic_hits$rank <- unlist(lapply(seq_along(region_strand), function(index) {
      rank <- seq_len(region_lengths[[index]])
      if (identical(region_strand[[index]], "-")) {
        rank <- rev(rank)
      }
      rank
    }))
  }

  overlaps <- if (same_strand) {
    GenomicRanges::findOverlaps(peaks, genomic_hits)
  } else {
    unstranded_hits <- genomic_hits
    GenomicRanges::strand(unstranded_hits) <- "*"
    GenomicRanges::findOverlaps(peaks, unstranded_hits)
  }

  if (length(overlaps) == 0) {
    return(NA)
  }

  overlaps <- overlaps[peak_distribution_first_hit_index(S4Vectors::queryHits(overlaps))]
  query_index <- S4Vectors::queryHits(overlaps)
  subject_index <- S4Vectors::subjectHits(overlaps)
  hits <- genomic_hits[subject_index]
  gene_id <- hits$tx_name

  if (type == "Intron") {
    genomic_region_lengths <- table(genomic_hits$tx_name)
    annotation <- paste(type, " (", gene_id, ", intron ", hits$rank, " of ", genomic_region_lengths[gene_id], ")", sep = "")
  } else if (type == "CDS") {
    genomic_region_lengths <- table(genomic_hits$tx_name)
    annotation <- paste(type, " (", gene_id, ", CDS ", hits$rank, " of ", genomic_region_lengths[gene_id], ", ", hits$length, "bp)", sep = "")
  } else if (type == "fiveUTR") {
    annotation <- paste0("UTR5 (", gene_id, ")")
  } else if (type == "threeUTR") {
    annotation <- paste0("UTR3 (", gene_id, ")")
  } else if (type == "Stopcodon") {
    annotation <- paste0("Stop Codon (", gene_id, ")")
  } else if (type == "Startcodon") {
    annotation <- paste0("Start Codon (", gene_id, ")")
  } else {
    annotation <- paste0(type, " (", gene_id, ")")
  }

  list(queryIndex = query_index, annotation = annotation, gene = gene_id)
}

peak_distribution_update_annotations <- function(peaks, genomic_region, type, annotation_state, same_strand = FALSE) {
  hits <- peak_distribution_collect_annotation_hits(peaks, genomic_region, type, same_strand = same_strand)
  if (length(hits) > 1) {
    hit_index <- hits$queryIndex
    annotation_state$annotation[hit_index] <- hits$annotation
    annotation_state$detailGenomicAnnotation[hit_index, type] <- TRUE
  }
  annotation_state
}

peak_distribution_write_detail_annotations <- function(annotation_result, source_file, output_dir, sample_index) {
  if (is.null(output_dir) || !nzchar(output_dir)) {
    return(NULL)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  delimiter <- peak_distribution_detect_delimiter(source_file)
  bed_table <- data.table::fread(source_file, header = FALSE, sep = delimiter)
  original_cols <- colnames(bed_table)
  merged_table <- data.table::copy(bed_table)

  annotation_dt <- annotation_result$annotation_table
  stop_codon_dt <- annotation_result$stop_codon_table
  start_codon_dt <- annotation_result$start_codon_table

  data.table::setnames(annotation_dt, c("start_end", "annotation"), c("start_end", "region_type"))
  annotation_dt[, c("start", "end") := data.table::tstrsplit(start_end, "-", fixed = TRUE)]
  annotation_dt[, c("start_end", "end") := NULL]
  annotation_dt[, start := as.integer(start)]
  merged_table <- annotation_dt[merged_table, on = .(start >= V2, start <= V3), nomatch = 0]
  data.table::setnames(merged_table, c("start", "start.1"), c("V2", "V3"))
  region_cols <- c(original_cols, "region_type")
  merged_table <- merged_table[, ..region_cols]

  data.table::setnames(stop_codon_dt, c("start_end", "annotation"), c("start_end", "stop_codon"))
  stop_codon_dt[, c("start", "end") := data.table::tstrsplit(start_end, "-", fixed = TRUE)]
  stop_codon_dt[, c("start_end", "end") := NULL]
  stop_codon_dt[, start := as.integer(start)]
  merged_table <- stop_codon_dt[merged_table, on = .(start >= V2, start <= V3), nomatch = 0]
  data.table::setnames(merged_table, c("start", "start.1"), c("V2", "V3"))
  stop_cols <- c(original_cols, "region_type", "stop_codon")
  merged_table <- merged_table[, ..stop_cols]

  data.table::setnames(start_codon_dt, c("start_end", "annotation"), c("start_end", "start_codon"))
  start_codon_dt[, c("start", "end") := data.table::tstrsplit(start_end, "-", fixed = TRUE)]
  start_codon_dt[, c("start_end", "end") := NULL]
  start_codon_dt[, start := as.integer(start)]
  merged_table <- start_codon_dt[merged_table, on = .(start >= V2, start <= V3), nomatch = 0]
  data.table::setnames(merged_table, c("start", "start.1"), c("V2", "V3"))
  start_cols <- c(original_cols, "region_type", "stop_codon", "start_codon")
  merged_table <- merged_table[, ..start_cols]

  output_file <- file.path(output_dir, sprintf("peak_details_%s.bed", sample_index))
  data.table::fwrite(merged_table, output_file, sep = "\t")
  output_file
}

peak_distribution_annotate_peaks <- function(peaks, txdb_features, source_file, detail_output_dir = NULL, sample_index = 1, genomic_annotation_priority = peak_distribution_annotation_priority()) {
  annotation <- rep(NA_character_, length(peaks))
  flag <- rep(FALSE, length(peaks))
  detail_annotation <- data.frame(
    genic = flag,
    Intergenic = flag,
    Promoter = flag,
    fiveUTR = flag,
    threeUTR = flag,
    CDS = flag,
    Intron = flag,
    Stopcodon = flag
  )

  annotation_state <- list(annotation = annotation, detailGenomicAnnotation = detail_annotation)
  stop_state <- list(annotation = annotation, detailGenomicAnnotation = detail_annotation)
  exon_state <- list(annotation = annotation, detailGenomicAnnotation = detail_annotation)
  transcript_state <- list(annotation = annotation, detailGenomicAnnotation = detail_annotation)
  start_state <- list(annotation = annotation, detailGenomicAnnotation = detail_annotation)

  exon_state <- peak_distribution_update_annotations(peaks, txdb_features$exons, "CDS", exon_state, same_strand = FALSE)

  for (annotation_priority in rev(genomic_annotation_priority)) {
    if (annotation_priority == "Intergenic") {
      annotation[is.na(annotation)] <- "Intergenic"
      annotation_state$annotation <- annotation
    } else if (annotation_priority == "Intron") {
      annotation_state <- peak_distribution_update_annotations(peaks, txdb_features$introns, "Intron", annotation_state, same_strand = FALSE)
    } else if (annotation_priority == "CDS") {
      annotation_state <- peak_distribution_update_annotations(peaks, txdb_features$cds, "CDS", annotation_state, same_strand = FALSE)
    } else if (annotation_priority == "3UTR") {
      annotation_state <- peak_distribution_update_annotations(peaks, txdb_features$threeUTRs, "threeUTR", annotation_state, same_strand = FALSE)
    } else if (annotation_priority == "5UTR") {
      annotation_state <- peak_distribution_update_annotations(peaks, txdb_features$fiveUTRs, "fiveUTR", annotation_state, same_strand = FALSE)
    } else if (annotation_priority == "Stopcodon") {
      stop_state <- peak_distribution_update_annotations(peaks, txdb_features$stopCodon, "Stopcodon", stop_state, same_strand = FALSE)
    } else if (annotation_priority == "Promoter") {
      annotation_state <- peak_distribution_update_annotations(peaks, txdb_features$promoters, "Promoter", annotation_state, same_strand = FALSE)
    } else if (annotation_priority == "Startcodon") {
      start_state <- peak_distribution_update_annotations(peaks, txdb_features$startCodon, "Startcodon", start_state, same_strand = FALSE)
    } else if (annotation_priority == "Transcript") {
      transcript_state <- peak_distribution_update_annotations(peaks, txdb_features$transcripts, "Transcript", transcript_state, same_strand = FALSE)
    }

    annotation <- annotation_state$annotation
    detail_annotation <- annotation_state$detailGenomicAnnotation
  }

  genic_index <- which(apply(detail_annotation[, c("CDS", "Intron", "fiveUTR", "threeUTR")], 1, any))
  detail_annotation[-genic_index, "Intergenic"] <- TRUE
  detail_annotation[genic_index, "genic"] <- TRUE

  start_end <- paste0(BiocGenerics::start(peaks), "-", BiocGenerics::end(peaks))
  annotation_table <- data.table::as.data.table(detail_annotation)
  annotation_table[, start_end := start_end]
  annotation_table[, annotation := annotation]
  annotation_table <- annotation_table[, (1:8) := NULL]

  stop_codon_table <- data.table::as.data.table(stop_state$detailGenomicAnnotation)
  stop_codon_table[, start_end := start_end]
  stop_codon_table[, annotation := stop_state$annotation]
  stop_codon_table <- stop_codon_table[, (1:8) := NULL]

  start_codon_table <- data.table::as.data.table(start_state$detailGenomicAnnotation)
  start_codon_table[, start_end := start_end]
  start_codon_table[, annotation := start_state$annotation]
  start_codon_table <- start_codon_table[, (1:9) := NULL]

  detail_file <- peak_distribution_write_detail_annotations(
    annotation_result = list(
      annotation_table = data.table::copy(annotation_table),
      stop_codon_table = data.table::copy(stop_codon_table),
      start_codon_table = data.table::copy(start_codon_table)
    ),
    source_file = source_file,
    output_dir = detail_output_dir,
    sample_index = sample_index
  )

  list(
    annotation = annotation,
    annotation_exon = exon_state$annotation,
    annotation_stopcodon = stop_state$annotation,
    annotation_startcodon = start_state$annotation,
    annotation_transcript = transcript_state$annotation,
    detailGenomicAnnotation = detail_annotation,
    annotation_table = annotation_table,
    stop_codon_table = stop_codon_table,
    start_codon_table = start_codon_table,
    detail_file = detail_file
  )
}
