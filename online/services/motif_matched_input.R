motif_write_shared_background <- function(sequences, path) {
  counts <- colSums(Biostrings::alphabetFrequency(sequences, baseOnly = TRUE)[, c("A", "C", "G", "T"), drop = FALSE])
  frequencies <- (counts + 1) / (sum(counts) + 4)
  motif_write_unix_lines(path, c("# Shared zero-order query background; add-one pseudocount per base",
    paste(names(frequencies), sprintf("%.17g", frequencies))))
  list(source = "query.fa", order = 0L, pseudocountPerBase = 1L,
    sharedByQueryAndControl = TRUE, frequencies = as.list(frequencies), file = basename(path))
}

motif_prepare_matched_input <- function(regions, result, annotation_bundle, exclusion, package_dir, progress = NULL) {
  genome <- motif_load_seqdb(result$species)
  chromosome_lengths <- motif_reference_lengths(genome)
  chrom <- as.character(GenomicRanges::seqnames(regions))
  valid_genome <- chrom %in% names(chromosome_lengths) & IRanges::start(regions) >= 1L &
    IRanges::end(regions) <= chromosome_lengths[chrom]
  valid_genome[is.na(valid_genome)] <- FALSE
  known_strand <- as.character(GenomicRanges::strand(regions)) %in% c("+", "-")
  wide_enough <- IRanges::width(regions) >= result$motif_width
  audit <- data.frame(input_index = seq_along(regions), chrom = chrom,
    start = IRanges::start(regions) - 1L, end = IRanges::end(regions),
    strand = as.character(GenomicRanges::strand(regions)),
    reason = ifelse(!known_strand, "unknown_RNA_strand",
      ifelse(!valid_genome, "outside_selected_genome",
        ifelse(!wide_enough, "shorter_than_motif_width", "pending_match"))))
  write_audit <- function() utils::write.table(audit, file.path(package_dir, "input_audit.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  keep <- which(audit$reason == "pending_match")
  if (!length(keep)) {
    write_audit()
    stop("No full peaks are eligible. Matched analysis requires BED RNA strands (+/-), valid genome coordinates, and peaks at least as wide as the motif.")
  }
  query <- regions[keep]
  names(query) <- sprintf("query_%d", keep)
  sequences <- motif_reference_sequences(genome, query)
  names(sequences) <- names(query)
  annotation <- motif_matching_annotation(annotation_bundle$txdb, names(chromosome_lengths))
  exclusion <- exclusion[as.character(GenomicRanges::seqnames(exclusion)) %in% names(chromosome_lengths)]
  exclusion <- GenomicRanges::reduce(exclusion, ignore.strand = TRUE)
  matched <- motif_match_background(query, sequences, annotation, genome, exclusion, progress)
  audit$reason[keep] <- matched$queries$reason
  write_audit()
  utils::write.table(matched$queries, file.path(package_dir, "matching_queries.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  controls <- matched$controls
  if (!nrow(controls)) stop("No non-peak background matched the fixed gene, feature and GC criteria. See input_audit.tsv for excluded peaks.")
  controls$chrom <- as.character(GenomicRanges::seqnames(query))[controls$query]
  controls$strand <- as.character(GenomicRanges::strand(query))[controls$query]
  controls$control_id <- sprintf("control_%d", seq_len(nrow(controls)))
  pair_table <- controls[, setdiff(names(controls), "sequence"), with = FALSE]
  utils::write.table(pair_table, file.path(package_dir, "matching_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  first <- controls[!duplicated(controls$query)]
  control_regions <- GenomicRanges::GRanges(first$chrom, IRanges::IRanges(first$start, first$end), strand = first$strand)
  names(control_regions) <- first$control_id
  control_sequences <- Biostrings::DNAStringSet(first$sequence)
  names(control_sequences) <- names(control_regions)
  query <- query[first$query]; sequences <- sequences[first$query]
  if (!all(IRanges::width(query) == Biostrings::width(control_sequences))) stop("Matched motif input lost length pairing.")
  summary <- list(inputCount = length(regions), eligibleCount = length(keep),
    matchedCount = length(query), matchedFraction = length(query) / length(regions),
    excludedReasons = as.list(table(audit$reason[audit$reason != "matched"])),
    config = matched$config, exclusionScope = "Union of all saved BED files, before subset selection; either strand",
    annotationMd5 = as.list(tools::md5sum(annotation_bundle$txdb_path)),
    selection = "No motif counts used; only peaks with at least one accepted background enter discovery",
    orientation = "RNA strand from BED; unknown strands excluded; no reverse-complement search")
  list(query_regions = query, query_sequences = sequences, control_regions = control_regions,
    control_sequences = control_sequences, matching_summary = summary,
    extra_files = c("input_audit.tsv", "matching_queries.tsv", "matching_pairs.tsv"))
}
