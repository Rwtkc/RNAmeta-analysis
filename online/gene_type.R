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

gene_type_load_annotation_bundle <- function(species) {
  if (exists("gene_type_cached_annotation_bundle", mode = "function")) {
    return(gene_type_cached_annotation_bundle(species))
  }

  gene_type_load_annotation_bundle_uncached(species)
}

gene_type_load_annotation_bundle_uncached <- function(species) {
  bundle <- peak_distribution_resolve_annotation_bundle(species)
  txdb <- AnnotationDbi::loadDb(bundle$txdb_path)
  GenomeInfoDb::seqlevels(txdb) <- sub("^([^chr])", "chr\\1", GenomeInfoDb::seqlevels(txdb))

  gff <- peak_distribution_load_rda_object(bundle$gff_path)

  if ("gene_biotype" %in% names(S4Vectors::mcols(gff))) {
    names(S4Vectors::mcols(gff))[names(S4Vectors::mcols(gff)) == "gene_biotype"] <- "gene_type"
  }

  if ("transcript_biotype" %in% names(S4Vectors::mcols(gff))) {
    names(S4Vectors::mcols(gff))[names(S4Vectors::mcols(gff)) == "transcript_biotype"] <- "transcript_type"
  }

  GenomeInfoDb::seqlevels(gff) <- sub("^([^chr])", "chr\\1", GenomeInfoDb::seqlevels(gff))

  list(
    species = species,
    species_abb = bundle$species_abb,
    txdb = txdb,
    gff = gff,
    txdb_path = bundle$txdb_path,
    gff_path = bundle$gff_path
  )
}

gene_type_prepare_feature_ranges <- function(gff_data, txdb) {
  transcript_type <- as.character(S4Vectors::mcols(gff_data)$transcript_type)
  keep <- !is.na(transcript_type) & nzchar(transcript_type)

  valid_seqlevels <- intersect(GenomeInfoDb::seqlevels(gff_data), peak_distribution_valid_chromosomes(txdb))
  if (length(valid_seqlevels)) {
    seq_keep <- as.character(GenomeInfoDb::seqnames(gff_data)) %in% valid_seqlevels
    keep <- keep & seq_keep
  }

  feature_ranges <- gff_data[keep]
  if (!length(feature_ranges)) {
    return(feature_ranges)
  }

  if (length(valid_seqlevels)) {
    feature_ranges <- GenomeInfoDb::keepSeqlevels(feature_ranges, valid_seqlevels, pruning.mode = "coarse")
  }

  S4Vectors::mcols(feature_ranges) <- S4Vectors::DataFrame(
    transcript_type = as.character(S4Vectors::mcols(feature_ranges)$transcript_type)
  )

  feature_ranges
}

gene_type_feature_ranges <- function(species, annotation_bundle = NULL) {
  if (exists("gene_type_cached_feature_ranges", mode = "function")) {
    return(gene_type_cached_feature_ranges(species))
  }

  if (is.null(annotation_bundle)) {
    annotation_bundle <- gene_type_load_annotation_bundle(species)
  }

  gene_type_prepare_feature_ranges(annotation_bundle$gff, annotation_bundle$txdb)
}

gene_type_read_peaks <- function(source_file, txdb, sample_name = "Group1") {
  peaks <- peak_distribution_read_peaks(source_file, txdb, sample_name = sample_name)
  GenomicRanges::resize(peaks, width = 1, fix = "center")
}

gene_type_find_overlapping_type_values <- function(query_regions, feature_ranges) {
  overlaps <- GenomicRanges::findOverlaps(query_regions, feature_ranges)

  if (!length(overlaps)) {
    return(character())
  }

  transcript_type <- S4Vectors::mcols(feature_ranges)$transcript_type[S4Vectors::subjectHits(overlaps)]
  gene_type_normalize_type_values(transcript_type)
}

gene_type_default_controls <- function() {
  list()
}

gene_type_normalize_type_values <- function(type_values) {
  type_values <- as.character(type_values)
  type_values <- type_values[!is.na(type_values) & nzchar(type_values)]
  type_values
}

gene_type_runtime_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$annotation_bundle <- new.env(hash = TRUE, parent = emptyenv())
  env$feature_ranges <- new.env(hash = TRUE, parent = emptyenv())
  env$analysis_result <- new.env(hash = TRUE, parent = emptyenv())
  env
})

gene_type_runtime_cache_key <- function(...) {
  paste(vapply(list(...), as.character, character(1)), collapse = "::")
}

gene_type_cached_annotation_bundle <- function(species) {
  key <- gene_type_runtime_cache_key(species)

  if (exists(key, envir = gene_type_runtime_cache_env$annotation_bundle, inherits = FALSE)) {
    return(get(key, envir = gene_type_runtime_cache_env$annotation_bundle, inherits = FALSE))
  }

  bundle <- gene_type_load_annotation_bundle_uncached(species)
  assign(key, bundle, envir = gene_type_runtime_cache_env$annotation_bundle)
  bundle
}

gene_type_cached_feature_ranges <- function(species) {
  key <- gene_type_runtime_cache_key(species)

  if (exists(key, envir = gene_type_runtime_cache_env$feature_ranges, inherits = FALSE)) {
    return(get(key, envir = gene_type_runtime_cache_env$feature_ranges, inherits = FALSE))
  }

  bundle <- gene_type_cached_annotation_bundle(species)
  feature_ranges <- gene_type_prepare_feature_ranges(bundle$gff, bundle$txdb)
  assign(key, feature_ranges, envir = gene_type_runtime_cache_env$feature_ranges)
  feature_ranges
}
