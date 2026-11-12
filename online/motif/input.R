motif_default_controls <- function() {
  list(
    engine = "meme",
    analysis_mode = "centered",
    sample_name = "",
    sites_subset = "all",
    motif_width = 5L,
    sequence_window = 9L,
    max_motifs = 6L,
    job_lookup = ""
  )
}

motif_analysis_mode_choices <- function() {
  c("Peak-centred window" = "centered", "Full peaks with matched background" = "matched_full")
}

motif_validate_analysis_mode <- function(value = NULL) {
  if (is.null(value)) return("centered")
  if (length(value) != 1L || is.na(value) || !value %in% unname(motif_analysis_mode_choices())) {
    stop("Select a supported motif analysis mode.")
  }
  as.character(value)
}

motif_validate_analysis_controls <- function(controls) {
  controls <- utils::modifyList(motif_default_controls(), controls)
  controls$analysis_mode <- motif_validate_analysis_mode(controls$analysis_mode)
  if (!controls$engine %in% unname(motif_engine_choices())) stop("Select a supported motif engine.")
  if (identical(controls$analysis_mode, "centered")) {
    controls$sequence_window <- motif_validate_sequence_window(controls$sequence_window,
      motif_effective_width(controls$engine, controls$motif_width))
  }
  controls
}

motif_engine_choices <- function() {
  c(
    "MEME" = "meme",
    "STREME" = "streme",
    "HOMER" = "homer"
  )
}

motif_sites_subset_choices <- function() {
  c(
    "All peaks" = "all",
    "Top 1000 by score" = "score_1000",
    "Top 2000 by score" = "score_2000",
    "Top 3000 by score" = "score_3000",
    "Top 1000 by peak length" = "length_1000",
    "Top 2000 by peak length" = "length_2000",
    "Top 3000 by peak length" = "length_3000"
  )
}

motif_safe_integer <- function(value, fallback = 0L) {
  numeric_value <- suppressWarnings(as.integer(value))
  if (!is.finite(numeric_value) || length(numeric_value) != 1L || is.na(numeric_value)) {
    return(as.integer(fallback))
  }

  as.integer(numeric_value)
}

motif_min_width <- function(engine) {
  engine_key <- tolower(as.character(engine)[[1]])

  if (identical(engine_key, "streme")) {
    return(3L)
  }

  if (identical(engine_key, "meme")) {
    return(2L)
  }

  1L
}

motif_effective_width <- function(engine, motif_width) {
  max(motif_min_width(engine), motif_safe_integer(motif_width, motif_default_controls()$motif_width))
}

motif_normalize_subset_value <- function(subset_value = NULL) {
  raw_value <- subset_value
  if (is.null(raw_value)) {
    raw_value <- motif_default_controls()$sites_subset
  }
  raw_value <- as.character(raw_value)[[1]]

  if (!nzchar(raw_value) || identical(raw_value, "All")) {
    return("all")
  }

  if (raw_value %in% c("1000", "2000", "3000")) {
    return(sprintf("score_%s", raw_value))
  }

  valid_values <- unname(motif_sites_subset_choices())
  if (raw_value %in% valid_values) {
    return(raw_value)
  }

  motif_default_controls()$sites_subset
}

motif_subset_spec <- function(subset_value = NULL) {
  normalized_value <- motif_normalize_subset_value(subset_value)

  if (identical(normalized_value, "all")) {
    return(list(
      value = "all",
      mode = "all",
      size = Inf
    ))
  }

  parts <- strsplit(normalized_value, "_", fixed = TRUE)[[1]]
  list(
    value = normalized_value,
    mode = parts[[1]],
    size = motif_safe_integer(parts[[2]], 0L)
  )
}

motif_subset_peak_ranges <- function(peak_gr, subset_value = "All") {
  if (is.null(peak_gr) || !length(peak_gr)) {
    return(peak_gr)
  }

  subset_spec <- motif_subset_spec(subset_value)
  if (identical(subset_spec$mode, "all")) {
    return(peak_gr)
  }

  subset_n <- motif_safe_integer(subset_spec$size, length(peak_gr))
  if (!is.finite(subset_n) || subset_n <= 0) {
    return(peak_gr)
  }

  ordered_peaks <- peak_gr
  if (identical(subset_spec$mode, "score") && "score" %in% names(S4Vectors::mcols(peak_gr))) {
    score_values <- suppressWarnings(as.numeric(S4Vectors::mcols(peak_gr)$score))
    if (any(!is.na(score_values))) {
      ordered_peaks <- peak_gr[order(score_values, decreasing = TRUE, na.last = TRUE)]
    }
  } else if (identical(subset_spec$mode, "length")) {
    peak_lengths <- as.numeric(GenomicRanges::width(peak_gr))
    ordered_peaks <- peak_gr[order(peak_lengths, decreasing = TRUE, na.last = TRUE)]
  }

  ordered_peaks[seq_len(min(length(ordered_peaks), subset_n))]
}

motif_progress_update <- function(progress, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  safe_value <- as.integer(round(max(0, min(100, value))))
  progress$set(
    value = safe_value,
    message = "Running Motif Analysis",
    detail = sprintf("%s%% | %s", safe_value, detail)
  )

  invisible(NULL)
}

motif_sequence_window_width <- function() {
  motif_default_controls()$sequence_window
}

motif_validate_sequence_window <- function(value = NULL, motif_width = 1L) {
  if (is.null(value)) value <- motif_sequence_window_width()
  window <- suppressWarnings(as.numeric(value))
  if (length(window) != 1L || is.na(window) || !is.finite(window) ||
      window < 1 || window > .Machine$integer.max || window != floor(window)) {
    stop("Sequence window must be a positive integer in nt.")
  }
  if (window < motif_width) {
    stop(sprintf("Sequence window (%d nt) must be at least the motif width (%d nt).", as.integer(window), motif_width))
  }
  as.integer(window)
}

motif_genome_specs <- function() {
  list(
    hg38 = list(package = "BSgenome.Hsapiens.UCSC.hg38", object = "Hsapiens"),
    hg19 = list(package = "BSgenome.Hsapiens.UCSC.hg19", object = "Hsapiens"),
    mm10 = list(package = "BSgenome.Mmusculus.UCSC.mm10", object = "Mmusculus"),
    mm9 = list(package = "BSgenome.Mmusculus.UCSC.mm9", object = "Mmusculus"),
    rn6 = list(package = "BSgenome.Rnorvegicus.UCSC.rn6", object = "Rnorvegicus"),
    cel_WBcel235 = list(package = "BSgenome.Celegans.lq.celWBcel235", object = "Celegans"),
    sce_R64 = list(package = "BSgenome.Scerevisiae.lq.sceR64", object = "Scerevisiae"),
    dme_BDGP6 = list(package = "BSgenome.Dmelanogaster.lq.dmeBDGP6", object = "Dmelanogaster"),
    ara_TAIR10 = list(package = "BSgenome.Athaliana.lq.araTAIR10", object = "Athaliana"),
    bdi_v3 = list(package = "BSgenome.Bdistachyon.lq.bdiv3", object = "Bdistachyon"),
    dre_GRCz11 = list(package = "BSgenome.Drerio.lq.dreGRCz11", object = "Drerio"),
    gga_GRCg6a = list(package = "BSgenome.Ggallus.lq.ggaGRCg6a", object = "Ggallus"),
    ggo_gorGor4 = list(package = "BSgenome.Ggorilla.lq.ggogorGor4", object = "Ggorilla"),
    gma_v2 = list(package = "BSgenome.Gmax.lq.gmav2", object = "Gmax"),
    mmu_Mmul_8 = list(package = "BSgenome.Mmulatta.lq.mmuMmul8", object = "Mmulatta"),
    mtr_MedtrA17_4 = list(package = "BSgenome.Mtruncatula.lq.mtrMedtrA174", object = "Mtruncatula"),
    osa_IRGSP_1 = list(package = "BSgenome.Osativa.lq.osaIRGSP1", object = "Osativa"),
    ptc_v3 = list(package = "BSgenome.Ptrichocarpa.lq.ptcv3", object = "Ptrichocarpa"),
    ptr_v3 = list(package = "BSgenome.Ptroglodytes.lq.ptrv3", object = "Ptroglodytes"),
    sbi_v3 = list(package = "BSgenome.Sbicolor.lq.sbiv3", object = "Sbicolor"),
    ssc_Sscrofa11 = list(package = "BSgenome.Sscrofa.lq.sscSscrofa11", object = "Sscrofa"),
    zma_RefGen_v4 = list(package = "BSgenome.Zmays.lq.zmaRefGenv4", object = "Zmays"),
    osa_rufipogon = list(package = "BSgenome.Orufipogon.lq.oryzaORW1943", object = "Orufipogon"),
    vvi_12X = list(package = "BSgenome.Vvinifera.lq.vvi12X", object = "Vvinifera"),
    ahy_tifrunner2 = list(package = "BSgenome.Ahypogaea.lq.arahyTifrunner.gnm2", object = "Ahypogaea")
  )
}

motif_species_genome_code <- function(species) {
  species_value <- as.character(species)[[1]]
  available_codes <- names(motif_genome_specs())

  if (!species_value %in% available_codes && is.null(motif_fasta_path(species_value))) {
    stop(sprintf(
      "Unsupported motif species code: %s",
      species_value
    ))
  }

  species_value
}

motif_load_seqdb <- function(species) {
  genome_code <- motif_species_genome_code(species)
  fasta <- motif_fasta_path(genome_code)
  if (!is.null(fasta)) {
    index <- Rsamtools::scanFaIndex(fasta)
    original <- as.character(GenomicRanges::seqnames(index))
    normalized <- ifelse(startsWith(original, "chr"), original, paste0("chr", original))
    if (anyDuplicated(normalized)) stop("FASTA sequence names are ambiguous after chromosome normalization.")
    return(list(fasta_path = fasta,
      lengths = stats::setNames(IRanges::width(index), normalized),
      original_names = stats::setNames(original, normalized)))
  }
  spec <- motif_genome_specs()[[genome_code]]

  if (is.null(spec)) {
    stop(sprintf("Motif sequence extraction is not configured for species: %s", species))
  }

  if (!requireNamespace(spec$package, quietly = TRUE)) {
    stop(sprintf(
      "Required BSgenome package '%s' is not installed in this R environment.",
      spec$package
    ))
  }

  getExportedValue(spec$package, spec$object)
}

motif_fasta_path <- function(species) {
  roots <- unique(c(Sys.getenv("RNAMETA_FASTA_ROOT", unset = ""), app_path("data", "fa")))
  roots <- roots[nzchar(roots)]
  candidates <- unlist(lapply(roots, function(root) file.path(root, paste0(species, c(".fa", ".fasta")))))
  found <- candidates[file.exists(candidates) & file.exists(paste0(candidates, ".fai"))]
  if (length(found)) normalizePath(found[[1]], winslash = "/", mustWork = TRUE) else NULL
}

motif_reference_lengths <- function(genome) {
  if (is.list(genome) && !is.null(genome$fasta_path)) return(genome$lengths)
  GenomeInfoDb::seqlengths(genome)
}

motif_reference_sequences <- function(genome, regions) {
  if (!is.list(genome) || is.null(genome$fasta_path)) return(Biostrings::getSeq(genome, regions))
  original <- unname(genome$original_names[as.character(GenomicRanges::seqnames(regions))])
  if (anyNA(original)) stop("Sequence extraction requested an unknown FASTA chromosome.")
  query <- GenomicRanges::GRanges(original, IRanges::ranges(regions), strand = GenomicRanges::strand(regions))
  sequences <- Biostrings::getSeq(Rsamtools::FaFile(genome$fasta_path), param = query)
  names(sequences) <- names(regions)
  sequences
}

motif_extract_sequences <- function(query_regions, species, prefix = "seq") {
  seq_db <- motif_load_seqdb(species)
  chr_len <- motif_reference_lengths(seq_db)
  valid_seqlevels <- names(chr_len)
  query_regions <- query_regions[as.character(GenomicRanges::seqnames(query_regions)) %in% valid_seqlevels]
  query_regions <- GenomeInfoDb::keepSeqlevels(
    query_regions,
    value = intersect(GenomeInfoDb::seqlevels(query_regions), valid_seqlevels),
    pruning.mode = "coarse"
  )

  genome_range <- GenomicRanges::GRanges(
    seqnames = names(chr_len),
    ranges = IRanges::IRanges(1, width = chr_len)
  )
  query_regions <- IRanges::subsetByOverlaps(
    query_regions,
    genome_range,
    type = "within",
    ignore.strand = TRUE
  )

  if (!length(query_regions)) {
    stop("No motif regions remained after filtering against the selected genome.")
  }

  sequences <- motif_reference_sequences(seq_db, query_regions)

  names(sequences) <- sprintf("%s_%d", prefix, seq_along(sequences))
  names(query_regions) <- names(sequences)
  list(regions = query_regions, sequences = sequences)
}

motif_shuffle_regions <- function(peak_gr, species) {
  seq_db <- motif_load_seqdb(species)
  chromosome_lengths <- motif_reference_lengths(seq_db)
  peak_gr <- peak_gr[as.character(GenomicRanges::seqnames(peak_gr)) %in% names(chromosome_lengths)]

  if (!length(peak_gr)) {
    stop("No motif regions remained after matching BED peaks to the selected genome.")
  }

  chr_lens <- chromosome_lengths[as.character(GenomeInfoDb::seqnames(peak_gr))]
  max_starts <- chr_lens - IRanges::width(peak_gr) + 1L
  if (anyNA(max_starts) || any(max_starts < 1)) stop("Motif control regions exceed chromosome lengths.")
  starts <- vapply(max_starts, function(limit) sample.int(limit, 1L), integer(1))
  GenomicRanges::GRanges(
    seqnames = as.character(GenomeInfoDb::seqnames(peak_gr)),
    ranges = IRanges::IRanges(starts, width = IRanges::width(peak_gr)),
    strand = GenomicRanges::strand(peak_gr)
  )
}

motif_center_regions <- function(regions, width = motif_sequence_window_width()) {
  if (is.null(regions) || !length(regions)) return(regions)
  GenomicRanges::resize(regions, width = motif_validate_sequence_window(width), fix = "center")
}

motif_write_bed <- function(gr, path) {
  bed_df <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(gr)),
    chromStart = pmax(IRanges::start(gr) - 1L, 0L),
    chromEnd = IRanges::end(gr),
    name = if (is.null(names(gr))) sprintf("region_%d", seq_along(gr)) else names(gr),
    score = 0,
    strand = as.character(GenomicRanges::strand(gr)),
    stringsAsFactors = FALSE
  )

  utils::write.table(
    bed_df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

motif_write_unix_lines <- function(path, lines) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(paste(c(lines, ""), collapse = "\n"))), con)
}

motif_sequence_width_summary <- function(sequences) {
  sequence_widths <- as.integer(Biostrings::width(sequences))
  sequence_widths <- sequence_widths[is.finite(sequence_widths) & !is.na(sequence_widths)]

  if (!length(sequence_widths)) {
    return(list(
      min_width = NA_integer_,
      max_width = NA_integer_
    ))
  }

  list(
    min_width = min(sequence_widths),
    max_width = max(sequence_widths)
  )
}

motif_sequence_engine_min_width <- function(engine) {
  engine_key <- tolower(as.character(engine)[[1]])

  if (identical(engine_key, "streme")) {
    return(3L)
  }

  if (identical(engine_key, "meme")) {
    return(2L)
  }

  1L
}

motif_sequence_width_recommendation <- function(engine, sequence_label, width_summary) {
  engine_key <- toupper(as.character(engine)[[1]])
  min_allowed_width <- motif_sequence_engine_min_width(engine)

  if (is.na(width_summary$min_width) || is.na(width_summary$max_width)) {
    return("No sequence-length summary was available for this motif job.")
  }

  if (width_summary$max_width <= 1L) {
    return(sprintf(
      paste(
        "The extracted %s sequences are too short for motif discovery (%d-%d nt).",
        "Please upload a wider peak BED file, or a BED expanded around the summit before running motif analysis."
      ),
      sequence_label,
      width_summary$min_width,
      width_summary$max_width
    ))
  }

  if (width_summary$min_width < min_allowed_width) {
    return(sprintf(
      "%s requires motif width >= %d nt, but observed %s sequence lengths were %d-%d nt. Try a different BED subset or another motif engine.",
      engine_key,
      min_allowed_width,
      sequence_label,
      width_summary$min_width,
      width_summary$max_width
    ))
  }

  sprintf(
    "Recommended motif width to keep all %s sequences for %s: %d-%d nt.",
    sequence_label,
    engine_key,
    min_allowed_width,
    width_summary$min_width
  )
}

motif_filter_sequences_by_width <- function(regions, sequences, motif_width, sequence_label = "query", engine = "homer") {
  if (length(regions) != length(sequences)) stop("Motif regions and sequences must be paired before width filtering.")
  motif_width <- as.integer(motif_width[[1]])
  width_summary <- motif_sequence_width_summary(sequences)
  keep_index <- Biostrings::width(sequences) >= motif_width

  if (!any(keep_index)) {
    stop(sprintf(
      paste(
        "All %s sequences were shorter than the selected motif width (%d nt).",
        "Observed %s sequence lengths: %d-%d nt.",
        "%s"
      ),
      sequence_label,
      motif_width,
      sequence_label,
      width_summary$min_width,
      width_summary$max_width,
      motif_sequence_width_recommendation(engine, sequence_label, width_summary)
    ))
  }

  list(
    regions = regions[keep_index],
    sequences = sequences[keep_index]
  )
}
