structure_extract_sequence <- function(upload_context, selection, window_size = "fullsize", preferred_transcript_id = NULL) {
  spec <- structure_transcript_reference_spec(upload_context$species)
  if (is.null(spec) || !structure_transcript_reference_ready(spec)) {
    stop("The selected species does not have ready cDNA FASTA and GFF3 references for Structure.")
  }

  interval <- selection$interval
  transcript <- structure_map_bed_to_transcript(selection, spec, preferred_transcript_id = preferred_transcript_id)
  fasta_path <- app_path("data", "fa", spec$cdna)
  seqname <- structure_match_transcript_seqname(
    transcript_id = transcript$transcript_id,
    transcript_name = transcript$transcript_name,
    fasta_path = fasta_path
  )
  fasta_index <- Rsamtools::scanFaIndex(fasta_path)
  seq_widths <- IRanges::width(fasta_index)
  seq_names <- as.character(GenomeInfoDb::seqnames(fasta_index))
  seq_length <- as.integer(seq_widths[match(seqname, seq_names)])

  range <- GenomicRanges::GRanges(seqnames = seqname, ranges = IRanges::IRanges(start = 1L, end = seq_length))
  dna <- Rsamtools::scanFa(Rsamtools::FaFile(fasta_path), param = range)[[1]]
  if (length(dna) == 0) {
    stop("No transcript sequence was extracted for the selected BED interval.")
  }

  sequence <- toupper(as.character(dna))
  sequence <- chartr("T", "U", sequence)

  list(
    sequence = sequence,
    locus = sprintf("%s:%s-%s", transcript$transcript_id, 1L, nchar(sequence)),
    bed_locus = sprintf("%s:%s-%s", interval$chrom, transcript$bed_start_1based, transcript$bed_end_1based),
    strand = transcript$strand %||% ".",
    window_size = "fullsize",
    interval_length = as.integer(interval$end - interval$start),
    used_bed_interval = FALSE,
    matched_bed_row = isTRUE(selection$matched_bed_row),
    transcript = transcript,
    highlight_ranges = transcript$highlight_ranges,
    highlight_start = transcript$highlight_start,
    highlight_end = transcript$highlight_end,
    transcript_length = transcript$transcript_length,
    transcript_seqname = seqname
  )
}

structure_rnafold_path <- function() {
  candidates <- c(
    Sys.getenv("RNAMETA_RNAFOLD", ""),
    Sys.which("RNAfold")
  )
  candidates <- candidates[nzchar(candidates)]
  matches <- candidates[file.exists(path.expand(candidates))]

  if (length(matches)) {
    return(path.expand(matches[[1]]))
  }

  ""
}

structure_fold_with_rnafold <- function(sequence) {
  rnafold <- structure_rnafold_path()
  if (!nzchar(rnafold)) {
    return(NULL)
  }

  input <- tempfile(pattern = "rnameta_structure_", fileext = ".fa")
  writeLines(c(">rnameta", sequence), input)
  on.exit(unlink(input, force = TRUE), add = TRUE)

  output <- tryCatch(
    system2(rnafold, args = "--noPS", stdin = input, stdout = TRUE, stderr = TRUE),
    error = function(...) character()
  )
  structure_line <- output[grepl("^[().]+\\s", output)]
  if (!length(structure_line)) {
    details <- paste(utils::tail(output, 6L), collapse = " | ")
    folded <- NA_character_
    attr(folded, "error") <- if (nzchar(details)) details else "RNAfold produced no parseable output."
    return(folded)
  }

  sub("\\s.*$", "", structure_line[[1]])
}

structure_fold_sequence <- function(sequence) {
  folded <- structure_fold_with_rnafold(sequence)
  if (!is.null(folded) && nchar(folded) == nchar(sequence)) {
    return(list(structure = folded, method = "RNAfold"))
  }

  rnafold_error <- attr(folded, "error") %||% ""
  stop(sprintf(
    "RNAfold did not return a valid dot-bracket structure. %s",
    rnafold_error
  ))
}

structure_fasta_seqnames <- function(fasta_path) {
  index <- Rsamtools::scanFaIndex(fasta_path)
  as.character(GenomeInfoDb::seqnames(index))
}

structure_match_fasta_seqname <- function(chrom, fasta_path) {
  seqnames <- structure_fasta_seqnames(fasta_path)
  candidates <- unique(c(
    chrom,
    paste0("Chr", chrom),
    paste0("chr", chrom),
    sub("^chr", "", chrom, ignore.case = TRUE)
  ))
  matches <- candidates[candidates %in% seqnames]
  if (length(matches)) {
    return(matches[[1]])
  }

  stop(sprintf("Chromosome '%s' was not found in the selected reference FASTA.", chrom))
}

structure_transcript_reference_specs <- function() {
  specs <- list(
    ara_TAIR10 = list(
      id = "ara_TAIR10",
      name = "Arabidopsis thaliana (TAIR10)",
      cdna = "Arabidopsis_thaliana.TAIR10.cdna.all.fa",
      cdna_index = "Arabidopsis_thaliana.TAIR10.cdna.all.fa.fai",
      gff3 = "Arabidopsis_thaliana.TAIR10.51.gff3.gz"
    ),
    hg38.gencode_v46 = list(
      id = "hg38.gencode_v46",
      name = "Homo sapiens (GRCh38, GENCODE v46)",
      cdna = "hg38.transcripts.fa",
      cdna_index = "hg38.transcripts.fa.fai",
      gff3 = "gencode.v46.annotation.gff3.gz"
    ),
    hg19.gencode_v19 = list(
      id = "hg19",
      name = "Homo sapiens (GRCh37, GENCODE v19)",
      cdna = "hg19.transcripts.fa",
      cdna_index = "hg19.transcripts.fa.fai",
      gff3 = "gencode.v19.annotation.gff3.gz"
    ),
    osa_IRGSP_1 = list(
      id = "osa_IRGSP_1",
      name = "Oryza sativa (IRGSP 1.0)",
      cdna = "Oryza_sativa.IRGSP-1.0.51.transcripts.fa",
      cdna_index = "Oryza_sativa.IRGSP-1.0.51.transcripts.fa.fai",
      gff3 = "Oryza_sativa.IRGSP-1.0.51.chr.gff3.gz"
    )
  )
  specs$hg38 <- specs$hg38.gencode_v46
  specs$hg19 <- specs$hg19.gencode_v19
  specs
}

structure_transcript_reference_spec <- function(species) {
  specs <- structure_transcript_reference_specs()
  species <- as.character(species %||% "")
  if (!species %in% names(specs)) {
    return(NULL)
  }

  specs[[species]]
}

structure_transcript_reference_ready <- function(spec) {
  if (is.null(spec)) {
    return(FALSE)
  }

  file.exists(app_path("data", "fa", spec$cdna)) &&
    file.exists(app_path("data", "fa", spec$cdna_index)) &&
    file.exists(app_path("data", "fa", spec$gff3))
}

structure_transcript_cache_root <- function() {
  path <- app_path("cache", "structure_transcript_maps")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

structure_transcript_cache_path <- function(spec) {
  file.path(structure_transcript_cache_root(), paste0(spec$id, ".rds"))
}

structure_strip_gff_prefix <- function(value) {
  value <- trimws(as.character(value %||% ""))
  value[is.na(value)] <- ""
  encoded <- grepl("%", value, fixed = TRUE)
  value[encoded] <- utils::URLdecode(value[encoded])
  gsub("^(gene|transcript|mRNA|rna|exon|cds):", "", value, ignore.case = TRUE)
}

structure_id_keys <- function(value) {
  value <- structure_strip_gff_prefix(value)
  value <- value[nzchar(value)]
  if (!length(value)) {
    return(character())
  }

  unique(c(value, sub("\\.[0-9]+$", "", value)))
}

structure_gff_attribute_value <- function(attributes, key) {
  pattern <- sprintf("(^|;)%s=([^;]+)", key)
  match <- regexpr(pattern, attributes, ignore.case = TRUE)
  values <- rep("", length(attributes))
  keep <- !is.na(match) & match > 0L
  values[keep] <- sub(sprintf("^;?%s=", key), "", regmatches(attributes, match), ignore.case = TRUE)
  values
}

structure_gff_parent_values <- function(attributes) {
  parent <- structure_gff_attribute_value(attributes, "Parent")
  strsplit(parent, ",", fixed = TRUE)
}

structure_read_gff3_table <- function(gff3_path) {
  con <- if (grepl("\\.gz$", gff3_path, ignore.case = TRUE)) gzfile(gff3_path, "rt") else file(gff3_path, "rt")
  on.exit(close(con), add = TRUE)
  chunks <- list()
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines)) next
    chunk <- data.table::fread(
      text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t", quote = "", fill = TRUE,
      col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
    )
    chunks[[length(chunks) + 1L]] <- chunk[tolower(type) %in% c("exon", "mrna", "transcript", "ncrna", "lnc_rna", "rrna", "trna", "snrna", "snorna", "mirna") | grepl("(^|;)ID=transcript:", attributes)]
  }
  gff <- data.table::rbindlist(chunks)
  if (!nrow(gff)) {
    stop("The Structure GFF3 reference has no annotation rows.")
  }
  gff[, chrom := browser_canonical_chr_vector(chrom)]
  gff[, start := suppressWarnings(as.integer(start))]
  gff[, end := suppressWarnings(as.integer(end))]
  gff[is.finite(start) & is.finite(end) & end >= start & nzchar(chrom)]
}

structure_build_transcript_map <- function(spec) {
  gff3_path <- app_path("data", "fa", spec$gff3)
  gff <- structure_read_gff3_table(gff3_path)

  transcripts <- gff[tolower(type) %in% c("mrna", "transcript", "ncrna", "lnc_rna", "rrna", "trna", "snrna", "snorna", "mirna") | grepl("(^|;)ID=transcript:", attributes, ignore.case = TRUE)]
  if (!nrow(transcripts)) {
    transcripts <- gff[grepl("(^|;)ID=transcript:", attributes, ignore.case = TRUE)]
  }
  if (nrow(transcripts)) {
    transcripts[, transcript_id := structure_strip_gff_prefix(structure_gff_attribute_value(attributes, "ID"))]
    transcripts[, gene_id := structure_strip_gff_prefix(structure_gff_attribute_value(attributes, "Parent"))]
    transcripts[, transcript_name := structure_strip_gff_prefix(structure_gff_attribute_value(attributes, "Name"))]
    transcripts[, transcript_biotype := structure_gff_attribute_value(attributes, "biotype")]
    transcripts[!nzchar(transcript_name), transcript_name := structure_gff_attribute_value(attributes, "transcript_name")]
    transcripts[!nzchar(transcript_biotype), transcript_biotype := structure_gff_attribute_value(attributes, "transcript_type")]
    transcripts <- transcripts[nzchar(transcript_id)]
    transcripts <- transcripts[, .(
      gene_id = gene_id[[1]],
      transcript_name = transcript_name[[1]],
      transcript_biotype = transcript_biotype[[1]],
      tx_start = min(start),
      tx_end = max(end),
      chrom = chrom[[1]],
      strand = strand[[1]]
    ), by = .(transcript_id)]
  } else {
    transcripts <- data.table::data.table(
      transcript_id = character(),
      gene_id = character(),
      transcript_name = character(),
      transcript_biotype = character(),
      tx_start = integer(),
      tx_end = integer(),
      chrom = character(),
      strand = character()
    )
  }

  exons <- gff[tolower(type) == "exon"]
  if (!nrow(exons)) {
    stop("The Structure GFF3 reference has no exon rows.")
  }

  parent_values <- structure_gff_parent_values(exons$attributes)
  exons <- exons[rep(seq_len(nrow(exons)), lengths(parent_values))]
  exons[, transcript_id := structure_strip_gff_prefix(unlist(parent_values, use.names = FALSE))]
  exons[, exon_id := structure_strip_gff_prefix(structure_gff_attribute_value(attributes, "ID"))]
  exons[, exon_rank := suppressWarnings(as.integer(structure_gff_attribute_value(attributes, "rank")))]
  exons[is.na(exon_rank), exon_rank := suppressWarnings(as.integer(structure_gff_attribute_value(attributes, "exon_number")))]
  exons <- exons[nzchar(transcript_id), .(chrom, start, end, strand, transcript_id, exon_id, exon_rank)]

  if (!nrow(exons)) {
    stop("The Structure GFF3 reference has no transcript-linked exons.")
  }

  if (nrow(transcripts)) {
    exons <- merge(
      exons,
      transcripts[, .(transcript_id, gene_id, transcript_name, transcript_biotype, tx_start, tx_end)],
      by = "transcript_id",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    exons[, `:=`(
      gene_id = "",
      transcript_name = "",
      transcript_biotype = "",
      tx_start = NA_integer_,
      tx_end = NA_integer_
    )]
  }
  exons[, exon_width := end - start + 1L]
  attr(exons, "annotation_parser_version") <- 2L
  exons[]
}

structure_load_transcript_map <- function(spec) {
  cache_path <- structure_transcript_cache_path(spec)
  source_paths <- file.path(app_path("data", "fa"), c(spec$gff3, spec$cdna, spec$cdna_index))
  cache_ready <- file.exists(cache_path) &&
    all(file.exists(source_paths)) &&
    file.info(cache_path)$mtime >= max(file.info(source_paths)$mtime, na.rm = TRUE)

  if (isTRUE(cache_ready)) {
    cached_map <- readRDS(cache_path)
    if (identical(attr(cached_map, "annotation_parser_version"), 2L)) return(cached_map)
  }

  map <- structure_build_transcript_map(spec)
  saveRDS(map, cache_path)
  map
}

structure_transcript_fasta_seqnames <- function(fasta_path) {
  as.character(GenomeInfoDb::seqnames(Rsamtools::scanFaIndex(fasta_path)))
}

structure_match_transcript_seqname <- function(transcript_id, transcript_name, fasta_path) {
  seqnames <- structure_transcript_fasta_seqnames(fasta_path)
  seqname_ids <- sub("\\|.*$", "", seqnames)
  candidates <- unique(c(
    structure_id_keys(transcript_id),
    structure_id_keys(transcript_name)
  ))
  matches <- seqnames[seqname_ids %in% candidates]
  if (length(matches)) {
    return(matches[[1]])
  }

  stop(sprintf("Transcript '%s' was not found in the selected cDNA FASTA.", transcript_id))
}

structure_map_bed_to_transcript <- function(selection, spec, preferred_transcript_id = NULL) {
  interval <- selection$interval
  exons <- structure_load_transcript_map(spec)
  bed_start <- as.integer(interval$start) + 1L
  bed_end <- as.integer(interval$end)
  bed_name_keys <- structure_id_keys(interval$name %||% "")

  overlaps <- exons[
    chrom == interval$chrom &
      end >= bed_start &
      start <= bed_end
  ]
  if (!nrow(overlaps)) {
    stop("The selected BED interval does not overlap any annotated transcript exon.")
  }

  overlaps[, overlap_start := pmax(start, bed_start)]
  overlaps[, overlap_end := pmin(end, bed_end)]
  overlaps[, overlap_bases := pmax(0L, overlap_end - overlap_start + 1L)]
  overlaps <- overlaps[overlap_bases > 0L]

  if (!nrow(overlaps)) {
    stop("The selected BED interval does not overlap any annotated transcript exon.")
  }

  candidates <- overlaps[, .(
    overlap_bases = sum(overlap_bases),
    chrom = chrom[[1]],
    strand = strand[[1]],
    gene_id = gene_id[[1]] %||% "",
    transcript_name = transcript_name[[1]] %||% "",
    transcript_biotype = transcript_biotype[[1]] %||% "",
    tx_start = suppressWarnings(min(tx_start, na.rm = TRUE)),
    tx_end = suppressWarnings(max(tx_end, na.rm = TRUE))
  ), by = transcript_id]

  interval_strand <- as.character(interval$strand %||% ".")[[1]]
  if (interval_strand %in% c("+", "-")) {
    candidates[, strand_score := as.integer(strand == interval_strand)]
  } else {
    candidates[, strand_score := 0L]
  }
  candidates[, name_score := vapply(seq_len(.N), function(i) {
    keys <- unique(c(
      structure_id_keys(transcript_id[[i]]),
      structure_id_keys(gene_id[[i]]),
      structure_id_keys(transcript_name[[i]])
    ))
    as.integer(length(intersect(keys, bed_name_keys)) > 0L)
  }, integer(1))]
  data.table::setorder(candidates, -name_score, -strand_score, -overlap_bases, transcript_id)

  preferred_transcript_id <- trimws(as.character(preferred_transcript_id %||% ""))
  selected_transcript <- if (nzchar(preferred_transcript_id)) {
    candidates[transcript_id == preferred_transcript_id]
  } else {
    candidates[1]
  }
  if (!nrow(selected_transcript)) {
    stop("The selected transcript is no longer available for this BED interval. Map the BED files again.")
  }
  selected_transcript <- selected_transcript[1]
  transcript_exons <- exons[transcript_id == selected_transcript$transcript_id]
  if (!nrow(transcript_exons)) {
    stop("Transcript exon coordinates could not be resolved for Structure.")
  }

  if (identical(selected_transcript$strand[[1]], "-")) {
    data.table::setorder(transcript_exons, -start, -end)
  } else {
    data.table::setorder(transcript_exons, start, end)
  }

  transcript_exons[, cumulative_before := cumsum(data.table::shift(exon_width, fill = 0L))]
  hit_exons <- transcript_exons[end >= bed_start & start <= bed_end]
  hit_exons[, overlap_start := pmax(start, bed_start)]
  hit_exons[, overlap_end := pmin(end, bed_end)]

  if (identical(selected_transcript$strand[[1]], "-")) {
    hit_exons[, highlight_start := cumulative_before + (end - overlap_end + 1L)]
    hit_exons[, highlight_end := cumulative_before + (end - overlap_start + 1L)]
  } else {
    hit_exons[, highlight_start := cumulative_before + (overlap_start - start + 1L)]
    hit_exons[, highlight_end := cumulative_before + (overlap_end - start + 1L)]
  }

  hit_exons <- hit_exons[highlight_end >= highlight_start]
  if (!nrow(hit_exons)) {
    stop("The selected BED interval could not be projected onto transcript coordinates.")
  }

  highlight_ranges <- lapply(seq_len(nrow(hit_exons)), function(i) {
    list(
      start = as.integer(hit_exons$highlight_start[[i]]),
      end = as.integer(hit_exons$highlight_end[[i]])
    )
  })
  transcript_length <- sum(transcript_exons$exon_width)

  list(
    transcript_id = selected_transcript$transcript_id[[1]],
    transcript_name = selected_transcript$transcript_name[[1]] %||% "",
    gene_id = selected_transcript$gene_id[[1]] %||% "",
    biotype = selected_transcript$transcript_biotype[[1]] %||% "",
    chrom = selected_transcript$chrom[[1]],
    strand = selected_transcript$strand[[1]],
    transcript_start = suppressWarnings(as.integer(selected_transcript$tx_start[[1]])),
    transcript_end = suppressWarnings(as.integer(selected_transcript$tx_end[[1]])),
    transcript_length = as.integer(transcript_length),
    bed_start_1based = bed_start,
    bed_end_1based = bed_end,
    highlight_ranges = highlight_ranges,
    highlight_start = min(vapply(highlight_ranges, `[[`, integer(1), "start")),
    highlight_end = max(vapply(highlight_ranges, `[[`, integer(1), "end")),
    overlap_bases = as.integer(sum(hit_exons$overlap_end - hit_exons$overlap_start + 1L)),
    candidate_count = nrow(candidates)
  )
}
