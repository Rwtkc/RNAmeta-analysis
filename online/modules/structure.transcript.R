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
  list(
    ara_TAIR10 = list(
      id = "ara_TAIR10",
      name = "Arabidopsis thaliana (TAIR10)",
      cdna = "Arabidopsis_thaliana.TAIR10.cdna.all.fa",
      cdna_index = "Arabidopsis_thaliana.TAIR10.cdna.all.fa.fai",
      gff3 = "Arabidopsis_thaliana.TAIR10.51.gff3.gz"
    ),
    osa_IRGSP_1 = list(
      id = "osa_IRGSP_1",
      name = "Oryza sativa (IRGSP 1.0)",
      cdna = "Oryza_sativa.IRGSP-1.0.cdna.all.fa",
      cdna_index = "Oryza_sativa.IRGSP-1.0.cdna.all.fa.fai",
      gff3 = "Oryza_sativa.IRGSP-1.0.51.chr.gff3.gz"
    )
  )
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
  value <- utils::URLdecode(value)
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
  match <- regexec(pattern, attributes, ignore.case = TRUE)
  parts <- regmatches(attributes, match)
  vapply(parts, function(item) {
    if (length(item) >= 3L) item[[3]] else ""
  }, character(1))
}

structure_gff_parent_values <- function(attributes) {
  parent <- structure_gff_attribute_value(attributes, "Parent")
  strsplit(parent, ",", fixed = TRUE)
}

structure_read_gff3_table <- function(gff3_path) {
  con <- if (grepl("\\.gz$", gff3_path, ignore.case = TRUE)) {
    gzfile(gff3_path, open = "rt")
  } else {
    file(gff3_path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) {
    stop("The Structure GFF3 reference has no annotation rows.")
  }

  gff <- data.table::fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    data.table = TRUE,
    col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  )
  gff[, chrom := browser_canonical_chr_vector(chrom)]
  gff[, start := suppressWarnings(as.integer(start))]
  gff[, end := suppressWarnings(as.integer(end))]
  gff[is.finite(start) & is.finite(end) & end >= start & nzchar(chrom)]
}

structure_build_transcript_map <- function(spec) {
  gff3_path <- app_path("data", "fa", spec$gff3)
  gff <- structure_read_gff3_table(gff3_path)

  transcripts <- gff[tolower(type) %in% c("mrna", "transcript", "ncrna", "lnc_rna", "rrna", "trna", "snrna", "snorna", "mirna")]
  if (!nrow(transcripts)) {
    transcripts <- gff[grepl("(^|;)ID=transcript:", attributes, ignore.case = TRUE)]
  }
  if (nrow(transcripts)) {
    transcripts[, transcript_id := structure_strip_gff_prefix(structure_gff_attribute_value(attributes, "ID"))]
    transcripts[, gene_id := structure_strip_gff_prefix(structure_gff_attribute_value(attributes, "Parent"))]
    transcripts[, transcript_name := structure_strip_gff_prefix(structure_gff_attribute_value(attributes, "Name"))]
    transcripts[, transcript_biotype := structure_gff_attribute_value(attributes, "biotype")]
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
  exons[]
}

structure_load_transcript_map <- function(spec) {
  cache_path <- structure_transcript_cache_path(spec)
  source_paths <- file.path(app_path("data", "fa"), c(spec$gff3, spec$cdna, spec$cdna_index))
  cache_ready <- file.exists(cache_path) &&
    all(file.exists(source_paths)) &&
    file.info(cache_path)$mtime >= max(file.info(source_paths)$mtime, na.rm = TRUE)

  if (isTRUE(cache_ready)) {
    return(readRDS(cache_path))
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
  candidates <- unique(c(
    structure_id_keys(transcript_id),
    structure_id_keys(transcript_name)
  ))
  matches <- candidates[candidates %in% seqnames]
  if (length(matches)) {
    return(matches[[1]])
  }

  stop(sprintf("Transcript '%s' was not found in the selected cDNA FASTA.", transcript_id))
}

structure_map_bed_to_transcript <- function(selection, spec) {
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

  selected_transcript <- candidates[1]
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
