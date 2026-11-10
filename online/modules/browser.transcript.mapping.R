browser_transcript_reference_spec <- function(species) {
  if (!exists("structure_transcript_reference_spec", mode = "function") ||
      !exists("structure_transcript_reference_ready", mode = "function") ||
      !exists("structure_load_transcript_map", mode = "function")) {
    return(NULL)
  }

  spec <- structure_transcript_reference_spec(species)
  if (is.null(spec) || !structure_transcript_reference_ready(spec)) {
    return(NULL)
  }

  spec
}

browser_empty_transcript_table <- function() {
  data.table::data.table(
    row_key = character(),
    source_file = character(),
    source_row = integer(),
    transcript_interval = character(),
    transcript_id = character(),
    gene_id = character(),
    transcript_length = integer(),
    igv_locus = character(),
    overlap_bases = integer()
  )
}

browser_react_transcript_table_payload <- function(table) {
  empty_payload <- list(
    columns = list(),
    rows = list(),
    total = 0L
  )
  if (is.null(table) || !nrow(table)) {
    return(empty_payload)
  }

  bed_cols <- grep("^bed_col_[0-9]+$", names(table), value = TRUE)
  bed_cols <- bed_cols[order(as.integer(sub("^bed_col_", "", bed_cols)))]
  first_bed_cols <- bed_cols[seq_len(min(3L, length(bed_cols)))]
  remaining_bed_cols <- if (length(bed_cols) > 3L) bed_cols[4:length(bed_cols)] else character()
  display_cols <- c(first_bed_cols, "transcript_interval", "transcript_id", "igv_link", remaining_bed_cols)
  display_cols <- display_cols[display_cols %in% c(names(table), "igv_link")]
  bed_headers <- browser_bed_column_names(length(bed_cols))
  headers <- c(
    bed_headers[seq_along(first_bed_cols)],
    "Transcript interval",
    "Transcript ID",
    "IGV link",
    if (length(remaining_bed_cols)) bed_headers[seq(4L, length(bed_cols))] else character()
  )
  headers <- headers[seq_along(display_cols)]

  columns <- lapply(seq_along(display_cols), function(index) {
    list(id = display_cols[[index]], label = headers[[index]])
  })
  rows <- lapply(seq_len(nrow(table)), function(index) {
    row <- table[index]
    values <- lapply(display_cols, function(column) {
      if (identical(column, "igv_link")) {
        has_locus <- nzchar(as.character(row$transcript_id[[1]] %||% "")) &&
          nzchar(as.character(row$igv_locus[[1]] %||% ""))
        return(if (has_locus) "Open in IGV" else "")
      }

      value <- row[[column]][[1]] %||% ""
      if (length(value) != 1L || is.na(value)) "" else as.character(value)
    })
    names(values) <- display_cols
    values$row_key <- as.character(row$row_key[[1]] %||% "")
    values$can_open <- nzchar(as.character(row$transcript_id[[1]] %||% "")) &&
      nzchar(as.character(row$igv_locus[[1]] %||% ""))
    values
  })

  list(columns = columns, rows = rows, total = as.integer(nrow(table)))
}

browser_transcript_locus <- function(transcript_id, transcript_length) {
  transcript_id <- as.character(transcript_id %||% "")[[1]]
  transcript_length <- suppressWarnings(as.integer(transcript_length))
  if (!nzchar(transcript_id) || !is.finite(transcript_length) || transcript_length < 1L) {
    return("")
  }

  sprintf("%s:1-%s", transcript_id, transcript_length)
}

browser_format_transcript_segments <- function(starts, ends) {
  starts <- suppressWarnings(as.integer(starts))
  ends <- suppressWarnings(as.integer(ends))
  keep <- is.finite(starts) & is.finite(ends) & starts >= 1L & ends >= starts
  if (!any(keep)) {
    return("")
  }

  starts <- starts[keep]
  ends <- ends[keep]
  ord <- order(starts, ends)
  starts <- starts[ord]
  ends <- ends[ord]

  paste(sprintf("%s-%s", starts, ends), collapse = ";")
}

browser_parse_transcript_segments <- function(interval_string) {
  interval_string <- trimws(as.character(interval_string %||% "")[[1]])
  if (!nzchar(interval_string)) {
    return(data.table::data.table(start = integer(), end = integer()))
  }

  pieces <- strsplit(interval_string, "[;,|]+", perl = TRUE)[[1]]
  pieces <- trimws(pieces)
  pieces <- pieces[nzchar(pieces)]
  if (!length(pieces)) {
    return(data.table::data.table(start = integer(), end = integer()))
  }

  parsed <- regmatches(pieces, regexec("^([0-9]+)-([0-9]+)$", pieces))
  starts <- integer()
  ends <- integer()
  for (match in parsed) {
    if (length(match) == 3L) {
      starts <- c(starts, suppressWarnings(as.integer(match[[2]])))
      ends <- c(ends, suppressWarnings(as.integer(match[[3]])))
    }
  }

  segments <- data.table::data.table(start = starts, end = ends)
  if (!nrow(segments)) {
    return(segments)
  }

  segments <- segments[is.finite(start) & is.finite(end) & end >= start]
  data.table::setorder(segments, start, end)
  segments[]
}

browser_prepare_transcript_exons <- function(transcript_spec) {
  exons <- data.table::copy(structure_load_transcript_map(transcript_spec))
  if (!nrow(exons)) {
    return(exons)
  }

  exons[, chrom := browser_canonical_chr_vector(chrom)]
  exons[, exon_width := as.integer(end - start + 1L)]
  exons[, sort_start := ifelse(strand == "-", -start, start)]
  exons[, sort_end := ifelse(strand == "-", -end, end)]
  data.table::setorder(exons, transcript_id, sort_start, sort_end)
  exons[, cumulative_before := cumsum(data.table::shift(exon_width, fill = 0L)), by = transcript_id]
  exons[, transcript_length := as.integer(sum(exon_width)), by = transcript_id]
  exons[]
}

browser_map_bed_to_transcript_rows <- function(bed, file_entry, file_index, transcript_exons) {
  if (ncol(bed) < 3 || !nrow(bed) || !nrow(transcript_exons)) {
    return(browser_empty_transcript_table())
  }

  original_col_count <- ncol(bed)
  original_cols <- data.table::data.table(row_id = seq_len(nrow(bed)))
  for (col_index in seq_len(original_col_count)) {
    original_cols[[paste0("bed_col_", col_index)]] <- as.character(bed[[col_index]])
  }

  bed_ranges <- data.table::data.table(
    row_id = seq_len(nrow(bed)),
    chrom = browser_canonical_chr_vector(bed[[1]]),
    start = suppressWarnings(as.integer(bed[[2]])) + 1L,
    end = suppressWarnings(as.integer(bed[[3]])),
    bed_name = if (original_col_count >= 4) as.character(bed[[4]]) else "",
    bed_strand = if (original_col_count >= 6) as.character(bed[[6]]) else "."
  )
  bed_ranges <- bed_ranges[is.finite(start) & is.finite(end) & end >= start & nzchar(chrom)]

  if (!nrow(bed_ranges)) {
    return(browser_empty_transcript_table())
  }

  overlap_exons <- data.table::copy(transcript_exons)
  data.table::setkey(overlap_exons, chrom, start, end)
  data.table::setkey(bed_ranges, chrom, start, end)
  hits <- data.table::foverlaps(
    bed_ranges,
    overlap_exons,
    by.x = c("chrom", "start", "end"),
    by.y = c("chrom", "start", "end"),
    type = "any",
    nomatch = 0L
  )

  if (!nrow(hits)) {
    return(browser_empty_transcript_table())
  }

  bed_start_col <- if ("i.start" %in% names(hits)) "i.start" else "start"
  bed_end_col <- if ("i.end" %in% names(hits)) "i.end" else "end"
  hits[, bed_start := get(bed_start_col)]
  hits[, bed_end := get(bed_end_col)]
  hits[, overlap_start := pmax(start, bed_start)]
  hits[, overlap_end := pmin(end, bed_end)]
  hits[, overlap_bases := pmax(0L, overlap_end - overlap_start + 1L)]
  hits <- hits[overlap_bases > 0L & nzchar(transcript_id)]

  if (!nrow(hits)) {
    return(browser_empty_transcript_table())
  }

  hits[, highlight_start := ifelse(
    strand == "-",
    cumulative_before + (end - overlap_end + 1L),
    cumulative_before + (overlap_start - start + 1L)
  )]
  hits[, highlight_end := ifelse(
    strand == "-",
    cumulative_before + (end - overlap_start + 1L),
    cumulative_before + (overlap_end - start + 1L)
  )]
  hits[, transcript_interval_start := as.integer(pmin(highlight_start, highlight_end))]
  hits[, transcript_interval_end := as.integer(pmax(highlight_start, highlight_end))]

  mapped <- hits[, .(
    overlap_bases = as.integer(sum(overlap_bases, na.rm = TRUE)),
    chrom = chrom[[1]],
    strand = strand[[1]],
    gene_id = gene_id[[1]] %||% "",
    transcript_name = transcript_name[[1]] %||% "",
    transcript_length = as.integer(transcript_length[[1]])
  ), by = .(row_id, transcript_id, transcript_interval_start, transcript_interval_end)]

  bed_meta <- bed_ranges[, .(row_id, bed_name, bed_strand)]
  mapped <- merge(mapped, bed_meta, by = "row_id", all.x = TRUE, sort = FALSE)
  mapped[, name_score := vapply(seq_len(.N), function(index) {
    if (!exists("structure_id_keys", mode = "function")) {
      return(0L)
    }
    bed_keys <- structure_id_keys(bed_name[[index]] %||% "")
    transcript_keys <- unique(c(
      structure_id_keys(transcript_id[[index]]),
      structure_id_keys(gene_id[[index]]),
      structure_id_keys(transcript_name[[index]])
    ))
    as.integer(length(intersect(bed_keys, transcript_keys)) > 0L)
  }, integer(1))]
  mapped[, strand_score := as.integer(bed_strand %in% c("+", "-") & bed_strand == strand)]
  mapped[, transcript_overlap_bases := as.integer(sum(overlap_bases, na.rm = TRUE)), by = .(row_id, transcript_id)]
  data.table::setorder(mapped, row_id, -name_score, -strand_score, -transcript_overlap_bases, transcript_id, transcript_interval_start, transcript_interval_end)
  mapped[, segment_index := seq_len(.N), by = .(row_id, transcript_id)]
  mapped[, transcript_interval := sprintf("%s-%s", transcript_interval_start, transcript_interval_end)]
  mapped[, igv_locus := mapply(browser_transcript_locus, transcript_id, transcript_length, USE.NAMES = FALSE)]
  mapped[, row_key := sprintf("f%s_r%s_tx%s_seg%s", file_index, row_id, transcript_id, segment_index)]
  mapped[, source_file := file_entry$file_name %||% basename(file_entry$file_path)]
  mapped[, source_row := as.integer(row_id)]

  mapped <- merge(mapped, original_cols, by.x = "row_id", by.y = "row_id", all.x = TRUE, sort = FALSE)
  data.table::setorder(mapped, row_id, -name_score, -strand_score, -transcript_overlap_bases, transcript_id, segment_index)

  keep <- c(
    "row_key", "source_file", "source_row", "transcript_interval", "transcript_id",
    "gene_id", "transcript_length", "igv_locus", "overlap_bases",
    paste0("bed_col_", seq_len(original_col_count))
  )
  mapped[, keep, with = FALSE]
}
