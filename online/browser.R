browser_read_bed_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(data.table::data.table())
  }

  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    lines <- readLines(connection, warn = FALSE)
    if (!length(lines)) {
      return(data.table::data.table())
    }

    return(data.table::fread(
      text = paste(lines, collapse = "\n"),
      header = FALSE,
      sep = "\t",
      fill = TRUE,
      quote = "",
      data.table = TRUE
    ))
  }

  data.table::fread(
    path,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    quote = "",
    data.table = TRUE
  )
}

browser_canonical_chr_vector <- function(chromosomes) {
  if (is.null(chromosomes)) {
    return(character())
  }

  chr <- trimws(as.character(chromosomes))
  chr[is.na(chr)] <- ""
  normalized <- tolower(sub("^chr", "", chr, ignore.case = TRUE))
  result <- chr

  nuclear <- normalized %in% as.character(1:22)
  result[nuclear] <- normalized[nuclear]
  result[normalized %in% c("x", "y")] <- toupper(normalized[normalized %in% c("x", "y")])
  result[normalized %in% c("m", "mt", "mitochondria", "mitochondrion")] <- "Mt"
  result[normalized %in% c("c", "pt", "chloroplast", "plastid")] <- "Pt"

  result
}

browser_annotation_region_label <- function(annotation_view) {
  switch(
    as.character(annotation_view %||% "transcript"),
    five_utr = "5'UTR",
    cds = "CDS",
    three_utr = "3'UTR",
    noncoding = "Non-coding",
    "Transcript"
  )
}

browser_clean_label_vector <- function(value) {
  value <- trimws(as.character(value %||% ""))
  value[is.na(value)] <- ""
  gsub("[[:space:];]+$", "", value)
}

browser_compact_overlap_label_vector <- function(annotation_name, annotation_view) {
  region_label <- browser_annotation_region_label(annotation_view)
  annotation_name <- browser_clean_label_vector(annotation_name)

  first_part <- browser_clean_label_vector(sub("\\|.*$", "", annotation_name))
  if (identical(as.character(annotation_view), "noncoding")) {
    return(ifelse(nzchar(first_part), sprintf("%s|%s", first_part, region_label), region_label))
  }

  has_third_part <- grepl("^[^|]*\\|[^|]*\\|", annotation_name)
  third_part <- rep("", length(annotation_name))
  third_part[has_third_part] <- sub("^[^|]*\\|[^|]*\\|([^|]*).*$", "\\1", annotation_name[has_third_part])
  gene_id <- browser_clean_label_vector(ifelse(nzchar(third_part), third_part, first_part))

  ifelse(nzchar(gene_id), sprintf("%s|%s", gene_id, region_label), region_label)
}

browser_upload_ranges <- function(bed) {
  if (ncol(bed) < 3 || !nrow(bed)) {
    return(data.table::data.table())
  }

  bed_ranges <- data.table::data.table(
    chrom = browser_canonical_chr_vector(bed[[1]]),
    start = suppressWarnings(as.integer(bed[[2]])),
    end = suppressWarnings(as.integer(bed[[3]])),
    row_id = seq_len(nrow(bed))
  )

  bed_ranges[is.finite(start) & is.finite(end) & end > start]
}

browser_annotation_ranges_cache <- new.env(parent = emptyenv())

browser_annotation_cache_key <- function(annotation_path, annotation_view) {
  file_info <- file.info(annotation_path)
  paste(
    normalizePath(annotation_path, winslash = "/", mustWork = TRUE),
    as.character(annotation_view),
    as.integer(as.numeric(file_info$mtime)),
    file_info$size,
    sep = "::"
  )
}

browser_annotation_ranges <- function(annotation_path, annotation_view) {
  if (is.null(annotation_path) || !nzchar(annotation_path) || !file.exists(annotation_path)) {
    return(data.table::data.table())
  }

  cache_key <- browser_annotation_cache_key(annotation_path, annotation_view)
  if (exists(cache_key, envir = browser_annotation_ranges_cache, inherits = FALSE)) {
    return(get(cache_key, envir = browser_annotation_ranges_cache, inherits = FALSE))
  }

  annotation <- browser_read_bed_table(annotation_path)
  if (ncol(annotation) < 3 || !nrow(annotation)) {
    annotation_ranges <- data.table::data.table()
  } else {
    annotation_names <- if (ncol(annotation) >= 4) as.character(annotation[[4]]) else rep("", nrow(annotation))
    annotation_ranges <- data.table::data.table(
      chrom = browser_canonical_chr_vector(annotation[[1]]),
      start = suppressWarnings(as.integer(annotation[[2]])),
      end = suppressWarnings(as.integer(annotation[[3]])),
      annotation_transcript_id = sub("\\|.*$", "", annotation_names),
      overlap_label = browser_compact_overlap_label_vector(annotation_names, annotation_view)
    )
    annotation_ranges <- annotation_ranges[is.finite(start) & is.finite(end) & end > start]
    data.table::setkey(annotation_ranges, chrom, start, end)
  }

  assign(cache_key, annotation_ranges, envir = browser_annotation_ranges_cache)
  annotation_ranges
}

browser_empty_transcript_table <- function() {
  data.table::data.table(
    row_key = character(),
    source_file = character(),
    source_file_index = integer(),
    source_row = integer(),
    transcript_interval = character(),
    transcript_id = character(),
    gene_id = character(),
    transcript_length = integer(),
    igv_locus = character(),
    overlap_bases = integer()
  )
}

browser_transcript_locus <- function(transcript_id, transcript_length) {
  transcript_id <- as.character(transcript_id %||% "")[[1]]
  transcript_length <- suppressWarnings(as.integer(transcript_length))
  if (!nzchar(transcript_id) || !is.finite(transcript_length) || transcript_length < 1L) {
    return("")
  }

  sprintf("%s:1-%s", transcript_id, transcript_length)
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

browser_transcript_id_key_rows <- function(values, identifiers, identifier_column) {
  if (!exists("structure_strip_gff_prefix", mode = "function")) {
    result <- data.table::data.table(match_key = character())
    result[[identifier_column]] <- identifiers[FALSE]
    return(result)
  }

  values <- structure_strip_gff_prefix(values)
  keep <- nzchar(values)
  if (!any(keep)) {
    result <- data.table::data.table(match_key = character())
    result[[identifier_column]] <- identifiers[FALSE]
    return(result)
  }

  selected_values <- values[keep]
  selected_ids <- identifiers[keep]
  key_rows <- data.table::rbindlist(list(
    data.table::data.table(
      match_key = selected_values,
      identifier = selected_ids
    ),
    data.table::data.table(
      match_key = sub("\\.[0-9]+$", "", selected_values),
      identifier = selected_ids
    )
  ), use.names = TRUE)
  data.table::setnames(key_rows, "identifier", identifier_column)
  unique(key_rows)
}

browser_build_transcript_id_key_index <- function(transcript_exons) {
  empty_index <- data.table::data.table(
    match_key = character(),
    transcript_id = character()
  )
  if (
    !exists("structure_id_keys", mode = "function") ||
      is.null(transcript_exons) ||
      !nrow(transcript_exons) ||
      !"transcript_id" %in% names(transcript_exons)
  ) {
    return(empty_index)
  }

  metadata <- unique(transcript_exons[, intersect(
    c("transcript_id", "gene_id", "transcript_name"),
    names(transcript_exons)
  ), with = FALSE])
  if (!nrow(metadata)) {
    return(empty_index)
  }

  key_rows <- lapply(c("transcript_id", "gene_id", "transcript_name"), function(column) {
    values <- if (column %in% names(metadata)) metadata[[column]] else rep("", nrow(metadata))
    browser_transcript_id_key_rows(values, metadata$transcript_id, "transcript_id")
  })
  index <- data.table::rbindlist(key_rows, use.names = TRUE, fill = TRUE)
  if (!nrow(index)) {
    return(empty_index)
  }

  unique(index[, .(match_key, transcript_id)])
}

browser_prepare_transcript_exons <- function(transcript_spec, selected_transcript_id = NULL) {
  exons <- structure_load_transcript_map(transcript_spec)
  if (!is.null(selected_transcript_id)) {
    exons <- exons[transcript_id == as.character(selected_transcript_id)[[1L]]]
  }
  exons <- data.table::copy(exons)
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
  data.table::setkey(exons, chrom, start, end)
  attr(exons, "browser_transcript_id_key_index") <- browser_build_transcript_id_key_index(exons)
  exons[]
}

browser_map_bed_to_transcript_rows <- function(
  bed,
  file_entry,
  file_index,
  transcript_exons,
  transcript_id_key_index = NULL
) {
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

  overlap_exons <- transcript_exons
  if (!identical(data.table::key(overlap_exons), c("chrom", "start", "end"))) {
    overlap_exons <- data.table::copy(transcript_exons)
    data.table::setkey(overlap_exons, chrom, start, end)
  }
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
  mapped[, name_score := 0L]
  if (is.null(transcript_id_key_index)) {
    transcript_id_key_index <- attr(transcript_exons, "browser_transcript_id_key_index")
  }
  if (
    exists("structure_id_keys", mode = "function") &&
      !is.null(transcript_id_key_index) &&
      nrow(transcript_id_key_index)
  ) {
    bed_id_key_index <- browser_transcript_id_key_rows(
      mapped$bed_name,
      mapped$row_id,
      "row_id"
    )
    if (nrow(bed_id_key_index)) {
      name_hits <- unique(merge(
        bed_id_key_index,
        transcript_id_key_index,
        by = "match_key",
        allow.cartesian = TRUE
      )[, .(row_id, transcript_id)])
      if (nrow(name_hits)) {
        mapped[name_hits, on = .(row_id, transcript_id), name_score := 1L]
      }
    }
  }
  mapped[, strand_score := as.integer(bed_strand %in% c("+", "-") & bed_strand == strand)]
  mapped[, transcript_overlap_bases := as.integer(sum(overlap_bases, na.rm = TRUE)), by = .(row_id, transcript_id)]
  data.table::setorder(mapped, row_id, -name_score, -strand_score, -transcript_overlap_bases, transcript_id, transcript_interval_start, transcript_interval_end)
  mapped[, segment_index := seq_len(.N), by = .(row_id, transcript_id)]
  mapped[, transcript_interval := sprintf("%s-%s", transcript_interval_start, transcript_interval_end)]
  mapped[, igv_locus := mapply(browser_transcript_locus, transcript_id, transcript_length, USE.NAMES = FALSE)]
  mapped[, row_key := sprintf("f%s_r%s_tx%s_seg%s", file_index, row_id, transcript_id, segment_index)]
  mapped[, source_file := file_entry$file_name %||% basename(file_entry$file_path)]
  mapped[, source_file_index := as.integer(file_index)]
  mapped[, source_row := as.integer(row_id)]

  mapped <- merge(mapped, original_cols, by.x = "row_id", by.y = "row_id", all.x = TRUE, sort = FALSE)
  data.table::setorder(mapped, row_id, -name_score, -strand_score, -transcript_overlap_bases, transcript_id, segment_index)

  keep <- c(
    "row_key", "source_file", "source_file_index", "source_row", "transcript_interval", "transcript_id",
    "gene_id", "transcript_length", "igv_locus", "overlap_bases",
    paste0("bed_col_", seq_len(original_col_count))
  )
  mapped[, keep, with = FALSE]
}

browser_transcript_match_key <- function(transcript_id) {
  transcript_id <- trimws(as.character(transcript_id %||% ""))
  tolower(sub("\\.[0-9]+$", "", transcript_id))
}

browser_project_annotation_to_transcript <- function(annotation_ranges, transcript_exons, transcript_id) {
  empty <- data.table::data.table(start = integer(), end = integer())
  selected_id <- as.character(transcript_id %||% "")[[1]]
  if (!nzchar(selected_id) || !nrow(annotation_ranges) || !nrow(transcript_exons)) {
    return(empty)
  }

  exons <- data.table::copy(transcript_exons[transcript_id == selected_id])
  if (!nrow(exons)) {
    return(empty)
  }

  features <- data.table::copy(annotation_ranges[annotation_transcript_id == selected_id])
  if (!nrow(features)) {
    target_key <- browser_transcript_match_key(selected_id)
    features <- data.table::copy(annotation_ranges[
      browser_transcript_match_key(annotation_transcript_id) == target_key
    ])
  }
  if (!nrow(features)) {
    return(empty)
  }

  features <- features[, .(
    chrom,
    start = as.integer(start) + 1L,
    end = as.integer(end)
  )]
  features <- features[is.finite(start) & is.finite(end) & end >= start]
  if (!nrow(features)) {
    return(empty)
  }

  exons <- exons[, .(
    chrom,
    start = as.integer(start),
    end = as.integer(end),
    strand,
    cumulative_before = as.integer(cumulative_before)
  )]
  data.table::setkey(exons, chrom, start, end)
  data.table::setkey(features, chrom, start, end)
  hits <- data.table::foverlaps(
    features,
    exons,
    by.x = c("chrom", "start", "end"),
    by.y = c("chrom", "start", "end"),
    type = "any",
    nomatch = 0L
  )
  if (!nrow(hits)) {
    return(empty)
  }

  feature_start_col <- if ("i.start" %in% names(hits)) "i.start" else "start"
  feature_end_col <- if ("i.end" %in% names(hits)) "i.end" else "end"
  hits[, feature_start := get(feature_start_col)]
  hits[, feature_end := get(feature_end_col)]
  hits[, overlap_start := pmax(start, feature_start)]
  hits[, overlap_end := pmin(end, feature_end)]
  hits <- hits[overlap_end >= overlap_start]
  if (!nrow(hits)) {
    return(empty)
  }

  hits[, transcript_start := ifelse(
    strand == "-",
    cumulative_before + (end - overlap_end + 1L),
    cumulative_before + (overlap_start - start + 1L)
  )]
  hits[, transcript_end := ifelse(
    strand == "-",
    cumulative_before + (end - overlap_start + 1L),
    cumulative_before + (overlap_end - start + 1L)
  )]

  projected <- hits[, .(
    start = as.integer(pmin(transcript_start, transcript_end)),
    end = as.integer(pmax(transcript_start, transcript_end))
  )]
  projected <- unique(projected[start >= 1L & end >= start])
  data.table::setorder(projected, start, end)
  projected[]
}

browser_intersect_transcript_ranges <- function(ranges, segments) {
  empty <- data.table::data.table(start = integer(), end = integer())
  if (!nrow(ranges) || !nrow(segments)) {
    return(empty)
  }

  annotation_ranges <- data.table::copy(ranges)[, .(
    chrom = "transcript",
    start = as.integer(start),
    end = as.integer(end)
  )]
  bed_ranges <- data.table::copy(segments)[, .(
    chrom = "transcript",
    start = as.integer(start),
    end = as.integer(end)
  )]
  annotation_ranges <- annotation_ranges[is.finite(start) & is.finite(end) & end >= start]
  bed_ranges <- bed_ranges[is.finite(start) & is.finite(end) & end >= start]
  if (!nrow(annotation_ranges) || !nrow(bed_ranges)) {
    return(empty)
  }

  data.table::setkey(annotation_ranges, chrom, start, end)
  data.table::setkey(bed_ranges, chrom, start, end)
  hits <- data.table::foverlaps(
    annotation_ranges,
    bed_ranges,
    by.x = c("chrom", "start", "end"),
    by.y = c("chrom", "start", "end"),
    type = "any",
    nomatch = 0L
  )
  if (!nrow(hits)) {
    return(empty)
  }

  annotation_start_col <- if ("i.start" %in% names(hits)) "i.start" else "start"
  annotation_end_col <- if ("i.end" %in% names(hits)) "i.end" else "end"
  hits[, .(
    start = as.integer(pmax(start, get(annotation_start_col))),
    end = as.integer(pmin(end, get(annotation_end_col)))
  )][end >= start][, unique(.SD)][order(start, end)]
}
