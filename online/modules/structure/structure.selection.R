structure_parse_coordinate <- function(coordinate) {
  value <- trimws(as.character(coordinate %||% ""))
  match <- regexec("^([^:[:space:]]+):([0-9,]+)-([0-9,]+)$", value)
  parts <- regmatches(value, match)[[1]]

  if (length(parts) != 4L) {
    stop("Enter a coordinate like chr1:72338-74096 using the chromosome/start/end text from the original BED row.")
  }

  start <- suppressWarnings(as.integer(gsub(",", "", parts[[3]])))
  end <- suppressWarnings(as.integer(gsub(",", "", parts[[4]])))
  if (!is.finite(start) || !is.finite(end) || end <= start) {
    stop("Coordinate start/end must be numeric and end must be greater than start.")
  }

  list(
    raw = value,
    chrom = parts[[2]],
    canonical_chrom = browser_canonical_chr(parts[[2]]),
    start = start,
    end = end
  )
}

structure_match_coordinate_row <- function(file_entry, coordinate) {
  bed <- structure_bed_table(file_entry, limit = Inf)
  if (!nrow(bed)) {
    return(NULL)
  }

  match_row <- bed[
    chrom == coordinate$canonical_chrom &
      start == coordinate$start &
      end == coordinate$end
  ]

  if (!nrow(match_row)) {
    return(NULL)
  }

  match_row[1]
}

structure_resolve_coordinate_selection <- function(upload_context, sample_index, coordinate) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (!length(saved_files)) {
    stop("Save a BED file in Upload / Run before rendering Structure.")
  }

  sample_index <- suppressWarnings(as.integer(sample_index))
  if (!is.finite(sample_index)) {
    sample_index <- 1L
  }
  sample_index <- min(max(sample_index, 1L), length(saved_files))
  file_entry <- saved_files[[sample_index]]
  parsed_coordinate <- structure_parse_coordinate(coordinate)
  matched_row <- structure_match_coordinate_row(file_entry, parsed_coordinate)
  if (is.null(matched_row)) {
    stop("Coordinate was not found in the uploaded BED file. Enter a coordinate that exactly matches one BED row.")
  }

  list(
    file = file_entry,
    interval = matched_row,
    coordinate = parsed_coordinate,
    matched_bed_row = TRUE
  )
}

structure_mapping_row_value <- function(mapping_row, field, fallback = "") {
  value <- mapping_row[[field]] %||% fallback
  if (length(value) < 1L || is.na(value[[1]])) {
    return(fallback)
  }
  value[[1]]
}

structure_resolve_mapped_row_selection <- function(upload_context, mapping_row) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (!length(saved_files)) {
    stop("Save a BED file in Upload / Run before submitting RNAfold.")
  }

  file_index <- suppressWarnings(as.integer(structure_mapping_row_value(mapping_row, "source_file_index", NA_integer_)))
  if (!is.finite(file_index) || file_index < 1L || file_index > length(saved_files)) {
    stop("The selected mapping row no longer matches a saved BED file. Map the BED files again.")
  }

  chrom_raw <- trimws(as.character(structure_mapping_row_value(mapping_row, "bed_col_1")))
  start <- suppressWarnings(as.integer(structure_mapping_row_value(mapping_row, "bed_col_2", NA_integer_)))
  end <- suppressWarnings(as.integer(structure_mapping_row_value(mapping_row, "bed_col_3", NA_integer_)))
  chrom <- browser_canonical_chr(chrom_raw)
  if (!nzchar(chrom) || !is.finite(start) || !is.finite(end) || end <= start) {
    stop("The selected mapping row has invalid BED coordinates. Map the BED files again.")
  }

  interval <- data.table::data.table(
    row_index = as.integer(structure_mapping_row_value(mapping_row, "source_row", 1L)),
    chrom = chrom,
    start = start,
    end = end,
    name = as.character(structure_mapping_row_value(mapping_row, "bed_col_4", "")),
    strand = as.character(structure_mapping_row_value(mapping_row, "bed_col_6", "."))
  )

  transcript_interval <- trimws(as.character(structure_mapping_row_value(mapping_row, "transcript_interval", "")))
  interval_match <- regexec("^([0-9]+)-([0-9]+)$", transcript_interval)
  interval_parts <- regmatches(transcript_interval, interval_match)[[1]]
  if (length(interval_parts) != 3L) {
    stop("The selected row has no valid transcript interval. Map the BED files again.")
  }
  mapped_start <- suppressWarnings(as.integer(interval_parts[[2]]))
  mapped_end <- suppressWarnings(as.integer(interval_parts[[3]]))
  if (!is.finite(mapped_start) || !is.finite(mapped_end) || mapped_start < 1L || mapped_end < mapped_start) {
    stop("The selected row has no valid transcript interval. Map the BED files again.")
  }

  list(
    file = saved_files[[file_index]],
    interval = interval,
    coordinate = list(raw = sprintf("%s:%s-%s", chrom_raw, start, end), chrom = chrom_raw, canonical_chrom = chrom, start = start, end = end),
    matched_bed_row = TRUE,
    mapped_transcript_id = as.character(structure_mapping_row_value(mapping_row, "transcript_id", "")),
    mapped_gene_id = as.character(structure_mapping_row_value(mapping_row, "gene_id", "")),
    mapped_transcript_interval = list(start = mapped_start, end = mapped_end)
  )
}

structure_resolve_selected_interval <- function(upload_context, sample_index, interval_index) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (!length(saved_files)) {
    stop("Save a BED file in Upload / Run before rendering Structure.")
  }

  sample_index <- suppressWarnings(as.integer(sample_index))
  if (!is.finite(sample_index)) {
    sample_index <- 1L
  }
  sample_index <- min(max(sample_index, 1L), length(saved_files))
  file_entry <- saved_files[[sample_index]]

  interval_index <- suppressWarnings(as.integer(interval_index))
  if (!is.finite(interval_index)) {
    interval_index <- 1L
  }

  bed <- structure_bed_table(file_entry, limit = max(250L, interval_index))
  if (!nrow(bed)) {
    stop("The selected BED file has no valid intervals.")
  }

  match_row <- which(bed$row_index == interval_index)
  if (!length(match_row)) {
    interval_index <- min(max(interval_index, 1L), nrow(bed))
    match_row <- interval_index
  }

  list(
    file = file_entry,
    interval = bed[match_row[[1]]],
    interval_index = interval_index
  )
}
