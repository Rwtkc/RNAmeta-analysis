browser_build_transcript_table <- function(upload_context, cache_scope = NULL, progress = NULL) {
  species <- if (is.null(upload_context)) "" else upload_context$species
  saved_files <- if (is.null(upload_context)) list() else upload_run_context_saved_files(upload_context)
  transcript_spec <- browser_transcript_reference_spec(species)
  if (is.null(transcript_spec)) {
    stop("Transcript FASTA/GFF3 references are not ready for this species.")
  }

  if (!length(saved_files)) {
    return(browser_empty_transcript_table())
  }

  browser_progress_update(progress, 28, "Loading transcript annotation")
  transcript_exons <- browser_prepare_transcript_exons(transcript_spec)
  transcript_id_key_index <- attr(transcript_exons, "browser_transcript_id_key_index")
  browser_progress_update(progress, 42, "Reading saved BED intervals")

  tables <- lapply(seq_along(saved_files), function(file_index) {
    browser_progress_update(
      progress,
      48 + round(34 * (file_index - 1L) / length(saved_files)),
      sprintf("Mapping BED file %s of %s to transcript coordinates", file_index, length(saved_files))
    )
    bed <- browser_read_bed_table(saved_files[[file_index]]$file_path)
    mapped <- browser_map_bed_to_transcript_rows(
      bed = bed,
      file_entry = saved_files[[file_index]],
      file_index = file_index,
      transcript_exons = transcript_exons,
      transcript_id_key_index = transcript_id_key_index
    )
    browser_progress_update(
      progress,
      48 + round(34 * file_index / length(saved_files)),
      sprintf("Mapped BED file %s of %s", file_index, length(saved_files))
    )
    mapped
  })

  browser_progress_update(progress, 88, "Combining transcript-mapped intervals")
  table <- data.table::rbindlist(tables, fill = TRUE)
  if (!nrow(table)) {
    return(table)
  }

  table[, row_key := make.unique(row_key, sep = "_")]
  browser_progress_update(progress, 92, "Finalizing transcript mapping table")
  table[]
}

browser_transcript_table_page <- function(table, page = 1L, page_size = 5L) {
  if (is.null(table) || !nrow(table)) {
    return(table)
  }

  page <- max(1L, as.integer(page %||% 1L))
  page_size <- max(1L, as.integer(page_size %||% 5L))
  total_pages <- max(1L, ceiling(nrow(table) / page_size))
  page <- min(page, total_pages)
  start <- ((page - 1L) * page_size) + 1L
  end <- min(nrow(table), start + page_size - 1L)
  table[start:end]
}

browser_table_link_input_id <- function(row_key) {
  paste0("open_igv_", gsub("[^A-Za-z0-9_.-]+", "_", as.character(row_key %||% "")))
}

browser_bed_column_names <- function(n_cols) {
  standard_names <- c(
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts"
  )

  if (n_cols <= length(standard_names)) {
    return(standard_names[seq_len(n_cols)])
  }

  c(standard_names, paste0("extra_", seq_len(n_cols - length(standard_names))))
}

browser_table_page_buttons <- function(page, total_pages, ns = identity) {
  if (total_pages <= 1L) {
    return(NULL)
  }

  pages <- if (total_pages <= 8L) {
    seq_len(total_pages)
  } else {
    window_start <- max(1L, min(page - 4L, total_pages - 7L))
    window_end <- min(total_pages - 2L, window_start + 5L)
    unique(c(window_start:window_end, max(1L, total_pages - 1L):total_pages))
  }

  buttons <- list()
  previous_page <- 0L
  for (page_number in pages) {
    if (previous_page > 0L && page_number > previous_page + 1L) {
      buttons <- c(buttons, list(span(class = "browser-pager-ellipsis", "...")))
    }

    buttons <- c(buttons, list(tags$button(
      type = "button",
      class = if (identical(as.integer(page_number), as.integer(page))) {
        "browser-pager-page browser-pager-page--active"
      } else {
        "browser-pager-page"
      },
      onclick = sprintf(
        "Shiny.setInputValue('%s', %d, {priority: 'event'});",
        ns("table_page_jump"),
        page_number
      ),
      as.character(page_number)
    )))
    previous_page <- page_number
  }

  buttons
}

browser_render_transcript_table <- function(table, page = 1L, page_size = 5L, ns = identity, selected_key = NULL, search_value = "") {
  if (is.null(table) || !nrow(table)) {
    return(div(
      class = "browser-transcript-table-placeholder",
      tags$strong("No transcript overlap was found."),
      tags$p("Check that the saved BED coordinates use the same reference as the selected species.")
    ))
  }

  total_pages <- max(1L, ceiling(nrow(table) / page_size))
  page <- min(max(1L, as.integer(page %||% 1L)), total_pages)
  page_rows <- browser_transcript_table_page(table, page = page, page_size = page_size)
  bed_cols <- grep("^bed_col_[0-9]+$", names(table), value = TRUE)
  bed_order <- order(as.integer(sub("^bed_col_", "", bed_cols)))
  bed_cols <- bed_cols[bed_order]
  first_bed_cols <- bed_cols[seq_len(min(3L, length(bed_cols)))]
  remaining_bed_cols <- if (length(bed_cols) > 3L) bed_cols[4:length(bed_cols)] else character()
  display_cols <- c(first_bed_cols, "transcript_interval", "transcript_id", "igv_link", remaining_bed_cols)
  bed_headers <- browser_bed_column_names(length(bed_cols))
  headers <- c(
    bed_headers[seq_along(first_bed_cols)],
    "Transcript interval",
    "Transcript ID",
    "IGV link",
    if (length(remaining_bed_cols)) bed_headers[seq(4L, length(bed_cols))] else character()
  )

  rows <- lapply(seq_len(nrow(page_rows)), function(row_index) {
    row <- page_rows[row_index]
    selected <- identical(as.character(row$row_key[[1]]), as.character(selected_key %||% ""))
    cells <- lapply(display_cols, function(col_name) {
      if (identical(col_name, "igv_link")) {
        if (nzchar(row$transcript_id[[1]] %||% "") && nzchar(row$igv_locus[[1]] %||% "")) {
          return(tags$td(
            tags$a(
              href = "#",
              class = "browser-transcript-table__link",
              `data-row-key` = row$row_key[[1]],
              onclick = sprintf(
                "event.preventDefault(); event.stopPropagation(); Shiny.setInputValue('%s', '%s', {priority: 'event'}); return false;",
                ns("browser_table_click"),
                gsub("'", "\\\\'", row$row_key[[1]], fixed = TRUE)
              ),
              "Open in IGV"
            )
          ))
        }
        return(tags$td(class = "browser-transcript-table__muted", "No transcript"))
      }

      value <- row[[col_name]][[1]] %||% ""
      tags$td(htmltools::htmlEscape(as.character(value)))
    })

    tags$tr(
      class = if (selected) "browser-transcript-table__row browser-transcript-table__row--selected" else "browser-transcript-table__row",
      cells
    )
  })

  div(
    class = "browser-transcript-table-block",
    div(
      class = "browser-transcript-table__toolbar",
      div(class = "browser-panel__eyebrow", "Transcript Mapping Table"),
      div(
        class = "browser-transcript-table__search",
        textInput(
          inputId = ns("table_search"),
          label = NULL,
          value = search_value %||% "",
          placeholder = "Search chrom, transcript ID, BED columns..."
        )
      )
    ),
    div(
      class = "browser-transcript-table-wrap",
      tags$table(
        class = "browser-transcript-table",
        tags$thead(tags$tr(lapply(headers, tags$th))),
        tags$tbody(rows)
      )
    ),
    div(
      class = "browser-transcript-table__footer",
      span(
        class = "browser-transcript-table__count",
        sprintf("Showing %s-%s of %s transcript-mapped rows", ((page - 1L) * page_size) + 1L, ((page - 1L) * page_size) + nrow(page_rows), nrow(table))
      ),
      div(
        class = "browser-transcript-table__pager browser-transcript-table__pager--bottom",
        actionButton(ns("table_prev_bottom"), "Previous", class = "browser-pager-btn", disabled = if (page <= 1L) "disabled" else NULL),
        browser_table_page_buttons(page, total_pages, ns = ns),
        numericInput(
          ns("table_page_jump"),
          label = "Go to page",
          value = page,
          min = 1L,
          max = total_pages,
          step = 1L,
          width = "7rem"
        ),
        actionButton(ns("table_next_bottom"), "Next", class = "browser-pager-btn", disabled = if (page >= total_pages) "disabled" else NULL)
      )
    )
  )
}

browser_write_lines_if_changed <- function(lines, path) {
  lines <- if (is.null(lines)) character() else as.character(lines)
  existing <- if (file.exists(path)) readLines(path, warn = FALSE) else NULL

  if (!identical(existing, lines)) {
    writeLines(lines, path, useBytes = TRUE)
  }

  path
}

browser_write_table_if_changed <- function(table, path) {
  lines <- if (!is.null(table) && nrow(table)) {
    do.call(paste, c(unname(lapply(table, as.character)), sep = "\t"))
  } else {
    character()
  }

  browser_write_lines_if_changed(lines, path)
}

browser_selected_transcript_workspace <- function(selected_transcript, cache_scope = NULL) {
  if (is.null(selected_transcript) || !length(selected_transcript)) {
    return(NULL)
  }

  transcript_id <- as.character(selected_transcript$transcript_id %||% "")[[1]]
  transcript_length <- suppressWarnings(as.integer(selected_transcript$transcript_length %||% NA_integer_))
  segments <- browser_parse_transcript_segments(selected_transcript$transcript_interval %||% "")

  if (!nzchar(transcript_id) || !is.finite(transcript_length) || transcript_length < 1L || !nrow(segments)) {
    return(NULL)
  }

  segments[, start := pmax(1L, pmin(transcript_length, start))]
  segments[, end := pmax(start, pmin(transcript_length, end))]
  segments <- segments[is.finite(start) & is.finite(end) & end >= start]
  if (!nrow(segments)) {
    return(NULL)
  }

  interval_start <- min(segments$start)
  interval_end <- max(segments$end)

  target_dir <- browser_public_subdir("selected_transcripts", cache_scope)
  target_base <- gsub("[^A-Za-z0-9_.-]+", "_", paste(
    selected_transcript$row_key %||% transcript_id,
    transcript_id,
    sep = "_"
  ))
  sizes_path <- file.path(target_dir, sprintf("%s.chrom.sizes", target_base))
  bed_path <- file.path(target_dir, sprintf("%s.bed", target_base))

  browser_write_lines_if_changed(
    sprintf("%s\t%s", transcript_id, transcript_length),
    sizes_path
  )

  segments[, segment_index := seq_len(.N)]
  bed <- data.table::data.table(
    chrom = transcript_id,
    start = segments$start - 1L,
    end = segments$end,
    name = sprintf("%s:%s-%s", transcript_id, segments$start, segments$end),
    score = 0L,
    strand = "+"
  )
  browser_write_table_if_changed(bed, bed_path)

  list(
    transcript_id = transcript_id,
    transcript_length = transcript_length,
    chrom_sizes_path = sizes_path,
    bed_path = bed_path,
    segments = segments,
    locus = browser_transcript_locus(transcript_id, transcript_length)
  )
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

browser_selected_transcript_cytoband <- function(workspace, spec) {
  if (is.null(workspace) || is.null(spec)) {
    return(NULL)
  }

  transcript_spec <- browser_transcript_reference_spec(spec$id)
  if (is.null(transcript_spec)) {
    return(NULL)
  }

  transcript_exons <- browser_prepare_transcript_exons(transcript_spec, selected_transcript_id = workspace$transcript_id)
  if (!nrow(transcript_exons)) {
    return(NULL)
  }

  cds_view <- browser_annotation_view_spec(spec, "cds")
  cds_path <- browser_igv_asset_first(spec, cds_view$files)
  cds_ranges <- browser_project_annotation_to_transcript(
    annotation_ranges = browser_annotation_ranges(cds_path, "cds"),
    transcript_exons = transcript_exons,
    transcript_id = workspace$transcript_id
  )
  target_dir <- dirname(workspace$bed_path)
  target_base <- sub("\\.bed$", "", basename(workspace$bed_path))
  cytoband_path <- file.path(target_dir, sprintf("%s_cytoband.txt", target_base))
  band_rows <- data.table::data.table(
    chrom = character(),
    start = integer(),
    end = integer(),
    name = character(),
    stain = character()
  )

  if (nrow(cds_ranges)) {
    cds_start <- max(1L, min(workspace$transcript_length, min(cds_ranges$start)))
    cds_end <- min(workspace$transcript_length, max(cds_ranges$end))
    band_rows <- data.table::rbindlist(list(
      if (cds_start > 1L) data.table::data.table(
        chrom = workspace$transcript_id,
        start = 0L,
        end = cds_start - 1L,
        name = "5'UTR",
        stain = "gpos50"
      ) else NULL,
      data.table::data.table(
        chrom = workspace$transcript_id,
        start = cds_start - 1L,
        end = cds_end,
        name = "CDS",
        stain = "gneg"
      ),
      if (cds_end < workspace$transcript_length) data.table::data.table(
        chrom = workspace$transcript_id,
        start = cds_end,
        end = workspace$transcript_length,
        name = "3'UTR",
        stain = "gpos50"
      ) else NULL
    ), fill = TRUE)
  } else {
    band_rows <- data.table::data.table(
      chrom = workspace$transcript_id,
      start = 0L,
      end = as.integer(workspace$transcript_length),
      name = "Transcript",
      stain = "gpos50"
    )
  }

  browser_write_table_if_changed(band_rows, cytoband_path)
  list(
    path = cytoband_path,
    segments = lapply(seq_len(nrow(band_rows)), function(index) {
      list(
        label = band_rows$name[[index]],
        start = as.integer(band_rows$start[[index]]),
        end = as.integer(band_rows$end[[index]])
      )
    })
  )
}

browser_selected_transcript_fasta <- function(workspace, spec) {
  if (is.null(workspace) || is.null(spec)) {
    return(NULL)
  }

  transcript_spec <- browser_transcript_reference_spec(spec$id)
  if (is.null(transcript_spec)) {
    return(NULL)
  }

  target_base <- sub("\\.bed$", "", basename(workspace$bed_path))
  target_dir <- dirname(workspace$bed_path)
  fasta_path <- file.path(target_dir, sprintf("%s.fa", target_base))
  index_path <- paste0(fasta_path, ".fai")
  if (file.exists(fasta_path) && file.exists(index_path)) {
    return(list(path = fasta_path, index_path = index_path))
  }

  transcript_exons <- browser_prepare_transcript_exons(transcript_spec, selected_transcript_id = workspace$transcript_id)
  transcript_name <- if (nrow(transcript_exons) && "transcript_name" %in% names(transcript_exons)) {
    as.character(transcript_exons$transcript_name[[1]])
  } else {
    workspace$transcript_id
  }
  source_fasta <- app_path("data", "fa", transcript_spec$cdna)
  fasta_seqname <- structure_match_transcript_seqname(
    transcript_id = workspace$transcript_id,
    transcript_name = transcript_name,
    fasta_path = source_fasta
  )
  fasta_index <- Rsamtools::scanFaIndex(source_fasta)
  seq_names <- as.character(GenomeInfoDb::seqnames(fasta_index))
  sequence_length <- as.integer(IRanges::width(fasta_index)[match(fasta_seqname, seq_names)])
  sequence_range <- GenomicRanges::GRanges(
    seqnames = fasta_seqname,
    ranges = IRanges::IRanges(start = 1L, end = sequence_length)
  )
  sequence <- toupper(as.character(Rsamtools::scanFa(
    Rsamtools::FaFile(source_fasta),
    param = sequence_range
  )[[1]]))
  if (!nzchar(sequence)) {
    stop(sprintf("No cDNA sequence was found for transcript '%s'.", workspace$transcript_id))
  }

  line_width <- 60L
  starts <- seq.int(1L, nchar(sequence), by = line_width)
  sequence_lines <- substring(sequence, starts, pmin(starts + line_width - 1L, nchar(sequence)))
  header <- paste0(">", workspace$transcript_id)
  browser_write_lines_if_changed(c(header, sequence_lines), fasta_path)
  browser_write_lines_if_changed(
    paste(
      workspace$transcript_id,
      nchar(sequence),
      nchar(header, type = "bytes") + 1L,
      min(line_width, nchar(sequence)),
      min(line_width, nchar(sequence)) + 1L,
      sep = "\t"
    ),
    index_path
  )

  list(path = fasta_path, index_path = index_path)
}

browser_build_selected_transcript_tracks <- function(selected_transcript, spec = NULL, cache_scope = NULL) {
  workspace <- browser_selected_transcript_workspace(selected_transcript, cache_scope = cache_scope)
  if (is.null(workspace)) {
    return(NULL)
  }

  cytoband <- browser_selected_transcript_cytoband(
    workspace = workspace,
    spec = spec
  )
  fasta_files <- browser_selected_transcript_fasta(
    workspace = workspace,
    spec = spec
  )
  interval_track <- list(
    name = "BED interval",
    type = "annotation",
    format = "bed",
    url = browser_public_url(workspace$bed_path),
    displayMode = "EXPANDED",
    margin = 30L,
    height = max(80L, min(180L, 46L + (nrow(workspace$segments) * 18L))),
    color = "#5d7f6b",
    sourcePath = workspace$bed_path
  )
  list(
    tracks = list(interval_track),
    transcriptReference = c(workspace, list(
      cytoband_path = cytoband$path %||% NULL,
      transcript_segments = cytoband$segments %||% list(),
      fasta_path = fasta_files$path %||% NULL,
      fasta_index_path = fasta_files$index_path %||% NULL
    ))
  )
}

browser_uploaded_tracks <- function(saved_files, annotation_path = NULL, annotation_view = "transcript", annotation_label = "Transcript", use_filtered_cache = FALSE, cache_scope = NULL) {
  if (!length(saved_files)) {
    return(list())
  }

  colors <- c("#2f7fbf", "#2f9f73", "#b78f20", "#2447b8", "#c86d52", "#7e5aa8")

  tracks <- lapply(seq_along(saved_files), function(index) {
    if (identical(annotation_view, "transcript")) {
      public_path <- browser_public_upload_path(saved_files[[index]], cache_scope = cache_scope)
      track_name <- saved_files[[index]]$file_name
    } else {
      public_path <- browser_filtered_upload_path(
        saved_files[[index]],
        annotation_path = annotation_path,
        annotation_view = annotation_view,
        cache_scope = cache_scope
      )
      if (!isTRUE(use_filtered_cache) || is.null(public_path) || !file.exists(public_path)) {
        public_path <- browser_public_upload_path(saved_files[[index]], cache_scope = cache_scope)
      }
      track_name <- sprintf("%s (%s)", saved_files[[index]]$file_name, annotation_label)
    }

    if (is.null(public_path)) {
      return(NULL)
    }

    list(
      name = saved_files[[index]]$file_name,
      type = "annotation",
      format = "bed",
      url = browser_public_url(public_path),
      displayMode = "EXPANDED",
      height = 80,
      color = colors[((index - 1) %% length(colors)) + 1]
    )
  })

  Filter(Negate(is.null), tracks)
}
