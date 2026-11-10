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
  normalized <- tolower(sub("^chr", "", chr))
  result <- chr

  nuclear <- normalized %in% as.character(1:12)
  result[nuclear] <- normalized[nuclear]
  result[normalized %in% c("m", "mt", "mitochondria", "mitochondrion")] <- "Mt"
  result[normalized %in% c("c", "pt", "chloroplast", "plastid")] <- "Pt"

  result
}

browser_canonical_chr <- function(chromosome) {
  browser_canonical_chr_vector(chromosome %||% "")[[1]]
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

browser_compact_overlap_label <- function(annotation_name, annotation_view) {
  region_label <- browser_annotation_region_label(annotation_view)
  annotation_name <- browser_clean_label_vector(annotation_name)
  if (!length(annotation_name) || is.na(annotation_name[[1]])) {
    annotation_name <- ""
  }

  parts <- browser_clean_label_vector(strsplit(annotation_name[[1]], "\\|", fixed = FALSE)[[1]])
  parts <- parts[nzchar(parts)]

  if (identical(as.character(annotation_view), "noncoding")) {
    if (length(parts)) {
      return(sprintf("%s|%s", parts[[1]], region_label))
    }
    return(region_label)
  }

  gene_id <- if (length(parts) >= 3) parts[[3]] else if (length(parts) >= 1) parts[[1]] else ""
  if (nzchar(gene_id)) {
    return(sprintf("%s|%s", gene_id, region_label))
  }

  region_label
}

browser_select_filtered_upload_label <- function(overlap_labels, base_label, annotation_view, base_is_fallback = FALSE) {
  region_label <- browser_annotation_region_label(annotation_view)
  base_label <- browser_clean_label_vector(base_label)[[1]] %||% ""

  if (nzchar(base_label) && !isTRUE(base_is_fallback)) {
    return(sprintf("%s|%s", base_label, region_label))
  }

  overlap_labels <- browser_clean_label_vector(unlist(overlap_labels, use.names = FALSE))
  overlap_labels <- unique(overlap_labels[nzchar(overlap_labels)])
  if (length(overlap_labels)) {
    return(overlap_labels[[1]])
  }

  if (nzchar(base_label)) {
    return(sprintf("%s|%s", base_label, region_label))
  }

  region_label
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

browser_filtered_upload_path <- function(file_entry, annotation_path, annotation_view, cache_scope = NULL) {
  source_path <- file_entry$file_path
  if (
    is.null(source_path) || !nzchar(source_path) || !file.exists(source_path) ||
      is.null(annotation_path) || !nzchar(annotation_path) || !file.exists(annotation_path)
  ) {
    return(NULL)
  }

  filter_dir <- browser_public_subdir("filtered", cache_scope)

  source_md5 <- unname(tools::md5sum(source_path))
  annotation_md5 <- unname(tools::md5sum(annotation_path))
  source_mtime <- as.integer(as.numeric(file.info(source_path)$mtime))
  annotation_mtime <- as.integer(as.numeric(file.info(annotation_path)$mtime))
  target_name <- sprintf(
    "v7_%s_%s_%s_%s_%s_%s",
    gsub("[^A-Za-z0-9_.-]", "_", annotation_view),
    substr(gsub("[^A-Za-z0-9]", "", source_md5), 1, 12),
    substr(gsub("[^A-Za-z0-9]", "", annotation_md5), 1, 12),
    source_mtime,
    annotation_mtime,
    browser_sanitize_filename(file_entry$file_name)
  )
  file.path(filter_dir, target_name)
}

browser_filter_upload_by_annotation <- function(file_entry, annotation_path, annotation_view, bed = NULL, bed_ranges = NULL, annotation_ranges = NULL, cache_scope = NULL) {
  target_path <- browser_filtered_upload_path(file_entry, annotation_path, annotation_view, cache_scope = cache_scope)
  if (is.null(target_path)) {
    return(NULL)
  }

  if (file.exists(target_path)) {
    return(target_path)
  }

  source_path <- file_entry$file_path

  if (is.null(bed)) {
    bed <- browser_read_bed_table(source_path)
  }

  if (is.null(bed_ranges)) {
    bed_ranges <- browser_upload_ranges(bed)
  }

  if (is.null(annotation_ranges)) {
    annotation_ranges <- browser_annotation_ranges(annotation_path, annotation_view)
  }

  if (ncol(bed) < 3 || !nrow(bed) || !nrow(bed_ranges) || !nrow(annotation_ranges)) {
    file.create(target_path)
    return(target_path)
  }

  hits <- data.table::foverlaps(
    bed_ranges,
    annotation_ranges,
    by.x = c("chrom", "start", "end"),
    by.y = c("chrom", "start", "end"),
    type = "any",
    nomatch = 0L
  )
  keep_rows <- sort(unique(hits$row_id))

  bed_for_filter <- if ("row_id" %in% names(bed)) {
    data.table::copy(bed[, setdiff(names(bed), "row_id"), with = FALSE])
  } else {
    bed
  }
  filtered <- if (length(keep_rows)) data.table::copy(bed_for_filter[keep_rows]) else data.table::copy(bed_for_filter[0])

  if (length(keep_rows)) {
    name_col <- if (ncol(filtered) >= 4) names(filtered)[[4]] else "V4"
    if (!name_col %in% names(filtered)) {
      filtered[, (name_col) := ""]
    }

    overlap_labels <- hits[
      ,
      .(
        overlap_label = list(unique(overlap_label))
      ),
      by = row_id
    ]
    overlap_labels <- overlap_labels[match(keep_rows, row_id)]
    base_labels <- browser_clean_label_vector(filtered[[name_col]] %||% "")
    empty_labels <- !nzchar(base_labels)
    base_labels[empty_labels] <- sprintf(
      "%s:%s-%s",
      filtered[[1]][empty_labels],
      filtered[[2]][empty_labels],
      filtered[[3]][empty_labels]
    )
    filtered[[name_col]] <- mapply(
      browser_select_filtered_upload_label,
      overlap_labels$overlap_label,
      base_labels,
      MoreArgs = list(annotation_view = annotation_view),
      base_is_fallback = empty_labels,
      USE.NAMES = FALSE
    )
  }

  data.table::fwrite(filtered, file = target_path, sep = "\t", col.names = FALSE)
  target_path
}

browser_precompute_region_cache <- function(spec, saved_files, progress = NULL, cache_scope = NULL) {
  if (is.null(spec) || !length(saved_files)) {
    return(invisible(NULL))
  }

  views <- spec$annotation_views %||% list()
  region_view_ids <- setdiff(names(views), "transcript")
  total_steps <- max(1L, length(region_view_ids) * length(saved_files))
  completed <- 0L

  browser_progress_update(progress, 8, "Checking reference annotation tracks")

  upload_cache <- lapply(saved_files, function(file_entry) {
    bed <- browser_read_bed_table(file_entry$file_path)
    list(
      file_entry = file_entry,
      bed = bed,
      bed_ranges = browser_upload_ranges(bed)
    )
  })

  for (view_id in region_view_ids) {
    view <- browser_annotation_view_spec(spec, view_id)
    annotation_path <- browser_igv_asset_first(spec, view$files)
    if (is.null(annotation_path)) {
      next
    }

    annotation_ranges <- browser_annotation_ranges(annotation_path, view_id)

    for (upload_entry in upload_cache) {
      file_entry <- upload_entry$file_entry
      detail <- sprintf("Preparing %s BED track for %s", view$label, file_entry$file_name)
      browser_progress_update(progress, 10 + round(78 * completed / total_steps), detail)
      browser_filter_upload_by_annotation(
        file_entry = file_entry,
        annotation_path = annotation_path,
        annotation_view = view_id,
        bed = upload_entry$bed,
        bed_ranges = upload_entry$bed_ranges,
        annotation_ranges = annotation_ranges,
        cache_scope = cache_scope
      )
      completed <- completed + 1L
    }
  }

  browser_progress_update(progress, 92, "Writing IGV session inputs")
  invisible(NULL)
}

browser_public_url <- function(path) {
  relative <- sub(
    paste0("^", gsub("\\\\", "/", normalizePath(browser_public_dir(), winslash = "/", mustWork = TRUE)), "/?"),
    "",
    normalizePath(path, winslash = "/", mustWork = TRUE)
  )

  paste0("rnameta-igv/", relative)
}

browser_public_absolute_url <- function(url) {
  if (is.null(url) || !nzchar(url)) {
    return(url)
  }

  if (grepl("^(https?:)?//|^data:|^blob:|^/", url)) {
    return(url)
  }

  paste0("../", url)
}

browser_webapp_ready <- function() {
  dir.exists(app_path("external", "igv-webapp", "dist")) &&
    file.exists(app_path("external", "igv-webapp", "dist", "index.html"))
}

browser_initial_locus <- function(saved_files, default_locus) {
  if (!length(saved_files)) {
    return(default_locus)
  }

  first_file <- saved_files[[1]]$file_path
  if (is.null(first_file) || !file.exists(first_file)) {
    return(default_locus)
  }

  first_row <- tryCatch(
    data.table::fread(
      first_file,
      header = FALSE,
      sep = "\t",
      fill = TRUE,
      data.table = FALSE,
      quote = "",
      nrows = 1
    ),
    error = function(...) NULL
  )

  if (is.null(first_row) || !nrow(first_row) || ncol(first_row) < 3) {
    return(default_locus)
  }

  start <- suppressWarnings(as.integer(first_row[[2]][[1]]))
  end <- suppressWarnings(as.integer(first_row[[3]][[1]]))
  if (!is.finite(start) || !is.finite(end)) {
    return(default_locus)
  }

  center <- floor((start + end) / 2)
  left <- max(1L, center - 2500L)
  right <- max(left + 100L, center + 2500L)
  sprintf("%s:%d-%d", browser_canonical_chr(first_row[[1]][[1]]), left, right)
}
