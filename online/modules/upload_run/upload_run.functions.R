upload_run_detect_bed_header <- function(data_frame) {
  if (is.null(data_frame) || !nrow(data_frame) || ncol(data_frame) < 3) {
    return(FALSE)
  }

  first_values <- vapply(
    data_frame[1, seq_len(3), drop = TRUE],
    function(value) {
      gsub("[^a-z0-9]", "", tolower(trimws(as.character(value))))
    },
    character(1)
  )

  first_values[[1]] %in% c("chrom", "chr", "chromosome", "seqname", "seqnames") &&
    first_values[[2]] %in% c("chromstart", "start", "txstart") &&
    first_values[[3]] %in% c("chromend", "end", "txend")
}

upload_run_is_numeric_like <- function(value) {
  !is.na(suppressWarnings(as.numeric(trimws(as.character(value)))))
}

upload_run_looks_like_chromosome <- function(value) {
  normalized <- tolower(trimws(as.character(value)))
  grepl("^(chr)?([0-9]+|[xyzm]|mt|chloroplast|pt|c|scaffold|contig)", normalized)
}

upload_run_detect_bed_header_by_shape <- function(data_frame) {
  if (is.null(data_frame) || nrow(data_frame) < 2 || ncol(data_frame) < 3) {
    return(FALSE)
  }

  first_row <- data_frame[1, seq_len(3), drop = TRUE]
  next_rows <- data_frame[seq_len(min(nrow(data_frame), 5))[-1], seq_len(3), drop = FALSE]

  if (!nrow(next_rows)) {
    return(FALSE)
  }

  first_row_non_numeric <- !upload_run_is_numeric_like(first_row[[2]]) &&
    !upload_run_is_numeric_like(first_row[[3]])

  subsequent_numeric <- all(vapply(next_rows[[2]], upload_run_is_numeric_like, logical(1))) &&
    all(vapply(next_rows[[3]], upload_run_is_numeric_like, logical(1)))

  subsequent_chr_like <- sum(vapply(next_rows[[1]], upload_run_looks_like_chromosome, logical(1))) >= 1

  first_row_non_numeric && subsequent_numeric && subsequent_chr_like
}

upload_run_normalize_bed_file <- function(file_path) {
  if (is.null(file_path) || !nzchar(file_path) || !file.exists(file_path)) {
    stop("Uploaded BED file was not found.")
  }

  # Most uploaded BED files are already headerless. Avoid a full read/write cycle
  # on every Save unless the small initial sample confirms that a header exists.
  initial_rows <- data.table::fread(
    file_path,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    data.table = FALSE,
    quote = "",
    nrows = 5L
  )

  has_header <- upload_run_detect_bed_header(initial_rows) ||
    upload_run_detect_bed_header_by_shape(initial_rows)
  if (!has_header) {
    return(file_path)
  }

  bed_tbl <- data.table::fread(
    file_path,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    data.table = FALSE,
    quote = ""
  )
  bed_tbl <- bed_tbl[-1, , drop = FALSE]

  normalized_path <- tempfile(pattern = "upload_run_", fileext = ".bed")
  data.table::fwrite(bed_tbl, file = normalized_path, sep = "\t", col.names = FALSE)
  normalized_path
}

upload_run_saved_files <- function(uploaded_files) {
  if (is.null(uploaded_files) || !nrow(uploaded_files)) {
    return(list())
  }

  lapply(seq_len(nrow(uploaded_files)), function(index) {
    list(
      file_name = uploaded_files$name[[index]],
      file_path = upload_run_normalize_bed_file(uploaded_files$datapath[[index]]),
      source = "Upload"
    )
  })
}

upload_run_build_saved_context <- function(data_source, species, saved_files, preview = NULL, is_demo = FALSE) {
  file_count <- length(saved_files)
  first_file <- if (file_count) saved_files[[1]] else list(file_name = "", file_path = "")

  list(
    data_source = data_source,
    species = species,
    file_count = file_count,
    file_name = first_file$file_name,
    file_path = first_file$file_path,
    file_names = vapply(saved_files, function(file) file$file_name, character(1)),
    saved_files = saved_files,
    preview = preview,
    is_demo = is_demo
  )
}

upload_run_saved_file_labels <- function(saved_state) {
  if (is.null(saved_state) || is.null(saved_state$saved_files) || !length(saved_state$saved_files)) {
    return(character())
  }

  vapply(saved_state$saved_files, function(file) file$file_name, character(1))
}

upload_run_context_saved_files <- function(upload_context) {
  if (!is.null(upload_context$saved_files) && length(upload_context$saved_files)) {
    return(upload_context$saved_files)
  }

  if (!is.null(upload_context$file_path) && nzchar(upload_context$file_path)) {
    return(list(
      list(
        file_name = if (is.null(upload_context$file_name)) basename(upload_context$file_path) else upload_context$file_name,
        file_path = upload_context$file_path
      )
    ))
  }

  list()
}

upload_run_context_file_count <- function(upload_context) {
  length(upload_run_context_saved_files(upload_context))
}

upload_run_selected_file <- function(saved_state, index = 1) {
  if (is.null(saved_state) || is.null(saved_state$saved_files) || !length(saved_state$saved_files)) {
    return(NULL)
  }

  safe_index <- suppressWarnings(as.integer(index))
  if (!is.finite(safe_index)) {
    safe_index <- 1L
  }
  safe_index <- min(max(safe_index, 1L), length(saved_state$saved_files))
  saved_state$saved_files[[safe_index]]
}
