source(app_path("services", "motif_streme_parser.R"), local = TRUE, encoding = "UTF-8")

motif_normalize_pwm_rows <- function(pwm_matrix) {
  if (is.null(pwm_matrix) || !nrow(pwm_matrix) || !ncol(pwm_matrix)) {
    return(pwm_matrix)
  }

  row_sums <- rowSums(pwm_matrix)
  row_sums[row_sums <= 0] <- 1
  pwm_matrix / row_sums
}

motif_homer_header_info <- function(header_line, fallback_index = 1L) {
  fields <- strsplit(as.character(header_line), "\t", fixed = TRUE)[[1]]
  consensus <- sub("^>", "", fields[[1]])
  motif_name <- consensus

  if (length(fields) >= 2L && nzchar(fields[[2]])) {
    motif_name <- sub(",.*$", "", fields[[2]])
  }

  score_value <- if (length(fields) >= 3L) suppressWarnings(as.numeric(fields[[3]])) else NA_real_
  pvalue_field <- if (length(fields) >= 6L) fields[[6]] else ""
  pvalue_match <- regexec("P:([^,]+)", pvalue_field)
  pvalue_parts <- regmatches(pvalue_field, pvalue_match)[[1]]
  pvalue_value <- if (length(pvalue_parts) >= 2L) pvalue_parts[[2]] else ""

  list(
    name = sprintf("Motif %d | %s", as.integer(fallback_index), motif_name),
    consensus = consensus,
    score = score_value,
    pvalue = pvalue_value
  )
}

motif_homer_pwm_list <- function(lines) {
  pwm_table <- utils::read.table(
    text = paste(lines, collapse = "\n"),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  pwm_matrix <- as.matrix(pwm_table)
  storage.mode(pwm_matrix) <- "numeric"
  pwm_matrix <- motif_normalize_pwm_rows(pwm_matrix)

  list(
    A = as.numeric(pwm_matrix[, 1]),
    C = as.numeric(pwm_matrix[, 2]),
    G = as.numeric(pwm_matrix[, 3]),
    T = as.numeric(pwm_matrix[, 4])
  )
}

motif_parse_homer_motifs <- function(file_path, max_motifs = 6L) {
  if (!file.exists(file_path)) {
    return(list())
  }

  lines <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[!grepl("^;", lines)]
  header_index <- which(substr(lines, 1L, 1L) == ">")

  if (!length(header_index)) {
    return(list())
  }

  start_index <- header_index + 1L
  end_index <- c(header_index[-1L] - 1L, length(lines))
  motif_limit <- max(1L, motif_safe_integer(max_motifs, 6L))

  motifs <- lapply(seq_len(min(length(header_index), motif_limit)), function(index) {
    header_info <- motif_homer_header_info(lines[[header_index[[index]]]], index)
    pwm_lines <- lines[start_index[[index]]:end_index[[index]]]
    pwm_lines <- pwm_lines[nzchar(pwm_lines)]

    list(
      name = header_info$name,
      consensus = header_info$consensus,
      width = length(pwm_lines),
      score = header_info$score,
      pvalue = header_info$pvalue,
      pwm = motif_homer_pwm_list(pwm_lines)
    )
  })

  Filter(Negate(is.null), motifs)
}

motif_parse_homer_result_dir <- function(output_dir, max_motifs = 6L) {
  candidate_files <- c(
    file.path(output_dir, "homerMotifs.all.motifs"),
    file.path(output_dir, "homerMotifs.motifs12")
  )

  for (candidate_file in candidate_files) {
    motifs <- motif_parse_homer_motifs(candidate_file, max_motifs = max_motifs)
    if (length(motifs)) {
      return(motifs)
    }
  }

  motif_dir <- file.path(output_dir, "homerResults")
  if (!dir.exists(motif_dir)) {
    return(list())
  }

  motif_files <- list.files(
    motif_dir,
    pattern = "^motif[0-9]+\\.motif$",
    full.names = TRUE
  )

  if (!length(motif_files)) {
    return(list())
  }

  motif_files <- motif_files[order(
    suppressWarnings(as.integer(sub("^motif([0-9]+)\\.motif$", "\\1", basename(motif_files))))
  )]
  motif_limit <- max(1L, motif_safe_integer(max_motifs, 6L))
  motif_files <- utils::head(motif_files, motif_limit)

  motifs <- lapply(seq_along(motif_files), function(index) {
    parsed <- motif_parse_homer_motifs(motif_files[[index]], max_motifs = 1L)
    if (!length(parsed)) {
      return(NULL)
    }

    motif_entry <- parsed[[1]]
    motif_entry$name <- sprintf("Motif %d | %s", as.integer(index), motif_entry$consensus)
    motif_entry
  })

  Filter(Negate(is.null), motifs)
}

motif_meme_header_info <- function(header_line, fallback_index = 1L) {
  trimmed_line <- trimws(as.character(header_line))
  consensus_match <- regexec("^MOTIF\\s+(.+?)\\s+MEME-(\\d+)", trimmed_line)
  consensus_parts <- regmatches(trimmed_line, consensus_match)[[1]]
  consensus <- if (length(consensus_parts) >= 2L) consensus_parts[[2]] else sprintf("MEME-%d", fallback_index)

  evalue_match <- regexec("E-value\\s*=\\s*([^[:space:]]+)", trimmed_line)
  evalue_parts <- regmatches(trimmed_line, evalue_match)[[1]]
  evalue_value <- if (length(evalue_parts) >= 2L) evalue_parts[[2]] else ""

  score_match <- regexec("llr\\s*=\\s*([^[:space:]]+)", trimmed_line)
  score_parts <- regmatches(trimmed_line, score_match)[[1]]
  score_value <- if (length(score_parts) >= 2L) suppressWarnings(as.numeric(score_parts[[2]])) else NA_real_

  list(
    name = sprintf("Motif %d | %s", as.integer(fallback_index), consensus),
    consensus = consensus,
    score = score_value,
    evalue = evalue_value
  )
}

motif_meme_pwm_list <- function(lines) {
  pwm_lines <- trimws(lines)
  pwm_lines <- pwm_lines[nzchar(pwm_lines)]

  pwm_table <- utils::read.table(
    text = paste(pwm_lines, collapse = "\n"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  pwm_matrix <- as.matrix(pwm_table)
  storage.mode(pwm_matrix) <- "numeric"
  pwm_matrix <- motif_normalize_pwm_rows(pwm_matrix)

  list(
    A = as.numeric(pwm_matrix[, 1]),
    C = as.numeric(pwm_matrix[, 2]),
    G = as.numeric(pwm_matrix[, 3]),
    T = as.numeric(pwm_matrix[, 4])
  )
}

motif_parse_meme_motifs <- function(file_path, max_motifs = 6L) {
  if (!file.exists(file_path)) {
    return(list())
  }

  lines <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
  motif_header_index <- grep("^MOTIF\\s+.+\\s+MEME-\\d+\\s+width\\s*=", lines)

  if (!length(motif_header_index)) {
    return(list())
  }

  motif_limit <- max(1L, motif_safe_integer(max_motifs, 6L))

  motifs <- lapply(seq_len(min(length(motif_header_index), motif_limit)), function(index) {
    header_index <- motif_header_index[[index]]
    next_header_index <- if (index < length(motif_header_index)) motif_header_index[[index + 1L]] - 1L else length(lines)
    block_lines <- lines[header_index:next_header_index]
    header_info <- motif_meme_header_info(lines[[header_index]], index)

    probability_title_index <- grep("position-specific probability matrix", block_lines, fixed = TRUE)
    if (!length(probability_title_index)) {
      return(NULL)
    }

    probability_block <- block_lines[probability_title_index[[1]]:length(block_lines)]
    matrix_header_index <- grep("^letter-probability matrix:", trimws(probability_block))
    if (!length(matrix_header_index)) {
      return(NULL)
    }

    matrix_lines <- probability_block[(matrix_header_index[[1]] + 1L):length(probability_block)]
    matrix_lines <- trimws(matrix_lines)
    matrix_end_index <- match(TRUE, grepl("^-{5,}$", matrix_lines))
    if (!is.na(matrix_end_index)) {
      matrix_lines <- matrix_lines[seq_len(matrix_end_index - 1L)]
    }
    matrix_lines <- matrix_lines[grepl("^[0-9.[:space:]]+$", matrix_lines)]

    if (!length(matrix_lines)) {
      return(NULL)
    }

    list(
      name = header_info$name,
      consensus = header_info$consensus,
      width = length(matrix_lines),
      score = header_info$score,
      evalue = header_info$evalue,
      pwm = motif_meme_pwm_list(matrix_lines)
    )
  })

  Filter(Negate(is.null), motifs)
}

motif_detect_local_result <- function(engine, job_package_dir, max_motifs = 6L) {
  if (is.null(job_package_dir) || !nzchar(job_package_dir)) {
    return(NULL)
  }

  engine_key <- tolower(as.character(engine)[[1]])
  output_dir <- file.path(job_package_dir, motif_engine_output_dir_name(engine_key))

  if (!dir.exists(output_dir)) {
    return(NULL)
  }

  if (identical(engine_key, "homer")) {
    motifs <- motif_parse_homer_result_dir(output_dir, max_motifs = max_motifs)

    if (!length(motifs)) {
      return(NULL)
    }

    return(list(
      result_name = basename(output_dir),
      result_status = "local-output",
      result_message = "Loaded HOMER motifs directly from the server output directory.",
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      motifs = motifs
    ))
  }

  if (identical(engine_key, "meme")) {
    motif_file <- file.path(output_dir, "meme.txt")
    motifs <- motif_parse_meme_motifs(motif_file, max_motifs = max_motifs)

    if (!length(motifs)) {
      return(NULL)
    }

    return(list(
      result_name = basename(output_dir),
      result_status = "local-output",
      result_message = "Loaded MEME motifs directly from the server output directory.",
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      motifs = motifs
    ))
  }

  if (!identical(engine_key, "streme")) {
    return(NULL)
  }

  motif_file <- file.path(output_dir, "streme.txt")
  motifs <- motif_parse_streme_motifs(motif_file, max_motifs = max_motifs)

  if (!length(motifs)) {
    return(NULL)
  }

  list(
    result_name = basename(output_dir),
    result_status = "local-output",
    result_message = "Loaded STREME motifs directly from the server output directory.",
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    motifs = motifs
  )
}
