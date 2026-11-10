motif_streme_header_info <- function(header_line, meta_line = "", fallback_index = 1L) {
  trimmed_header <- trimws(as.character(header_line))
  motif_label <- sub("^MOTIF\\s+", "", trimmed_header)
  motif_label <- sub("\\s+STREME-\\d+$", "", motif_label)
  motif_label <- sub("^[0-9]+-", "", motif_label)
  if (!nzchar(motif_label)) {
    motif_label <- sprintf("STREME-%d", fallback_index)
  }

  value_match <- regexec("(E|P|S)=\\s*([^[:space:]]+)", trimws(as.character(meta_line)))
  value_parts <- regmatches(trimws(as.character(meta_line)), value_match)[[1]]
  value_value <- if (length(value_parts) >= 3L) value_parts[[3]] else ""

  list(
    name = sprintf("Motif %d | %s", as.integer(fallback_index), motif_label),
    consensus = motif_label,
    score = NA_real_,
    evalue = if (length(value_parts) >= 3L && value_parts[[2]] == "E") value_value else NULL,
    pvalue = if (length(value_parts) >= 3L && value_parts[[2]] == "P") value_value else NULL,
    statistic = if (length(value_parts) >= 3L && value_parts[[2]] == "S") value_value else NULL
  )
}

motif_streme_pwm_list <- function(lines) {
  pwm_lines <- trimws(lines)
  pwm_lines <- pwm_lines[nzchar(pwm_lines)]

  pwm_table <- utils::read.table(
    text = paste(pwm_lines, collapse = "\n"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  pwm_matrix <- as.matrix(pwm_table)
  storage.mode(pwm_matrix) <- "numeric"

  if (ncol(pwm_matrix) >= 5L) {
    pwm_matrix <- pwm_matrix[, 2:5, drop = FALSE]
  }

  pwm_matrix <- motif_normalize_pwm_rows(pwm_matrix)

  list(
    A = as.numeric(pwm_matrix[, 1]),
    C = as.numeric(pwm_matrix[, 2]),
    G = as.numeric(pwm_matrix[, 3]),
    T = as.numeric(pwm_matrix[, 4])
  )
}

motif_parse_streme_motifs <- function(file_path, max_motifs = 6L) {
  if (!file.exists(file_path)) {
    return(list())
  }

  lines <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
  motif_header_index <- grep("^MOTIF\\s+", lines)

  if (!length(motif_header_index)) {
    return(list())
  }

  motif_limit <- max(1L, motif_safe_integer(max_motifs, 6L))

  motifs <- lapply(seq_len(min(length(motif_header_index), motif_limit)), function(index) {
    header_index <- motif_header_index[[index]]
    next_header_index <- if (index < length(motif_header_index)) motif_header_index[[index + 1L]] - 1L else length(lines)
    block_lines <- lines[header_index:next_header_index]

    meta_line_index <- grep("w=\\s*[0-9]+", block_lines)
    if (!length(meta_line_index)) {
      return(NULL)
    }

    meta_line <- block_lines[[meta_line_index[[1]]]]
    width_value <- suppressWarnings(as.numeric(sub(".*w=\\s*([0-9]+).*", "\\1", trimws(meta_line))))
    if (!is.finite(width_value) || is.na(width_value) || width_value <= 0) {
      return(NULL)
    }

    numeric_candidate_lines <- block_lines[(meta_line_index[[1]] + 1L):length(block_lines)]
    numeric_candidate_lines <- trimws(numeric_candidate_lines)
    numeric_candidate_lines <- numeric_candidate_lines[grepl("^[0-9.[:space:]]+$", numeric_candidate_lines)]
    if (!length(numeric_candidate_lines)) {
      return(NULL)
    }

    pwm_lines <- head(numeric_candidate_lines, width_value)
    if (length(pwm_lines) < width_value) {
      return(NULL)
    }

    header_info <- motif_streme_header_info(block_lines[[1]], meta_line, index)

    list(
      name = header_info$name,
      consensus = header_info$consensus,
      width = as.integer(width_value),
      score = header_info$score,
      pvalue = header_info$pvalue,
      evalue = header_info$evalue,
      statistic = header_info$statistic,
      pwm = motif_streme_pwm_list(pwm_lines)
    )
  })

  Filter(Negate(is.null), motifs)
}
