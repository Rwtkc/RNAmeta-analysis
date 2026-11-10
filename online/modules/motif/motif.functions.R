`%||%` <- function(value, fallback) {
  if (is.null(value)) {
    return(fallback)
  }

  value
}

motif_default_controls <- function() {
  list(
    engine = "meme",
    analysis_mode = "centered",
    sample_name = "",
    sites_subset = "all",
    motif_width = 5L,
    sequence_window = 9L,
    max_motifs = 6L,
    job_lookup = ""
  )
}

motif_analysis_mode_choices <- function() {
  c("Peak-centred window" = "centered", "Full peaks with matched background" = "matched_full")
}

motif_validate_analysis_mode <- function(value = NULL) {
  if (is.null(value)) return("centered")
  if (length(value) != 1L || is.na(value) || !value %in% unname(motif_analysis_mode_choices())) {
    stop("Select a supported motif analysis mode.")
  }
  as.character(value)
}

motif_validate_analysis_controls <- function(controls) {
  controls <- utils::modifyList(motif_default_controls(), controls)
  controls$analysis_mode <- motif_validate_analysis_mode(controls$analysis_mode)
  if (!controls$engine %in% unname(motif_engine_choices())) stop("Select a supported motif engine.")
  if (identical(controls$analysis_mode, "centered")) {
    controls$sequence_window <- motif_validate_sequence_window(controls$sequence_window,
      motif_effective_width(controls$engine, controls$motif_width))
  }
  controls
}

motif_engine_choices <- function() {
  c(
    "MEME" = "meme",
    "STREME" = "streme",
    "HOMER" = "homer"
  )
}

motif_sites_subset_choices <- function() {
  c(
    "All peaks" = "all",
    "Top 1000 by score" = "score_1000",
    "Top 2000 by score" = "score_2000",
    "Top 3000 by score" = "score_3000",
    "Top 1000 by peak length" = "length_1000",
    "Top 2000 by peak length" = "length_2000",
    "Top 3000 by peak length" = "length_3000"
  )
}

motif_engine_label <- function(engine) {
  labels <- motif_engine_choices()
  match_index <- match(as.character(engine)[[1]], unname(labels))
  if (!is.na(match_index)) {
    return(names(labels)[[match_index]])
  }

  toupper(as.character(engine)[[1]])
}

motif_safe_integer <- function(value, fallback = 0L) {
  numeric_value <- suppressWarnings(as.integer(value))
  if (!is.finite(numeric_value) || length(numeric_value) != 1L || is.na(numeric_value)) {
    return(as.integer(fallback))
  }

  as.integer(numeric_value)
}

motif_choice_options <- function(choices) {
  if (!length(choices)) {
    return(list())
  }

  choice_names <- names(choices)
  if (is.null(choice_names) || !length(choice_names)) {
    choice_names <- rep("", length(choices))
  }

  lapply(seq_along(choices), function(index) {
    value <- as.character(unname(choices[[index]]))
    label <- choice_names[[index]]
    if (is.null(label) || !nzchar(label)) {
      label <- value
    }

    list(label = label, value = value)
  })
}

motif_min_width <- function(engine) {
  engine_key <- tolower(as.character(engine)[[1]])

  if (identical(engine_key, "streme")) {
    return(3L)
  }

  if (identical(engine_key, "meme")) {
    return(2L)
  }

  1L
}

motif_effective_width <- function(engine, motif_width) {
  max(motif_min_width(engine), motif_safe_integer(motif_width, motif_default_controls()$motif_width))
}

motif_result_parse_failure_message <- function(engine_label, engine, motif_width) {
  engine_name <- as.character(engine_label)[[1]]
  engine_key <- tolower(as.character(engine)[[1]])
  width_value <- motif_safe_integer(motif_width, 12L)

  if (identical(engine_key, "homer") && width_value <= 1L) {
    return(sprintf(
      paste(
        "%s finished, but RNAmeta Online could not read a supported motif result from the generated output.",
        "This often happens when motif width is too small or the BED intervals are extremely short.",
        "Try a larger motif width, or upload a wider peak BED file before running motif analysis again."
      ),
      engine_name
    ))
  }

  sprintf(
    paste(
      "%s finished, but RNAmeta Online could not read a supported motif result from the generated output.",
      "Try a different motif width or a different BED file and run the job again."
    ),
    engine_name
  )
}

motif_normalize_subset_value <- function(subset_value = NULL) {
  raw_value <- subset_value
  if (is.null(raw_value)) {
    raw_value <- motif_default_controls()$sites_subset
  }
  raw_value <- as.character(raw_value)[[1]]

  if (!nzchar(raw_value) || identical(raw_value, "All")) {
    return("all")
  }

  if (raw_value %in% c("1000", "2000", "3000")) {
    return(sprintf("score_%s", raw_value))
  }

  valid_values <- unname(motif_sites_subset_choices())
  if (raw_value %in% valid_values) {
    return(raw_value)
  }

  motif_default_controls()$sites_subset
}

motif_subset_spec <- function(subset_value = NULL) {
  normalized_value <- motif_normalize_subset_value(subset_value)

  if (identical(normalized_value, "all")) {
    return(list(
      value = "all",
      mode = "all",
      size = Inf
    ))
  }

  parts <- strsplit(normalized_value, "_", fixed = TRUE)[[1]]
  list(
    value = normalized_value,
    mode = parts[[1]],
    size = motif_safe_integer(parts[[2]], 0L)
  )
}

motif_subset_label <- function(subset_value = NULL) {
  normalized_value <- motif_normalize_subset_value(subset_value)
  choices <- motif_sites_subset_choices()
  match_index <- match(normalized_value, unname(choices))

  if (!is.na(match_index)) {
    return(names(choices)[[match_index]])
  }

  names(choices)[[1]]
}

motif_build_sample_choices <- function(upload_context) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (!length(saved_files)) {
    return(character())
  }

  file_names <- vapply(saved_files, function(file) file$file_name, character(1))
  stats::setNames(file_names, file_names)
}

motif_selected_saved_file <- function(upload_context, sample_name = NULL) {
  saved_files <- upload_run_context_saved_files(upload_context)
  if (!length(saved_files)) {
    return(NULL)
  }

  if (is.null(sample_name) || !nzchar(sample_name)) {
    return(saved_files[[1]])
  }

  matched <- Filter(function(file) identical(file$file_name, sample_name), saved_files)
  if (length(matched)) {
    return(matched[[1]])
  }

  saved_files[[1]]
}

motif_subset_peak_ranges <- function(peak_gr, subset_value = "All") {
  if (is.null(peak_gr) || !length(peak_gr)) {
    return(peak_gr)
  }

  subset_spec <- motif_subset_spec(subset_value)
  if (identical(subset_spec$mode, "all")) {
    return(peak_gr)
  }

  subset_n <- motif_safe_integer(subset_spec$size, length(peak_gr))
  if (!is.finite(subset_n) || subset_n <= 0) {
    return(peak_gr)
  }

  ordered_peaks <- peak_gr
  if (identical(subset_spec$mode, "score") && "score" %in% names(S4Vectors::mcols(peak_gr))) {
    score_values <- suppressWarnings(as.numeric(S4Vectors::mcols(peak_gr)$score))
    if (any(!is.na(score_values))) {
      ordered_peaks <- peak_gr[order(score_values, decreasing = TRUE, na.last = TRUE)]
    }
  } else if (identical(subset_spec$mode, "length")) {
    peak_lengths <- as.numeric(GenomicRanges::width(peak_gr))
    ordered_peaks <- peak_gr[order(peak_lengths, decreasing = TRUE, na.last = TRUE)]
  }

  ordered_peaks[seq_len(min(length(ordered_peaks), subset_n))]
}

motif_read_uploaded_bundle <- function(bundle_input) {
  if (is.null(bundle_input) || !nrow(bundle_input)) {
    return(NULL)
  }

  bundle_name <- bundle_input$name[[1]]
  bundle_path <- bundle_input$datapath[[1]]
  bundle_extension <- tolower(tools::file_ext(bundle_name))

  if (!identical(bundle_extension, "json")) {
    return(list(
      bundleName = bundle_name,
      bundleStatus = "unsupported",
      bundleMessage = "Only JSON motif bundles are currently supported in local replay mode.",
      motifs = list()
    ))
  }

  bundle_data <- tryCatch(
    jsonlite::fromJSON(bundle_path, simplifyVector = FALSE),
    error = function(err) {
      list(
        bundleName = bundle_name,
        bundleStatus = "invalid",
        bundleMessage = conditionMessage(err),
        motifs = list()
      )
    }
  )

  if (!is.null(bundle_data$bundleStatus)) {
    return(bundle_data)
  }

  bundle_engine <- bundle_data$engine
  bundle_sample_name <- bundle_data$sampleName
  bundle_source_file <- bundle_data$sourceFile
  bundle_motifs <- bundle_data$motifs

  if (is.null(bundle_motifs)) {
    bundle_motifs <- list()
  }

  list(
    bundleName = bundle_name,
    bundleStatus = "loaded",
    bundleMessage = "Motif result bundle loaded for local replay.",
    engine = bundle_engine,
    sampleName = bundle_sample_name,
    sourceFile = bundle_source_file,
    motifs = bundle_motifs
  )
}
