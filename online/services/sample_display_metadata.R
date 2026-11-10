rnameta_sample_original_name <- function(file_info) {
  file_name <- if (!is.null(file_info$file_name) && nzchar(file_info$file_name)) {
    file_info$file_name
  } else if (!is.null(file_info$file_path) && nzchar(file_info$file_path)) {
    basename(file_info$file_path)
  } else {
    ""
  }

  sample_name <- tools::file_path_sans_ext(basename(file_name))
  if (nzchar(sample_name)) {
    return(sample_name)
  }

  "Group1"
}

rnameta_build_sample_metadata <- function(saved_files, alias_prefix = "Sample", long_name_limit = 18) {
  original_names <- vapply(saved_files, rnameta_sample_original_name, character(1))
  display_names <- sprintf("%s %d", alias_prefix, seq_along(original_names))

  lapply(seq_along(saved_files), function(index) {
    list(
      display_name = display_names[[index]],
      original_name = original_names[[index]]
    )
  })
}

rnameta_sample_display_names <- function(sample_metadata) {
  vapply(sample_metadata, function(entry) entry$display_name, character(1))
}

rnameta_sample_original_names <- function(sample_metadata) {
  vapply(sample_metadata, function(entry) entry$original_name, character(1))
}

rnameta_sample_original_lookup <- function(sample_metadata) {
  stats::setNames(
    rnameta_sample_original_names(sample_metadata),
    rnameta_sample_display_names(sample_metadata)
  )
}

rnameta_lookup_original_name <- function(sample_metadata, display_name, fallback = display_name) {
  original_lookup <- rnameta_sample_original_lookup(sample_metadata)
  original_name <- unname(original_lookup[[display_name]])

  if (is.null(original_name) || !nzchar(original_name)) {
    fallback
  } else {
    original_name
  }
}

rnameta_data_export_settings <- function(format_value) {
  format <- tolower(if (is.null(format_value) || identical(format_value, "")) "csv" else format_value)
  format <- if (format %in% c("csv", "txt")) format else "csv"

  list(
    format = format,
    extension = if (identical(format, "txt")) "txt" else "csv",
    delimiter = if (identical(format, "txt")) "\t" else ","
  )
}

rnameta_data_export_mime_type <- function(format) {
  if (identical(format, "txt")) {
    "text/plain;charset=utf-8"
  } else {
    "text/csv;charset=utf-8"
  }
}

rnameta_data_export_content <- function(export_data, delimiter = ",") {
  lines <- character()
  connection <- textConnection("lines", "w", local = TRUE)
  on.exit(close(connection), add = TRUE)

  utils::write.table(
    export_data,
    file = connection,
    sep = delimiter,
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE
  )

  paste(lines, collapse = "\n")
}
