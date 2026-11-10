peak_distribution_export_defaults <- function() {
  list(
    format = "png",
    width = 3000,
    height = 1800,
    dpi = 300
  )
}

peak_distribution_export_config <- function(id) {
  ns <- NS(id)
  defaults <- peak_distribution_export_defaults()

  list(
    ids = list(
      format = ns("export_format"),
      width = ns("export_width"),
      height = ns("export_height"),
      dpi = ns("export_dpi"),
      trigger = ns("export_plot"),
      dataFormat = ns("data_export_format"),
      dataTrigger = ns("data_export")
    ),
    defaults = defaults,
    sectioned = TRUE,
    choices = list(
      format = list(
        list(label = "PNG", value = "png"),
        list(label = "PDF", value = "pdf")
      ),
      dataFormat = list(
        list(label = "CSV", value = "csv"),
        list(label = "TXT (tab-delimited)", value = "txt")
      )
    )
  )
}

peak_distribution_export_settings <- function(input) {
  defaults <- peak_distribution_export_defaults()

  width <- suppressWarnings(as.numeric(input$export_width))
  height <- suppressWarnings(as.numeric(input$export_height))
  dpi <- suppressWarnings(as.numeric(input$export_dpi))
  format <- tolower(if (is.null(input$export_format) || identical(input$export_format, "")) defaults$format else input$export_format)

  list(
    format = if (format %in% c("png", "pdf")) format else defaults$format,
    width = if (is.finite(width) && width > 0) width else defaults$width,
    height = if (is.finite(height) && height > 0) height else defaults$height,
    dpi = if (is.finite(dpi) && dpi > 0) dpi else defaults$dpi
  )
}

peak_distribution_export_filename <- function(extension) {
  sprintf("peak_distribution_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

peak_distribution_data_export_defaults <- function() {
  list(
    format = "csv"
  )
}

peak_distribution_data_export_settings <- function(input) {
  defaults <- peak_distribution_data_export_defaults()
  format <- tolower(if (is.null(input$data_export_format) || identical(input$data_export_format, "")) defaults$format else input$data_export_format)
  format <- if (format %in% c("csv", "txt")) format else defaults$format

  list(
    format = format,
    extension = if (identical(format, "txt")) "txt" else "csv",
    delimiter = if (identical(format, "txt")) "\t" else ","
  )
}

peak_distribution_data_export_filename <- function(extension) {
  sprintf("peak_distribution_plot_data_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

peak_distribution_data_export_mime_type <- function(format) {
  if (identical(format, "txt")) {
    "text/plain;charset=utf-8"
  } else {
    "text/csv;charset=utf-8"
  }
}

peak_distribution_data_export_table <- function(result) {
  export_data <- result$feature_counts
  if (is.null(export_data) || !nrow(export_data)) {
    return(data.frame())
  }

  sample_display_names <- result$sample_display_names
  sample_original_names <- result$sample_name
  name_lookup <- stats::setNames(sample_original_names, sample_display_names)
  source_files <- unname(name_lookup[as.character(export_data$Sample)])
  source_files[is.na(source_files)] <- as.character(export_data$Sample[is.na(source_files)])
  species <- if (is.null(result$species)) "" else result$species

  data.frame(
    Species = rep(species, nrow(export_data)),
    Sample = as.character(export_data$Sample),
    Source_File = source_files,
    Feature = as.character(export_data$Feature),
    Frequency = as.numeric(export_data$Frequency),
    stringsAsFactors = FALSE
  )
}

peak_distribution_data_export_content <- function(result, delimiter = ",") {
  export_data <- peak_distribution_data_export_table(result)
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
