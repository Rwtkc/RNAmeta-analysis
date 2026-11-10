peak_exon_type_export_defaults <- function() {
  list(
    format = "png",
    width = 3000,
    height = 1800,
    dpi = 300
  )
}

peak_exon_type_export_config <- function(id) {
  ns <- NS(id)
  defaults <- peak_exon_type_export_defaults()

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
    choices = list(
      format = list(
        list(label = "PNG", value = "png"),
        list(label = "PDF", value = "pdf")
      ),
      dataFormat = list(
        list(label = "CSV", value = "csv"),
        list(label = "TXT (tab-delimited)", value = "txt")
      )
    ),
    sectioned = TRUE
  )
}

peak_exon_type_export_settings <- function(input) {
  defaults <- peak_exon_type_export_defaults()
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

peak_exon_type_export_filename <- function(extension) {
  sprintf("peak_exon_type_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

peak_exon_type_data_export_filename <- function(extension) {
  sprintf("peak_exon_type_plot_data_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

peak_exon_type_data_export_table <- function(result) {
  export_data <- data.table::copy(result$peak_exon_type_data)
  original_lookup <- stats::setNames(result$sample_name, result$sample_display_names)
  export_data[, Species := result$species]
  export_data[, Source_File := unname(original_lookup[as.character(Sample)])]
  data.table::setcolorder(export_data, c("Species", "Sample", "Source_File", "exonLoc", "Count", "Percentage"))
  as.data.frame(export_data)
}
