meta_plot_export_defaults <- function() {
  list(
    format = "png",
    width = 3000,
    height = 1800,
    dpi = 300
  )
}

meta_plot_export_config <- function(id) {
  ns <- NS(id)
  defaults <- meta_plot_export_defaults()

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

meta_plot_export_settings <- function(input) {
  defaults <- meta_plot_export_defaults()

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

meta_plot_export_filename <- function(extension) {
  sprintf("meta_plot_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

meta_plot_data_export_filename <- function(extension) {
  sprintf("meta_plot_data_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

meta_plot_data_export_table <- function(result) {
  series <- result$chart_payload$series
  rows <- unlist(lapply(series, function(entry) {
    lapply(entry$values, function(point) {
      data.frame(
        Sample = entry$name,
        Source_File = if (is.null(entry$originalName)) entry$name else entry$originalName,
        Position = as.numeric(point$x),
        Density = as.numeric(point$density),
        Confidence_Down = if (is.null(point$confidenceDown)) NA_real_ else as.numeric(point$confidenceDown),
        Confidence_Up = if (is.null(point$confidenceUp)) NA_real_ else as.numeric(point$confidenceUp),
        stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)

  export_data <- do.call(rbind, rows)
  export_data$Species <- result$species
  export_data[, c("Species", "Sample", "Source_File", "Position", "Density", "Confidence_Down", "Confidence_Up")]
}
