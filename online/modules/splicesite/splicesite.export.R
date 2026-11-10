splicesite_export_defaults <- function() {
  list(
    format = "png",
    width = 3000,
    height = 1400,
    dpi = 300
  )
}

splicesite_export_config <- function(id) {
  ns <- NS(id)
  defaults <- splicesite_export_defaults()

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

splicesite_export_settings <- function(input) {
  defaults <- splicesite_export_defaults()
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

splicesite_export_filename <- function(extension) {
  sprintf("splicesite_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

splicesite_data_export_filename <- function(extension) {
  sprintf("splicesite_plot_data_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}
