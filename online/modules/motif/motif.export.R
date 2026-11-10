motif_export_defaults <- function() {
  list(
    format = "png",
    width = 2200,
    height = 1400,
    dpi = 300
  )
}

motif_export_config <- function(id) {
  ns <- NS(id)
  defaults <- motif_export_defaults()

  list(
    ids = list(
      format = ns("export_format"),
      width = ns("export_width"),
      height = ns("export_height"),
      dpi = ns("export_dpi"),
      trigger = ns("export_plot")
    ),
    defaults = defaults,
    sectioned = TRUE,
    choices = list(
      format = list(
        list(label = "PNG", value = "png"),
        list(label = "PDF", value = "pdf")
      )
    )
  )
}

motif_export_settings <- function(input) {
  defaults <- motif_export_defaults()

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

motif_export_filename <- function(engine, sample_name, format) {
  safe_engine <- gsub("[^A-Za-z0-9_-]+", "_", tolower(as.character(engine)[[1]]))
  safe_sample <- gsub("[^A-Za-z0-9._-]+", "_", basename(as.character(sample_name)[[1]]))
  sprintf(
    "motif_%s_%s_%s.zip",
    safe_engine,
    safe_sample,
    tolower(as.character(format)[[1]])
  )
}
