gene_matrix_export_defaults <- function() {
  list(
    format = "png",
    width = 10400,
    height = 7200,
    dpi = 900
  )
}

gene_matrix_export_config <- function(id) {
  ns <- NS(id)
  defaults <- gene_matrix_export_defaults()

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

gene_matrix_export_settings <- function(input) {
  defaults <- gene_matrix_export_defaults()
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

gene_matrix_export_filename <- function(extension) {
  sprintf("gene_matrix_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

gene_matrix_data_export_filename <- function(extension) {
  sprintf("gene_matrix_plot_data_%s.%s", format(Sys.Date(), "%Y%m%d"), extension)
}

gene_matrix_data_export_table <- function(result) {
  sample_metadata <- Map(
    function(display_name, original_name) {
      list(display_name = display_name, original_name = original_name)
    },
    result$sample_display_names,
    result$sample_name
  )
  rows <- gene_matrix_intersection_rows(result$gene_sets, result$sample_display_names, sample_metadata = sample_metadata)

  data.frame(
    Species = rep(result$species, length(rows)),
    Display_Sets = vapply(rows, function(entry) paste(unlist(entry$sets), collapse = " / "), character(1)),
    Source_Files = vapply(rows, function(entry) paste(unlist(entry$originalSets), collapse = " / "), character(1)),
    Intersection_Size = vapply(rows, function(entry) as.numeric(entry$size), numeric(1)),
    Genes = vapply(rows, function(entry) paste(unlist(entry$genes), collapse = "; "), character(1)),
    stringsAsFactors = FALSE
  )
}
