peak_exon_size_prepare_plot_data <- function(exon_size_data) {
  if (is.null(exon_size_data) || !nrow(exon_size_data)) {
    return(data.table::data.table(
      exonID = character(),
      tx_name = character(),
      rank = numeric(),
      length = numeric(),
      exonLoc = character(),
      Sample = character()
    ))
  }

  plot_data <- data.table::copy(exon_size_data)
  plot_data[, length := suppressWarnings(as.numeric(length))]
  plot_data <- plot_data[is.finite(length) & length > 0]
  plot_data[, exonLoc := factor(as.character(exonLoc), levels = peak_exon_location_levels())]
  plot_data[order(exonLoc, Sample, rank)]
}
