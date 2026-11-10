peak_exon_num_prepare_plot_data <- function(transcript_stats) {
  if (is.null(transcript_stats) || !nrow(transcript_stats)) {
    return(data.table::data.table(
      tx_name = character(),
      tx_len = numeric(),
      nexon = numeric(),
      Sample = character()
    ))
  }

  plot_data <- data.table::copy(transcript_stats)
  plot_data[, nexon := suppressWarnings(as.numeric(nexon))]
  plot_data <- plot_data[is.finite(nexon) & nexon > 0]
  plot_data[order(Sample, tx_name)]
}
