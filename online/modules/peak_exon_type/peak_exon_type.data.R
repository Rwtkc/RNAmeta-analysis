peak_exon_type_prepare_plot_data <- function(exon_size_data) {
  if (is.null(exon_size_data) || !nrow(exon_size_data)) {
    return(data.table::data.table(
      Sample = character(),
      exonLoc = character(),
      Count = numeric(),
      Percentage = numeric()
    ))
  }

  type_counts <- data.table::copy(exon_size_data)[
    ,
    .(Count = .N),
    by = .(Sample, exonLoc = as.character(exonLoc))
  ]
  type_counts[, total := sum(Count), by = Sample]
  type_counts[, Percentage := data.table::fifelse(total > 0, Count / total * 100, 0)]
  type_counts[, total := NULL]
  type_counts[, exonLoc := factor(exonLoc, levels = peak_exon_location_levels())]
  type_counts[order(Sample, exonLoc)]
}
