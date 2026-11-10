gene_type_build_chart_payload <- function(plot_data, sample_metadata = NULL) {
  if (!nrow(plot_data)) {
    return(NULL)
  }

  categories <- unique(as.character(plot_data$GeneType))
  samples <- unique(as.character(plot_data$Sample))
  palette <- peak_distribution_group_palette()
  max_value <- max(plot_data$Frequency, na.rm = TRUE)

  list(
    type = "grouped_bar",
    title = "Gene Type",
    xLabel = "",
    yLabel = "Frequency",
    categories = categories,
    hoverMinHeight = 18,
    showLabels = FALSE,
    yTicks = scales::pretty_breaks(n = 5)(c(0, max_value * 1.12)),
    yDomain = c(0, max_value * 1.12),
    series = lapply(seq_along(samples), function(index) {
      sample_name <- samples[[index]]
      sample_data <- plot_data[plot_data$Sample == sample_name, ]
      values <- sample_data$Frequency[match(categories, as.character(sample_data$GeneType))]
      values[is.na(values)] <- 0
      original_name <- rnameta_lookup_original_name(sample_metadata, sample_name, fallback = sample_name)

      list(
        name = sample_name,
        originalName = original_name,
        color = palette[((index - 1) %% length(palette)) + 1],
        values = as.numeric(values)
      )
    })
  )
}
