peak_distribution_group_palette <- function() {
  c("#de7a61", "#3f3f5b", "#81b29d", "#ffb5b0")
}

peak_distribution_feature_palette <- function() {
  c(
    "Promoter" = "#accbe8",
    "UTR5" = "#dda43d",
    "Start Codon" = "#519475",
    "CDS" = "#5d2b84",
    "Stop Codon" = "#eddb31",
    "UTR3" = "#4c9140",
    "Intron" = "#5b98e8",
    "Intergenic" = "#be3731"
  )
}

peak_distribution_build_chart_payload <- function(plot_data, chart_type = "bar", show_bar_labels = FALSE, series_metadata = NULL) {
  if (!nrow(plot_data)) {
    return(NULL)
  }

  if (identical(chart_type, "bar")) {
    categories <- unique(as.character(plot_data$Feature))
    samples <- unique(as.character(plot_data$Sample))
    palette <- peak_distribution_group_palette()
    max_value <- max(plot_data$Frequency, na.rm = TRUE)

      return(list(
        type = "grouped_bar",
        title = "Peak Distribution",
        xLabel = "",
        yLabel = "Frequency",
      categories = categories,
      showLabels = isTRUE(show_bar_labels),
      yTicks = scales::pretty_breaks(n = 5)(c(0, max_value * 1.12)),
      yDomain = c(0, max_value * 1.12),
      series = lapply(seq_along(samples), function(index) {
        sample_name <- samples[[index]]
        sample_data <- plot_data[plot_data$Sample == sample_name, ]
        values <- sample_data$Frequency[match(categories, as.character(sample_data$Feature))]
        values[is.na(values)] <- 0
        metadata_index <- NULL
        if (!is.null(series_metadata) && length(series_metadata)) {
          metadata_index <- which(vapply(series_metadata, function(entry) identical(entry$display_name, sample_name), logical(1)))[1]
        }
        original_name <- if (!is.na(metadata_index) && length(metadata_index) == 1L) series_metadata[[metadata_index]]$original_name else sample_name

        list(
          name = sample_name,
          originalName = original_name,
          color = palette[((index - 1) %% length(palette)) + 1],
          values = as.numeric(values)
        )
      })
    ))
  }

  categories <- unique(as.character(plot_data$Sample))
  features <- unique(as.character(plot_data$Feature))
  totals <- tapply(plot_data$Frequency, plot_data$Sample, sum)
  palette <- peak_distribution_feature_palette()

  list(
    type = "stacked_bar",
    title = "Peak Distribution Composition",
    xLabel = "",
    yLabel = "Percentage",
    categories = categories,
    showLabels = FALSE,
    yTicks = seq(0, 100, by = 20),
    yDomain = c(0, 100),
    series = lapply(features, function(feature_name) {
      feature_data <- plot_data[plot_data$Feature == feature_name, ]
      values <- feature_data$Frequency[match(categories, as.character(feature_data$Sample))]
      values[is.na(values)] <- 0
      percentages <- values / as.numeric(totals[categories]) * 100
      percentages[is.na(percentages)] <- 0

      list(
        name = feature_name,
        color = if (is.null(palette[[feature_name]])) "#859b7a" else unname(palette[[feature_name]]),
        values = as.numeric(percentages)
      )
    })
  )
}
