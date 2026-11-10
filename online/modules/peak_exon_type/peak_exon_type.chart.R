peak_exon_type_build_chart_payload <- function(plot_data, sample_metadata = NULL) {
  if (!nrow(plot_data)) {
    return(NULL)
  }

  categories <- unname(unique(as.character(plot_data$Sample)))
  exon_locations <- peak_exon_location_levels()
  exon_locations <- exon_locations[exon_locations %in% unique(as.character(plot_data$exonLoc))]
  palette <- c(
    "First" = "#a30545",
    "Middle" = "#f26e42",
    "Last" = "#7ecba4"
  )

  list(
    type = "stacked_bar",
    title = peak_exon_type_title(),
    xLabel = "",
    yLabel = "Percentage",
    tooltipSeriesLabel = "Exon Type",
    tooltipValueLabel = "Percentage",
    categories = categories,
    categoryOriginalNames = as.list(stats::setNames(
      vapply(categories, function(sample_name) rnameta_lookup_original_name(sample_metadata, sample_name, fallback = sample_name), character(1)),
      categories
    )),
    showLabels = FALSE,
    yTicks = seq(0, 100, by = 20),
    yDomain = c(0, 100),
    series = lapply(exon_locations, function(location_name) {
      location_data <- plot_data[as.character(exonLoc) == location_name]
        values <- location_data$Percentage[match(categories, as.character(location_data$Sample))]
        values[is.na(values)] <- 0

        list(
          name = location_name,
          color = unname(palette[[location_name]]),
          values = unname(as.numeric(values))
        )
      })
  )
}
