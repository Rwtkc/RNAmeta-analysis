peak_exon_size_build_chart_payload <- function(plot_data, sample_metadata = NULL) {
  if (!nrow(plot_data)) {
    return(NULL)
  }

  samples <- unname(unique(as.character(plot_data$Sample)))
  locations <- peak_exon_location_levels()
  locations <- unname(locations[locations %in% unique(as.character(plot_data$exonLoc))])
  palette <- peak_exon_group_palette()

  list(
    type = "boxplot_facet",
    title = peak_exon_size_title(),
    yLabel = "log2(Length bp)",
    scaleTransform = "log2",
    facets = lapply(locations, function(location_name) {
      facet_data <- plot_data[as.character(exonLoc) == location_name]

      list(
        name = location_name,
        groups = lapply(seq_along(samples), function(index) {
          sample_name <- samples[[index]]
          sample_data <- facet_data[Sample == sample_name]
          original_name <- rnameta_lookup_original_name(sample_metadata, sample_name, fallback = sample_name)

          list(
            name = sample_name,
            originalName = original_name,
            color = palette[((index - 1) %% length(palette)) + 1],
            values = unname(as.numeric(sample_data$length))
          )
        })
      )
    })
  )
}
