peak_gene_size_build_chart_payload <- function(plot_data, sample_metadata = NULL) {
  if (!nrow(plot_data)) {
    return(NULL)
  }

  samples <- unique(as.character(plot_data$Sample))
  palette <- peak_gene_size_group_palette()

  list(
    type = "boxplot",
    title = "Peak Gene Size",
    yLabel = "log2(txLength bp)",
    groups = lapply(seq_along(samples), function(index) {
      sample_name <- samples[[index]]
      sample_data <- plot_data[plot_data$Sample == sample_name, ]
      original_name <- rnameta_lookup_original_name(sample_metadata, sample_name, fallback = sample_name)

      list(
        name = sample_name,
        originalName = original_name,
        color = palette[((index - 1) %% length(palette)) + 1],
        values = as.numeric(sample_data$tx_len)
      )
    })
  )
}
