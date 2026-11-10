peak_distribution_build_plot <- function(plot_data, chart_type = "bar", show_bar_labels = FALSE, palette = c("#de7a61", "#3f3f5b", "#81b29d", "#ffb5b0")) {
  if (identical(chart_type, "bar")) {
    offset <- max(plot_data$Frequency, na.rm = TRUE) * 0.03
    plot_data$label_y <- plot_data$Frequency + offset
    legend_position <- if (length(unique(plot_data$Sample)) == 1) "none" else "right"

    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Feature, y = Frequency, fill = Sample)) +
      ggplot2::geom_bar(stat = "identity", color = "black", position = ggplot2::position_dodge(width = 1)) +
      ggplot2::scale_y_continuous(
        labels = function(y) paste0(round(y / 1000, 1), "k"),
        expand = ggplot2::expansion(mult = c(0, 0.2))
      ) +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::labs(y = "Frequency", x = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5),
        axis.text.y = ggplot2::element_text(angle = 0, hjust = 1),
        axis.title.y = ggplot2::element_text(vjust = 2.5, size = 15),
        panel.grid = ggplot2::element_blank(),
        legend.position = legend_position
      )

    if (isTRUE(show_bar_labels)) {
      plot <- plot + ggplot2::geom_text(
        ggplot2::aes(y = label_y, label = Frequency),
        position = ggplot2::position_dodge(width = 1),
        size = 3,
        show.legend = FALSE
      )
    }

    return(plot)
  }

  colors <- c(
    "Promoter" = "#accbe8",
    "UTR5" = "#dda43d",
    "Start Codon" = "#519475",
    "CDS" = "#5d2b84",
    "Stop Codon" = "#eddb31",
    "UTR3" = "#4c9140",
    "Intron" = "#5b98e8",
    "Intergenic" = "#be3731"
  )

  plot_data$TotalFrequency <- ave(plot_data$Frequency, plot_data$Sample, FUN = sum)
  plot_data$Percentage <- (plot_data$Frequency / plot_data$TotalFrequency) * 100

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Sample, y = Percentage, fill = Feature)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", color = "black", width = 0.3) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(scale = 1),
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::labs(y = "", x = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(angle = 0, hjust = 0.1),
      axis.title.x = ggplot2::element_text(vjust = -1, size = 15),
      axis.title.y = ggplot2::element_text(vjust = 2.5, size = 15),
      panel.grid.major = ggplot2::element_line(color = "gray", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_line(color = "gray", linewidth = 0.25),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}
