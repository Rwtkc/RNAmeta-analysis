peak_distribution_compute_feature_counts <- function(annotation_result) {
  annotation <- annotation_result$annotation
  annotation_stopcodon <- annotation_result$annotation_stopcodon
  annotation_startcodon <- annotation_result$annotation_startcodon

  annotation[grep("CDS \\(", annotation)] <- "CDS"
  annotation[grep("Intron \\(", annotation)] <- "Intron"
  annotation[grep("Promoter \\(", annotation)] <- "Promoter"
  annotation[grep("UTR5 \\(", annotation)] <- "UTR5"
  annotation[grep("UTR3 \\(", annotation)] <- "UTR3"

  annotation_stopcodon <- annotation_stopcodon[grep("Stop Codon \\(", annotation_stopcodon)]
  annotation_stopcodon[grep("Stop Codon \\(", annotation_stopcodon)] <- "Stop Codon"

  annotation_startcodon <- annotation_startcodon[grep("^Start Codon", annotation_startcodon)]
  annotation_startcodon[grep("^Start Codon", annotation_startcodon)] <- "Start Codon"

  annotation_table <- table(annotation)
  stop_table <- table(annotation_stopcodon)
  start_table <- table(annotation_startcodon)
  annotation_table <- c(annotation_table, start_table, stop_table)
  annotation_df <- data.table::as.data.table(annotation_table, keep.rownames = TRUE)
  data.table::setnames(annotation_df, c("Feature", "Frequency"))
  annotation_df$Feature <- factor(annotation_df$Feature, levels = peak_distribution_feature_order())
  data.table::setDT(annotation_df[order(annotation_df$Feature), ])
}

peak_distribution_prepare_plot_data <- function(feature_counts, selected_features = peak_distribution_feature_order(), sample_labels = NULL) {
  plot_data <- subset(feature_counts, Feature %in% selected_features)

  if (!is.null(sample_labels) && length(sample_labels) == length(unique(plot_data$Sample))) {
    plot_data$Sample <- factor(
      plot_data$Sample,
      levels = unique(plot_data$Sample),
      labels = sample_labels
    )
  }

  plot_data
}
