run_peak_exon_num_analysis <- function(upload_context, controls = peak_exon_num_default_controls(), progress = NULL) {
  controls <- utils::modifyList(peak_exon_num_default_controls(), controls)
  peak_exon_num_progress_update(progress, 8, "Preparing saved BED context for exon-count analysis")

  collected_data <- peak_exon_collect_analysis_data(
    upload_context = upload_context,
    module_label = peak_exon_num_title(),
    progress = progress,
    include_exon_size = FALSE,
    include_transcript_stats = TRUE
  )
  plot_data <- peak_exon_num_prepare_plot_data(collected_data$transcript_stats)
  peak_exon_num_progress_update(progress, 100, "Building chart payload and finalizing exon-count summaries")

  list(
    chart_payload = peak_exon_num_build_chart_payload(plot_data, sample_metadata = collected_data$sample_metadata),
    peak_exon_num_data = plot_data,
    species = collected_data$species,
    sample_name = collected_data$sample_original_names,
    sample_display_names = collected_data$sample_display_names,
    sample_count = collected_data$sample_count,
    interval_count = collected_data$interval_count,
    transcript_count = nrow(plot_data)
  )
}
