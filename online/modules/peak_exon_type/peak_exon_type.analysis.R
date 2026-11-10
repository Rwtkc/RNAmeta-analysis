run_peak_exon_type_analysis <- function(upload_context, controls = peak_exon_type_default_controls(), progress = NULL) {
  controls <- utils::modifyList(peak_exon_type_default_controls(), controls)
  peak_exon_type_progress_update(progress, 8, "Preparing saved BED context for exon-type analysis")

  collected_data <- peak_exon_collect_analysis_data(
    upload_context = upload_context,
    module_label = peak_exon_type_title(),
    progress = progress,
    include_exon_size = TRUE,
    include_transcript_stats = FALSE
  )
  plot_data <- peak_exon_type_prepare_plot_data(collected_data$exon_size_data)
  peak_exon_type_progress_update(progress, 100, "Building chart payload and finalizing exon-type summaries")

  list(
    chart_payload = peak_exon_type_build_chart_payload(plot_data, sample_metadata = collected_data$sample_metadata),
    peak_exon_type_data = plot_data,
    species = collected_data$species,
    sample_name = collected_data$sample_original_names,
    sample_display_names = collected_data$sample_display_names,
    sample_count = collected_data$sample_count,
    interval_count = collected_data$interval_count,
    exon_hit_count = sum(plot_data$Count, na.rm = TRUE)
  )
}
