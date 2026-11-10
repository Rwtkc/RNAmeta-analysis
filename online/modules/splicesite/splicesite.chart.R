splicesite_build_chart_payload <- function(density_data, sample_metadata, flank_size) {
  density_titles <- splicesite_panel_titles()
  density_panels <- lapply(names(density_titles), function(panel_key) {
    x_domain <- if (identical(panel_key, "5PSS")) c(-flank_size, flank_size * 2) else c(-flank_size * 2, flank_size)

    site_profile_density_panel(
      panel_title = density_titles[[panel_key]],
      panel_key = panel_key,
      density_data = density_data,
      sample_metadata = sample_metadata,
      x_domain = x_domain,
      x_label = "Region around the site (0 and 5' -> 3' direction)",
      y_label = "Density"
    )
  })

  list(
    title = "Splicesite",
    panels = density_panels,
    legend = site_profile_legend_entries(sample_metadata)
  )
}

splicesite_data_export_table <- function(result) {
  site_profile_density_export_table(
    density_data = result$density_data,
    sample_metadata = result$sample_metadata,
    module_name = "Splicesite"
  )
}
