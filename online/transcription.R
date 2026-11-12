transcription_default_controls <- function() {
  list(
    flank_size = 1000
  )
}

site_profile_prepare_transcription_windows <- function(species, txdb_features, txlens, flank_size = 1000) {
  site_profile_cached_boundary_windows("transcription", species, flank_size, function() {
    # One true transcript span per isoform, including its UTRs.
    transcripts <- BiocGenerics::unlist(range(txdb_features$exons))
    boundary_window <- function(fix) {
      anchor <- GenomicRanges::resize(transcripts, width = 1, fix = fix)
      symmetric <- GenomicRanges::resize(anchor, width = 2 * flank_size + 1, fix = "center")
      GenomicRanges::resize(symmetric, width = 2 * flank_size, fix = "start")
    }
    coverage_tss <- boundary_window("start")
    coverage_tes <- boundary_window("end")
    tx_names <- unique(as.character(data.table::as.data.table(txlens)$tx_name))
    heatmap_tss <- coverage_tss[names(coverage_tss) %in% tx_names]
    heatmap_tes <- coverage_tes[names(coverage_tes) %in% tx_names]

    list(
      coverage_tss = coverage_tss,
      coverage_tes = coverage_tes,
      heatmap_tss = heatmap_tss,
      heatmap_tes = heatmap_tes,
      x_positions = seq(-flank_size, flank_size - 1),
      x_domain = c(-flank_size, flank_size - 1),
      rel_offset = flank_size + 1L,
      flank_size = flank_size
    )
  })
}
