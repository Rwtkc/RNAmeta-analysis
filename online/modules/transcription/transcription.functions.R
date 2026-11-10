transcription_default_controls <- function() {
  list(
    flank_size = 1000
  )
}

transcription_panel_titles <- function() {
  c(
    "Transcription_start_site" = "Distribution on transcription start site",
    "Transcription_end_site" = "Distribution on transcription end site"
  )
}

transcription_has_boundary_signal <- function(density_hits) {
  !is.null(density_hits) && nrow(density_hits) > 0
}

transcription_matrix_has_signal <- function(tag_matrix) {
  !is.null(tag_matrix) &&
    length(tag_matrix) &&
    nrow(tag_matrix) > 0 &&
    any(is.finite(tag_matrix) & tag_matrix > 0)
}

transcription_empty_state_message <- function() {
  paste(
    "The analysis completed, but no peaks overlapped the transcription boundary windows.",
    "Try checking species selection, chromosome naming, or whether the current BED file is expected near TSS/TES regions."
  )
}

transcription_heatmap_empty_message <- function() {
  paste(
    "Boundary hits were detected, but the transcription heatmap panels contain no visible signal for the sampled windows.",
    "You can still inspect the density curves above, or try a different flank size or BED file if you expect stronger TSS/TES enrichment."
  )
}
