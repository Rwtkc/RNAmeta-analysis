peak_distribution_feature_order <- function() {
  c("Promoter", "UTR5", "Start Codon", "CDS", "Stop Codon", "UTR3", "Intron", "Intergenic")
}

peak_distribution_annotation_priority <- function() {
  c("Stopcodon", "Startcodon", "CDS", "5UTR", "3UTR", "Intron", "Transcript", "Promoter", "Intergenic")
}

peak_distribution_default_controls <- function() {
  list(
    selected_features = peak_distribution_feature_order(),
    detail_output_dir = NULL
  )
}

peak_distribution_feature_choices <- function() {
  lapply(
    peak_distribution_feature_order(),
    function(feature) {
      list(value = feature, label = feature)
    }
  )
}

peak_distribution_controls_config <- function(id, defaults = peak_distribution_default_controls()) {
  ns <- NS(id)

  list(
    ids = list(
      selected_features = ns("selected_features")
    ),
    defaults = list(
      selected_features = defaults$selected_features
    ),
    choices = list(
      features = peak_distribution_feature_choices()
    )
  )
}
