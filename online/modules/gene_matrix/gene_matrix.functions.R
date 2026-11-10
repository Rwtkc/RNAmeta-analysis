gene_matrix_saved_files <- function(upload_context) {
  if (!is.null(upload_context$saved_files) && length(upload_context$saved_files)) {
    return(upload_context$saved_files)
  }

  if (!is.null(upload_context$file_path) && nzchar(upload_context$file_path)) {
    return(list(list(
      file_name = if (is.null(upload_context$file_name)) basename(upload_context$file_path) else upload_context$file_name,
      file_path = upload_context$file_path
    )))
  }

  list()
}

gene_matrix_sample_count_ok <- function(sample_count) {
  isTRUE(sample_count >= 2 && sample_count <= 5)
}

gene_matrix_sample_labels <- function(saved_files) {
  rnameta_sample_display_names(rnameta_build_sample_metadata(saved_files))
}

gene_matrix_venn_palette <- function() {
  c("#FF0000", "#ff8080", "#a1c530", "#6676a0", "#33b39f")
}

gene_matrix_required_sample_message <- function() {
  "Gene Matrix requires 2 to 5 saved BED files in Upload / Run."
}

gene_matrix_plot_dimensions <- function() {
  list(width = 1800, height = 1350, res = 180)
}
