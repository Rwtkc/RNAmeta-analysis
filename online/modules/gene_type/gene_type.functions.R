gene_type_sample_label <- function(file_info) {
  rnameta_sample_original_name(file_info)
}

gene_type_default_controls <- function() {
  list()
}

gene_type_normalize_type_values <- function(type_values) {
  type_values <- as.character(type_values)
  type_values <- type_values[!is.na(type_values) & nzchar(type_values)]
  type_values
}
