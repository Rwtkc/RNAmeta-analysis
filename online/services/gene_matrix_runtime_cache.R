gene_matrix_runtime_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$transcript_meta <- new.env(hash = TRUE, parent = emptyenv())
  env$prepared_features <- new.env(hash = TRUE, parent = emptyenv())
  env
})

gene_matrix_runtime_cache_key <- function(...) {
  paste(vapply(list(...), as.character, character(1)), collapse = "::")
}

gene_matrix_cached_transcript_meta <- function(species, builder) {
  key <- gene_matrix_runtime_cache_key(species)

  if (exists(key, envir = gene_matrix_runtime_cache_env$transcript_meta, inherits = FALSE)) {
    return(get(key, envir = gene_matrix_runtime_cache_env$transcript_meta, inherits = FALSE))
  }

  transcript_meta <- builder()
  assign(key, transcript_meta, envir = gene_matrix_runtime_cache_env$transcript_meta)
  transcript_meta
}

gene_matrix_cached_prepared_features <- function(species, minimal_component_length, builder) {
  key <- gene_matrix_runtime_cache_key(species, minimal_component_length)

  if (exists(key, envir = gene_matrix_runtime_cache_env$prepared_features, inherits = FALSE)) {
    return(get(key, envir = gene_matrix_runtime_cache_env$prepared_features, inherits = FALSE))
  }

  prepared_features <- builder()
  assign(key, prepared_features, envir = gene_matrix_runtime_cache_env$prepared_features)
  prepared_features
}

gene_matrix_runtime_cache_status <- function() {
  list(
    transcript_meta_keys = ls(envir = gene_matrix_runtime_cache_env$transcript_meta),
    prepared_feature_keys = ls(envir = gene_matrix_runtime_cache_env$prepared_features)
  )
}

gene_matrix_runtime_cache_reset <- function() {
  rm(list = ls(envir = gene_matrix_runtime_cache_env$transcript_meta), envir = gene_matrix_runtime_cache_env$transcript_meta)
  rm(list = ls(envir = gene_matrix_runtime_cache_env$prepared_features), envir = gene_matrix_runtime_cache_env$prepared_features)
  invisible(TRUE)
}
