peak_distribution_runtime_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$annotation_bundle <- new.env(hash = TRUE, parent = emptyenv())
  env$txdb_features <- new.env(hash = TRUE, parent = emptyenv())
  env
})

peak_distribution_runtime_cache_key <- function(...) {
  paste(vapply(list(...), as.character, character(1)), collapse = "::")
}

peak_distribution_cached_annotation_bundle <- function(species) {
  key <- peak_distribution_runtime_cache_key(species)

  if (exists(key, envir = peak_distribution_runtime_cache_env$annotation_bundle, inherits = FALSE)) {
    return(get(key, envir = peak_distribution_runtime_cache_env$annotation_bundle, inherits = FALSE))
  }

  bundle <- peak_distribution_load_annotation_bundle(species)
  bundle$txlens_summary <- data.table::as.data.table(bundle$txlens)[, lapply(.SD, max), by = gene_id]
  assign(key, bundle, envir = peak_distribution_runtime_cache_env$annotation_bundle)
  bundle
}

peak_distribution_cached_txdb_features <- function(species, minimal_component_length = 100) {
  key <- peak_distribution_runtime_cache_key(species, minimal_component_length)

  if (exists(key, envir = peak_distribution_runtime_cache_env$txdb_features, inherits = FALSE)) {
    return(get(key, envir = peak_distribution_runtime_cache_env$txdb_features, inherits = FALSE))
  }

  bundle <- peak_distribution_cached_annotation_bundle(species)
  features <- peak_distribution_build_txdb_features(
    txdb = bundle$txdb,
    txlens_summary = bundle$txlens_summary,
    minimal_component_length = minimal_component_length
  )
  assign(key, features, envir = peak_distribution_runtime_cache_env$txdb_features)
  features
}

peak_distribution_runtime_cache_status <- function() {
  list(
    annotation_bundle_keys = ls(envir = peak_distribution_runtime_cache_env$annotation_bundle),
    txdb_feature_keys = ls(envir = peak_distribution_runtime_cache_env$txdb_features)
  )
}

peak_distribution_runtime_cache_reset <- function() {
  rm(list = ls(envir = peak_distribution_runtime_cache_env$annotation_bundle), envir = peak_distribution_runtime_cache_env$annotation_bundle)
  rm(list = ls(envir = peak_distribution_runtime_cache_env$txdb_features), envir = peak_distribution_runtime_cache_env$txdb_features)
  invisible(TRUE)
}
