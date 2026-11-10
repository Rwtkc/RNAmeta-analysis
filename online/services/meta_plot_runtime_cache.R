meta_plot_runtime_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$annotation_bundle <- new.env(hash = TRUE, parent = emptyenv())
  env$annotation_features <- new.env(hash = TRUE, parent = emptyenv())
  env$guitar_txdb <- new.env(hash = TRUE, parent = emptyenv())
  env$mapped_transcripts <- new.env(hash = TRUE, parent = emptyenv())
  env$sample_points <- new.env(hash = TRUE, parent = emptyenv())
  env
})

meta_plot_reset_runtime_caches <- function() {
  rm(list = ls(meta_plot_runtime_cache_env$annotation_bundle), envir = meta_plot_runtime_cache_env$annotation_bundle)
  rm(list = ls(meta_plot_runtime_cache_env$annotation_features), envir = meta_plot_runtime_cache_env$annotation_features)
  rm(list = ls(meta_plot_runtime_cache_env$guitar_txdb), envir = meta_plot_runtime_cache_env$guitar_txdb)
  rm(list = ls(meta_plot_runtime_cache_env$mapped_transcripts), envir = meta_plot_runtime_cache_env$mapped_transcripts)
  rm(list = ls(meta_plot_runtime_cache_env$sample_points), envir = meta_plot_runtime_cache_env$sample_points)
  invisible(TRUE)
}

meta_plot_file_fingerprint <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return("missing")
  }

  info <- file.info(path)
  paste(
    normalizePath(path, winslash = "/"),
    unname(info$size[[1]]),
    as.numeric(unname(info$mtime[[1]])),
    sep = "|"
  )
}

meta_plot_extract_controls_key <- function(species, controls) {
  paste(
    species,
    controls$pltTxType,
    controls$txfiveutrMinLength,
    controls$txcdsMinLength,
    controls$txthreeutrMinLength,
    controls$txlongNcrnaMinLength,
    controls$txlncrnaOverlapmrna,
    controls$txpromoterLength,
    controls$txtailLength,
    controls$stAmblguity,
    controls$txPrimaryOnly,
    sep = "|"
  )
}

meta_plot_sample_controls_key <- function(source_file, controls, extract_key) {
  paste(
    meta_plot_file_fingerprint(source_file),
    extract_key,
    controls$mapFilterTranscript,
    controls$stSampleModle,
    controls$stSampleNum,
    sep = "|"
  )
}

meta_plot_mapping_controls_key <- function(source_file, controls, extract_key) {
  paste(
    meta_plot_file_fingerprint(source_file),
    extract_key,
    controls$mapFilterTranscript,
    sep = "|"
  )
}

meta_plot_cached_annotation_bundle <- function(species, loader) {
  cache <- meta_plot_runtime_cache_env$annotation_bundle

  if (exists(species, envir = cache, inherits = FALSE)) {
    return(get(species, envir = cache, inherits = FALSE))
  }

  bundle <- loader()
  assign(species, bundle, envir = cache)
  bundle
}

meta_plot_cached_annotation_features <- function(species, builder) {
  cache <- meta_plot_runtime_cache_env$annotation_features

  if (exists(species, envir = cache, inherits = FALSE)) {
    return(get(species, envir = cache, inherits = FALSE))
  }

  features <- builder()
  assign(species, features, envir = cache)
  features
}

meta_plot_cached_guitar_txdb <- function(species, controls, builder) {
  key <- meta_plot_extract_controls_key(species, controls)
  cache <- meta_plot_runtime_cache_env$guitar_txdb

  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }

  guitar_txdb <- builder()
  assign(key, guitar_txdb, envir = cache)
  guitar_txdb
}

meta_plot_cached_mapped_transcripts <- function(source_file, controls, extract_key, mapper) {
  key <- meta_plot_mapping_controls_key(source_file, controls, extract_key)
  cache <- meta_plot_runtime_cache_env$mapped_transcripts

  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }

  mapped <- mapper()
  assign(key, mapped, envir = cache)
  mapped
}

meta_plot_cached_sample_points <- function(source_file, controls, extract_key, sampler) {
  key <- meta_plot_sample_controls_key(source_file, controls, extract_key)
  cache <- meta_plot_runtime_cache_env$sample_points

  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }

  sampled <- sampler()
  assign(key, sampled, envir = cache)
  sampled
}
