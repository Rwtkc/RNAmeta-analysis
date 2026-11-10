gene_type_runtime_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$annotation_bundle <- new.env(hash = TRUE, parent = emptyenv())
  env$feature_ranges <- new.env(hash = TRUE, parent = emptyenv())
  env$analysis_result <- new.env(hash = TRUE, parent = emptyenv())
  env
})

gene_type_runtime_cache_key <- function(...) {
  paste(vapply(list(...), as.character, character(1)), collapse = "::")
}

gene_type_cached_files_from_context <- function(upload_context) {
  if (!is.null(upload_context$saved_files) && length(upload_context$saved_files)) {
    return(upload_context$saved_files)
  }

  if (!is.null(upload_context$file_path) && nzchar(upload_context$file_path)) {
    return(list(
      list(
        file_name = if (is.null(upload_context$file_name)) basename(upload_context$file_path) else upload_context$file_name,
        file_path = upload_context$file_path
      )
    ))
  }

  list()
}

gene_type_cached_annotation_bundle <- function(species) {
  key <- gene_type_runtime_cache_key(species)

  if (exists(key, envir = gene_type_runtime_cache_env$annotation_bundle, inherits = FALSE)) {
    return(get(key, envir = gene_type_runtime_cache_env$annotation_bundle, inherits = FALSE))
  }

  bundle <- gene_type_load_annotation_bundle_uncached(species)
  assign(key, bundle, envir = gene_type_runtime_cache_env$annotation_bundle)
  bundle
}

gene_type_cached_feature_ranges <- function(species) {
  key <- gene_type_runtime_cache_key(species)

  if (exists(key, envir = gene_type_runtime_cache_env$feature_ranges, inherits = FALSE)) {
    return(get(key, envir = gene_type_runtime_cache_env$feature_ranges, inherits = FALSE))
  }

  bundle <- gene_type_cached_annotation_bundle(species)
  feature_ranges <- gene_type_prepare_feature_ranges(bundle$gff, bundle$txdb)
  assign(key, feature_ranges, envir = gene_type_runtime_cache_env$feature_ranges)
  feature_ranges
}

gene_type_cached_analysis_key <- function(upload_context) {
  saved_files <- gene_type_cached_files_from_context(upload_context)
  file_parts <- vapply(saved_files, function(file_entry) {
    file_info <- file.info(file_entry$file_path)
    file_size <- if (is.na(file_info$size)) 0 else file_info$size
    file_mtime <- if (is.na(file_info$mtime)) 0 else as.numeric(file_info$mtime)
    paste(
      normalizePath(file_entry$file_path, winslash = "/", mustWork = TRUE),
      as.numeric(file_size),
      file_mtime,
      sep = "|"
    )
  }, character(1))

  gene_type_runtime_cache_key(
    upload_context$species,
    paste(file_parts, collapse = "::")
  )
}

gene_type_cached_analysis_result <- function(upload_context) {
  key <- gene_type_cached_analysis_key(upload_context)

  if (!exists(key, envir = gene_type_runtime_cache_env$analysis_result, inherits = FALSE)) {
    return(NULL)
  }

  get(key, envir = gene_type_runtime_cache_env$analysis_result, inherits = FALSE)
}

gene_type_store_analysis_result <- function(upload_context, result) {
  key <- gene_type_cached_analysis_key(upload_context)
  assign(key, result, envir = gene_type_runtime_cache_env$analysis_result)
  result
}

gene_type_runtime_cache_status <- function() {
  list(
    annotation_bundle_keys = ls(envir = gene_type_runtime_cache_env$annotation_bundle),
    feature_range_keys = ls(envir = gene_type_runtime_cache_env$feature_ranges),
    analysis_result_keys = ls(envir = gene_type_runtime_cache_env$analysis_result)
  )
}

gene_type_runtime_cache_reset <- function() {
  rm(list = ls(envir = gene_type_runtime_cache_env$annotation_bundle), envir = gene_type_runtime_cache_env$annotation_bundle)
  rm(list = ls(envir = gene_type_runtime_cache_env$feature_ranges), envir = gene_type_runtime_cache_env$feature_ranges)
  rm(list = ls(envir = gene_type_runtime_cache_env$analysis_result), envir = gene_type_runtime_cache_env$analysis_result)
  invisible(TRUE)
}
