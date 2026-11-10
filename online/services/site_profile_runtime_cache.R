site_profile_runtime_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$boundary_windows <- new.env(hash = TRUE, parent = emptyenv())
  env$transcript_space_peaks <- new.env(hash = TRUE, parent = emptyenv())
  env$heatmap_panels <- new.env(hash = TRUE, parent = emptyenv())
  env
})

site_profile_runtime_cache_key <- function(...) {
  paste(vapply(list(...), as.character, character(1)), collapse = "::")
}

site_profile_file_fingerprint <- function(path) {
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

site_profile_cached_boundary_windows <- function(kind, species, flank_size, builder) {
  key <- site_profile_runtime_cache_key(kind, species, flank_size)

  if (exists(key, envir = site_profile_runtime_cache_env$boundary_windows, inherits = FALSE)) {
    return(get(key, envir = site_profile_runtime_cache_env$boundary_windows, inherits = FALSE))
  }

  windows <- builder()
  assign(key, windows, envir = site_profile_runtime_cache_env$boundary_windows)
  windows
}

site_profile_cached_transcript_space_peaks <- function(source_file, species, cache_tag = "default", builder) {
  key <- site_profile_runtime_cache_key(cache_tag, species, site_profile_file_fingerprint(source_file))

  if (exists(key, envir = site_profile_runtime_cache_env$transcript_space_peaks, inherits = FALSE)) {
    return(get(key, envir = site_profile_runtime_cache_env$transcript_space_peaks, inherits = FALSE))
  }

  mapped_peaks <- builder()
  assign(key, mapped_peaks, envir = site_profile_runtime_cache_env$transcript_space_peaks)
  mapped_peaks
}

site_profile_cached_heatmap_panels <- function(kind, source_file, species, flank_size, builder, cache_token = NULL) {
  cache_token <- if (is.null(cache_token) || !nzchar(cache_token)) "default" else cache_token
  key <- site_profile_runtime_cache_key(kind, species, flank_size, cache_token, site_profile_file_fingerprint(source_file))

  if (exists(key, envir = site_profile_runtime_cache_env$heatmap_panels, inherits = FALSE)) {
    return(get(key, envir = site_profile_runtime_cache_env$heatmap_panels, inherits = FALSE))
  }

  panels <- builder()
  assign(key, panels, envir = site_profile_runtime_cache_env$heatmap_panels)
  panels
}

site_profile_runtime_cache_status <- function() {
  list(
    boundary_window_keys = ls(envir = site_profile_runtime_cache_env$boundary_windows),
    transcript_space_peak_keys = ls(envir = site_profile_runtime_cache_env$transcript_space_peaks),
    heatmap_panel_keys = ls(envir = site_profile_runtime_cache_env$heatmap_panels)
  )
}

site_profile_runtime_cache_reset <- function() {
  rm(
    list = ls(envir = site_profile_runtime_cache_env$boundary_windows),
    envir = site_profile_runtime_cache_env$boundary_windows
  )
  rm(
    list = ls(envir = site_profile_runtime_cache_env$transcript_space_peaks),
    envir = site_profile_runtime_cache_env$transcript_space_peaks
  )
  rm(
    list = ls(envir = site_profile_runtime_cache_env$heatmap_panels),
    envir = site_profile_runtime_cache_env$heatmap_panels
  )
  invisible(TRUE)
}
