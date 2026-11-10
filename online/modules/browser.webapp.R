browser_annotation_view_spec <- function(spec, annotation_view) {
  views <- spec$annotation_views %||% list()
  annotation_view <- as.character(annotation_view %||% "transcript")

  if (!annotation_view %in% names(views)) {
    annotation_view <- "transcript"
  }

  view <- views[[annotation_view]]
  view$id <- annotation_view
  view
}

browser_build_annotation_track <- function(spec, annotation_view) {
  view <- browser_annotation_view_spec(spec, annotation_view)
  annotation_path <- browser_igv_asset_first(spec, view$files)
  if (is.null(annotation_path)) {
    return(NULL)
  }

  annotation_index <- browser_igv_asset_file(spec, view$index)
  track <- list(
    name = view$track_name,
    type = "annotation",
    format = "bed",
    url = browser_igv_asset_url(annotation_path),
    displayMode = "EXPANDED",
    height = view$height,
    color = view$color
  )

  if (!is.null(annotation_index)) {
    track$indexURL <- browser_igv_asset_url(annotation_index)
  }

  track$sourcePath <- annotation_path
  track
}

browser_build_config <- function(upload_context, annotation_view = "transcript", use_filtered_cache = FALSE, cache_scope = NULL, selected_transcript = NULL) {
  browser_register_resource_paths()

  species <- if (is.null(upload_context)) "" else upload_context$species
  spec <- browser_reference_spec(species)
  saved_files <- if (is.null(upload_context)) list() else upload_run_context_saved_files(upload_context)

  if (is.null(spec)) {
    return(list(
      status = "unsupported_species",
      message = "Save a BED file in Upload / Run before opening the genome browser.",
      species = species %||% "",
      cacheScope = browser_cache_scope(cache_scope),
      tracks = list()
    ))
  }

  annotation_view_spec <- browser_annotation_view_spec(spec, annotation_view)
  selected_track_bundle <- browser_build_selected_transcript_tracks(
    selected_transcript,
    spec = spec,
    cache_scope = cache_scope
  )
  transcript_reference <- if (!is.null(selected_track_bundle)) selected_track_bundle$transcriptReference else NULL

  tracks <- if (!is.null(selected_track_bundle)) selected_track_bundle$tracks else list()
  annotation_tracked_files <- unlist(lapply(tracks, function(track) {
    c(track$sourcePath %||% "", track$indexURL %||% "")
  }), use.names = FALSE)
  tracked_files <- c(
    if (!is.null(transcript_reference)) transcript_reference$chrom_sizes_path else "",
    annotation_tracked_files,
    vapply(saved_files, function(file) file$file_path, character(1))
  )
  tracked_files <- tracked_files[!is.na(tracked_files) & nzchar(tracked_files)]
  tracked_files <- vapply(tracked_files, function(path) {
    if (grepl("^rnameta-igv-assets/", path)) {
      return(file.path(app_path("data", "igv"), sub("^rnameta-igv-assets/", "", path)))
    }
    path
  }, character(1))
  tracked_files <- tracked_files[file.exists(tracked_files)]
  tracked_mtimes <- if (length(tracked_files)) {
    vapply(file.info(tracked_files)$mtime, as.numeric, numeric(1))
  } else {
    numeric()
  }

  list(
    status = if (length(tracks)) "ready" else "reference_only",
    message = if (length(tracks)) {
      "Loaded the selected transcript track for the mapped BED row."
    } else {
      "Reference and annotation are available. Map a BED row to a transcript before opening IGV."
    },
    species = spec$id,
    speciesLabel = spec$name,
    annotationView = annotation_view_spec$id,
    annotationViewLabel = annotation_view_spec$label,
    cacheScope = browser_cache_scope(cache_scope),
    fileCount = length(saved_files),
    locus = if (!is.null(transcript_reference)) transcript_reference$locus else selected_transcript$igv_locus %||% "",
    reference = if (!is.null(transcript_reference)) {
      list(
        id = transcript_reference$transcript_id,
        name = transcript_reference$transcript_id,
        format = "chromsizes",
        url = browser_public_url(transcript_reference$chrom_sizes_path),
        chromosomeOrder = list(transcript_reference$transcript_id)
      )
    } else {
      NULL
    },
    tracks = tracks,
    assetKey = paste(
      spec$id,
      length(saved_files),
      annotation_view_spec$id,
      selected_transcript$row_key %||% "",
      selected_transcript$transcript_id %||% "",
      selected_transcript$transcript_interval %||% "",
      paste(vapply(saved_files, function(file) file$file_name, character(1)), collapse = "|"),
      paste(tracked_mtimes, collapse = "|"),
      sep = "::"
    )
  )
}

browser_webapp_session_path <- function(config) {
  session_dir <- browser_public_subdir("sessions", config$cacheScope %||% NULL)

  raw_id <- paste(config$assetKey %||% as.numeric(Sys.time()), collapse = "")
  session_id <- substr(gsub("[^A-Za-z0-9_.-]", "_", raw_id), 1, 96)
  file.path(session_dir, sprintf("igv_webapp_%s.json", session_id))
}

browser_build_webapp_session <- function(config) {
  reference <- config$reference
  tracks <- config$tracks %||% list()

  if (is.null(reference) || !length(tracks)) {
    return(NULL)
  }

  reference$url <- browser_public_absolute_url(reference$url)
  if (!is.null(reference$aliasURL)) {
    reference$aliasURL <- browser_public_absolute_url(reference$aliasURL)
  }

  tracks <- lapply(tracks, function(track) {
    track$url <- browser_public_absolute_url(track$url)
    if (!is.null(track$indexURL)) {
      track$indexURL <- browser_public_absolute_url(track$indexURL)
    }
    track$sourcePath <- NULL
    track$transcriptReference <- NULL
    track
  })

  list(
    reference = reference,
    locus = config$locus %||% if (length(reference$chromosomeOrder)) reference$chromosomeOrder[[1]] else "1",
    loadDefaultGenomes = FALSE,
    showNavigation = TRUE,
    showRuler = TRUE,
    showSequence = FALSE,
    showSVGButton = FALSE,
    tracks = tracks
  )
}

browser_write_webapp_session <- function(config) {
  session_object <- browser_build_webapp_session(config)

  if (is.null(session_object)) {
    return(NULL)
  }

  session_path <- browser_webapp_session_path(config)
  jsonlite::write_json(
    session_object,
    session_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  session_path
}

browser_webapp_iframe_src <- function(config) {
  if (!browser_webapp_ready()) {
    return(NULL)
  }

  session_path <- browser_write_webapp_session(config)
  if (is.null(session_path)) {
    return(NULL)
  }

  session_url <- browser_public_absolute_url(browser_public_url(session_path))
  paste0(
    "rnameta-igv-webapp/index.html?build=20260803-locus-lock&sessionURL=",
    utils::URLencode(session_url, reserved = TRUE)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
