browser_public_dir <- function() {
  path <- app_path("cache", "igv_browser")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

browser_cache_scope <- function(value = NULL) {
  value <- as.character(value %||% "shared")[[1]]
  value <- gsub("[^A-Za-z0-9_.-]+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  if (!nzchar(value)) {
    value <- "shared"
  }

  substr(value, 1, 96)
}

browser_public_subdir <- function(kind, cache_scope = NULL) {
  path <- file.path(browser_public_dir(), kind, browser_cache_scope(cache_scope))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

browser_cleanup_cache_scope <- function(cache_scope) {
  scope <- browser_cache_scope(cache_scope)
  if (identical(scope, "shared")) {
    return(invisible(0L))
  }

  paths <- file.path(browser_public_dir(), c("uploads", "filtered", "sessions", "selected_transcripts"), scope)
  if (exists("rnameta_unlink_cache_paths", mode = "function")) {
    return(rnameta_unlink_cache_paths(paths, recursive = TRUE))
  }

  unlink(paths, recursive = TRUE, force = TRUE)
  invisible(sum(!file.exists(paths)))
}

browser_register_resource_paths <- function() {
  fa_dir <- app_path("data", "fa")
  igv_dir <- app_path("data", "igv")
  webapp_dir <- app_path("external", "igv-webapp", "dist")
  public_dir <- browser_public_dir()

  if (dir.exists(fa_dir) && !"rnameta-fa" %in% names(shiny::resourcePaths())) {
    shiny::addResourcePath("rnameta-fa", fa_dir)
  }

  if (dir.exists(igv_dir) && !"rnameta-igv-assets" %in% names(shiny::resourcePaths())) {
    shiny::addResourcePath("rnameta-igv-assets", igv_dir)
  }

  if (!"rnameta-igv" %in% names(shiny::resourcePaths())) {
    shiny::addResourcePath("rnameta-igv", public_dir)
  }

  if (dir.exists(webapp_dir) && !"rnameta-igv-webapp" %in% names(shiny::resourcePaths())) {
    shiny::addResourcePath("rnameta-igv-webapp", webapp_dir)
  }
}

browser_progress_start <- function(session) {
  progress <- shiny::Progress$new(session = session, min = 0, max = 100)
  progress$set(
    value = 0,
    message = "Preparing Genome Browser",
    detail = "0% | Initializing browser workspace"
  )
  progress
}

browser_progress_update <- function(progress, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  safe_value <- as.integer(round(max(0, min(100, value))))
  progress$set(
    value = safe_value,
    message = "Preparing Genome Browser",
    detail = sprintf("%s%% | %s", safe_value, detail)
  )

  invisible(NULL)
}

browser_progress_close <- function(progress) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  Sys.sleep(0.6)
  progress$close()
  invisible(NULL)
}
