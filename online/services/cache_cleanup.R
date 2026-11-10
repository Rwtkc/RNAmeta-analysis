rnameta_cache_cleanup_state <- new.env(parent = emptyenv())

rnameta_cache_root <- function() {
  path <- app_path("cache")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

rnameta_cache_path <- function(...) {
  file.path(rnameta_cache_root(), ...)
}

rnameta_is_cache_path <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(FALSE)
  }

  cache_root <- rnameta_cache_root()
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  identical(resolved, cache_root) || startsWith(resolved, paste0(cache_root, "/"))
}

rnameta_unlink_cache_path <- function(path, recursive = FALSE) {
  if (!rnameta_is_cache_path(path)) {
    return(invisible(FALSE))
  }

  unlink(path, recursive = recursive, force = TRUE)
  invisible(!file.exists(path))
}

rnameta_unlink_cache_paths <- function(paths, recursive = FALSE) {
  if (is.null(paths)) {
    return(invisible(0L))
  }

  paths <- unique(as.character(paths))
  paths <- paths[nzchar(paths)]
  if (!length(paths)) {
    return(invisible(0L))
  }

  removed <- 0L
  for (path in paths) {
    if (isTRUE(rnameta_unlink_cache_path(path, recursive = recursive))) {
      removed <- removed + 1L
    }
  }

  invisible(removed)
}

rnameta_is_temp_path <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(FALSE)
  }

  temp_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  identical(resolved, temp_root) || startsWith(resolved, paste0(temp_root, "/"))
}

rnameta_unlink_temp_paths <- function(paths, recursive = FALSE) {
  if (is.null(paths)) {
    return(invisible(0L))
  }

  paths <- unique(as.character(paths))
  paths <- paths[nzchar(paths)]
  if (!length(paths)) {
    return(invisible(0L))
  }

  removed <- 0L
  for (path in paths) {
    if (rnameta_is_temp_path(path)) {
      unlink(path, recursive = recursive, force = TRUE)
      if (!file.exists(path)) {
        removed <- removed + 1L
      }
    }
  }

  invisible(removed)
}

rnameta_prune_cache_children <- function(path, max_age_hours, recursive = TRUE) {
  if (is.null(path) || !dir.exists(path) || !rnameta_is_cache_path(path)) {
    return(0L)
  }

  cutoff <- Sys.time() - as.difftime(max_age_hours, units = "hours")
  children <- list.files(path, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (!length(children)) {
    return(0L)
  }

  removed <- 0L
  for (child in children) {
    info <- file.info(child)
    if (!is.na(info$mtime[[1]]) && info$mtime[[1]] < cutoff) {
      if (isTRUE(rnameta_unlink_cache_path(child, recursive = recursive))) {
        removed <- removed + 1L
      }
    }
  }

  removed
}

rnameta_cache_cleanup <- function() {
  removed <- c(
    motif_runtime_jobs = rnameta_prune_cache_children(rnameta_cache_path("motif_runtime_jobs"), 168),
    motif_jobs = rnameta_prune_cache_children(rnameta_cache_path("motif_jobs"), 168),
    structure_runtime_jobs = rnameta_prune_cache_children(rnameta_cache_path("structure_runtime_jobs"), 168)
  )

  invisible(removed)
}

rnameta_cache_cleanup_due <- function(interval_hours = 1) {
  now <- Sys.time()
  last_run <- get0("last_run", envir = rnameta_cache_cleanup_state, ifnotfound = NULL)
  if (!is.null(last_run) && difftime(now, last_run, units = "hours") < interval_hours) {
    return(invisible(FALSE))
  }

  assign("last_run", now, envir = rnameta_cache_cleanup_state)
  tryCatch(
    {
      rnameta_cache_cleanup()
      invisible(TRUE)
    },
    error = function(err) {
      message("RNAmeta cache cleanup skipped: ", conditionMessage(err))
      invisible(FALSE)
    }
  )
}
