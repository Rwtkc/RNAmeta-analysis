`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0L) fallback else value
}

app_path <- function(...) {
  file.path(getOption("rnameta.reference_root", "."), ...)
}
