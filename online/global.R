library(shiny)
library(htmltools)
library(jsonlite)

registerInputHandler(
  "r9::shiny.number",
  function(value, shinysession, name) {
    if (is.null(value) || identical(value, "")) {
      return(NULL)
    }

    suppressWarnings(as.numeric(value))
  },
  force = TRUE
)

options(
  shiny.autoreload = identical(Sys.getenv("RNAMETA_SHINY_AUTORELOAD", unset = "false"), "true"),
  shiny.fullstacktrace = TRUE,
  shiny.maxRequestSize = 200 * 1024^2
)

resolve_app_root <- function() {
  candidates <- c(
    normalizePath(".", winslash = "/", mustWork = TRUE),
    normalizePath("./new", winslash = "/", mustWork = FALSE)
  )

  matches <- candidates[vapply(
    candidates,
    function(path) {
      dir.exists(file.path(path, "modules", "welcome"))
    },
    logical(1)
  )]

  if (length(matches) == 0) {
    stop("Unable to locate the refactored Shiny app root.")
  }

  matches[[1]]
}

app_root <- resolve_app_root()

app_path <- function(...) {
  file.path(app_root, ...)
}
