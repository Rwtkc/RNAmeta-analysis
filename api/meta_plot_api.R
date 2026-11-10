#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

library(httpuv)
library(jsonlite)

`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0L) fallback else value
}

app_root <- Sys.getenv("RNAMETA_ONLINE_ROOT", "/opt/rnameta-online")
port <- suppressWarnings(as.integer(Sys.getenv("METAPLOT_API_PORT", "10812")))
if (!is.finite(port) || port < 1L || port > 65535L) port <- 10812L
if (!dir.exists(app_root)) stop(sprintf("RNAmeta-online root is missing: %s", app_root))

setwd(app_root)
source(file.path(app_root, "global.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "upload_run", "upload_run.shared.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "peak_distribution", "peak_distribution.data.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "services", "meta_plot_runtime_cache.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "meta_plot", "meta_plot.progress.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "meta_plot", "meta_plot.data.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "meta_plot", "meta_plot.annotation.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "meta_plot", "meta_plot.mapping.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "meta_plot", "meta_plot.plot.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "meta_plot", "meta_plot.chart.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "meta_plot", "meta_plot.export.R"), local = .GlobalEnv, encoding = "UTF-8")
source(file.path(app_root, "modules", "meta_plot", "meta_plot.analysis.R"), local = .GlobalEnv, encoding = "UTF-8")

upload_run_context_saved_files <- function(upload_context) upload_context$saved_files %||% list()

meta_plot_api_defaults <- list(
  pltTxType = "tx", headOrtail = "TRUE", txPrimaryOnly = "TRUE",
  txlncrnaOverlapmrna = "TRUE", mapFilterTranscript = "TRUE", enableCI = "FALSE",
  stSampleModle = "Equidistance", stSampleNum = 10L, stAmblguity = 5L,
  txfiveutrMinLength = 100L, txcdsMinLength = 100L, txthreeutrMinLength = 100L,
  txlongNcrnaMinLength = 100L, txpromoterLength = 1000L, txtailLength = 1000L,
  adjust = 1, overlapIndex = 1, siteLengthIndex = 1
)

json_response <- function(status, payload) {
  list(
    status = status,
    headers = list(
      "Content-Type" = "application/json; charset=utf-8",
      "Access-Control-Allow-Origin" = "*",
      "Access-Control-Allow-Headers" = "Content-Type",
      "Access-Control-Allow-Methods" = "GET, POST, OPTIONS"
    ),
    body = jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", dataframe = "rows", digits = 12)
  )
}

text_response <- function(status, body, content_type = "text/plain; charset=utf-8") {
  list(
    status = status,
    headers = list(
      "Content-Type" = content_type,
      "Access-Control-Allow-Origin" = "*",
      "Access-Control-Allow-Headers" = "Content-Type",
      "Access-Control-Allow-Methods" = "GET, POST, OPTIONS"
    ),
    body = body
  )
}

query_parameters <- function(request) {
  raw_path <- request$PATH_INFO %||% ""
  path_query <- if (grepl("\\?", raw_path)) sub("^[^?]*\\?", "", raw_path) else ""
  query <- request$QUERY_STRING %||% path_query
  if (!nzchar(query)) query <- path_query
  if (!nzchar(query)) return(list())
  parts <- strsplit(query, "&", fixed = TRUE)[[1]]
  values <- lapply(parts, function(part) {
    pair <- strsplit(part, "=", fixed = TRUE)[[1]]
    value <- if (length(pair) > 1L) paste(pair[-1L], collapse = "=") else ""
    utils::URLdecode(gsub("\\+", " ", value, fixed = FALSE))
  })
  keys <- vapply(parts, function(part) strsplit(part, "=", fixed = TRUE)[[1]][[1]], character(1))
  stats::setNames(values, keys)
}

database_catalog_payload <- function() {
  upload_run_database_control_catalog(upload_run_database_catalog(), include_files = TRUE)
}

database_file_request <- function(request) {
  payload <- request_payload(request)
  species <- single_string(payload$species)
  source_id <- single_string(payload$source_id)
  file_key <- single_string(payload$file_key)
  if (!nzchar(species) || !nzchar(source_id) || !nzchar(file_key)) {
    stop("Species, database, and BED file are required.", call. = FALSE)
  }
  selected <- upload_run_database_selected_files(species, source_id, file_key)
  path <- selected$files[[1]]$file_path
  content <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (nzchar(content)) content <- paste0(content, "\n")
  text_response(200L, content)
}

request_payload <- function(request) {
  body <- request$rook.input$read() %||% ""
  if (is.raw(body)) body <- rawToChar(body)
  if (length(body) > 1L) body <- paste(body, collapse = "")
  if (!nzchar(body)) stop("Request body is empty.", call. = FALSE)
  payload <- jsonlite::fromJSON(body, simplifyVector = FALSE)
  if (!is.list(payload)) stop("Request body must be a JSON object.", call. = FALSE)
  payload
}

single_string <- function(value, fallback = "") {
  if (is.null(value) || !length(value)) return(fallback)
  as.character(value)[[1]]
}

single_number <- function(value, fallback, minimum = -Inf, maximum = Inf, integer = FALSE) {
  parsed <- suppressWarnings(as.numeric(value %||% fallback))
  if (!length(parsed) || !is.finite(parsed[[1]]) || parsed[[1]] < minimum || parsed[[1]] > maximum) return(fallback)
  if (integer) as.integer(round(parsed[[1]])) else parsed[[1]]
}

normalize_controls <- function(input) {
  controls <- utils::modifyList(meta_plot_api_defaults, input %||% list())
  controls$pltTxType <- if (single_string(controls$pltTxType) %in% c("tx", "mrna", "ncrna")) single_string(controls$pltTxType) else "tx"
  controls$stSampleModle <- if (single_string(controls$stSampleModle) %in% c("Equidistance", "random")) single_string(controls$stSampleModle) else "Equidistance"
  for (name in c("headOrtail", "txPrimaryOnly", "txlncrnaOverlapmrna", "mapFilterTranscript", "enableCI")) {
    controls[[name]] <- if (identical(toupper(single_string(controls[[name]])), "TRUE")) "TRUE" else "FALSE"
  }
  controls$stSampleNum <- single_number(controls$stSampleNum, 10L, 1, 100, TRUE)
  controls$stAmblguity <- single_number(controls$stAmblguity, 5L, 0, 1000, TRUE)
  for (name in c("txfiveutrMinLength", "txcdsMinLength", "txthreeutrMinLength", "txlongNcrnaMinLength")) {
    controls[[name]] <- single_number(controls[[name]], 100L, 0, 10000000, TRUE)
  }
  controls$txpromoterLength <- single_number(controls$txpromoterLength, 1000L, 0, 10000000, TRUE)
  controls$txtailLength <- single_number(controls$txtailLength, 1000L, 0, 10000000, TRUE)
  controls$adjust <- single_number(controls$adjust, 1, 0.01, 100)
  controls$overlapIndex <- single_number(controls$overlapIndex, 1, 0, 10)
  controls$siteLengthIndex <- single_number(controls$siteLengthIndex, 1, 0, 10)
  controls
}

write_request_files <- function(files) {
  if (!is.list(files) || length(files) < 1L || length(files) > 6L) stop("Meta Plot accepts 1 to 6 BED files.", call. = FALSE)
  lapply(seq_along(files), function(index) {
    entry <- files[[index]]
    name <- basename(single_string(entry$name, sprintf("sample-%d.bed", index)))
    content <- single_string(entry$content)
    if (!nzchar(content)) stop(sprintf("%s is empty.", name), call. = FALSE)
    path <- tempfile(pattern = sprintf("meta_plot_input_%d_", index), fileext = ".bed")
    writeLines(content, con = path, useBytes = TRUE)
    list(file_name = name, file_path = path)
  })
}

species_payload <- function() {
  choices <- read_species_choices()
  lapply(choices, function(choice) list(value = choice$value, label = choice$label))
}

run_request <- function(payload) {
  species <- single_string(payload$species)
  available <- vapply(species_payload(), function(choice) choice$value, character(1))
  if (!nzchar(species) || !(species %in% available)) stop("Choose a supported species.", call. = FALSE)
  saved_files <- write_request_files(payload$files)
  source_paths <- vapply(saved_files, function(file) file$file_path, character(1))
  prepared_paths <- character()
  on.exit(unlink(unique(c(source_paths, prepared_paths)), force = TRUE), add = TRUE)
  result <- run_meta_plot_analysis(
    upload_context = list(species = species, saved_files = saved_files),
    controls = normalize_controls(payload$controls),
    progress = NULL
  )
  prepared_paths <- result$prepared_bed_files %||% character()
  list(
    chart = result$chart_payload,
    summary = list(
      species = result$species,
      sample_count = result$sample_count,
      interval_count = result$interval_count,
      sampled_point_count = result$sampled_point_count,
      transcript_hits = result$overlap_count
    ),
    data = meta_plot_data_export_table(result)
  )
}

handle_request <- function(request) {
  method <- toupper(request$REQUEST_METHOD %||% "GET")
  path <- sub("\\?.*$", "", request$PATH_INFO %||% "/")
  if (identical(method, "OPTIONS")) return(json_response(204L, list()))
  if (identical(method, "GET") && identical(path, "/health")) return(json_response(200L, list(status = "ok", service = "meta-plot-api")))
  if (identical(method, "GET") && identical(path, "/species")) return(json_response(200L, list(species = species_payload())))
  if (identical(method, "GET") && identical(path, "/database/catalog")) return(json_response(200L, list(databases = database_catalog_payload())))
  if (identical(method, "POST") && identical(path, "/database/file")) return(database_file_request(request))
  if (!identical(method, "POST") || !identical(path, "/run")) return(json_response(404L, list(error = "Not found.")))
  tryCatch(
    json_response(200L, run_request(request_payload(request))),
    error = function(error) {
      message(sprintf("%s | run | %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), conditionMessage(error)))
      json_response(400L, list(error = conditionMessage(error)))
    }
  )
}

app <- list(call = function(request) {
  tryCatch(handle_request(request), error = function(error) {
    message(sprintf("%s | request | %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), conditionMessage(error)))
    json_response(500L, list(error = "Meta Plot request failed."))
  })
})

message(sprintf("Meta Plot API listening on 0.0.0.0:%d", port))
httpuv::runServer("0.0.0.0", port, app)
