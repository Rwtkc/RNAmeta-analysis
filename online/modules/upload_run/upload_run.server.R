bed_column_names <- function(n_cols) {
  standard_names <- c(
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts"
  )

  if (n_cols <= length(standard_names)) {
    return(standard_names[seq_len(n_cols)])
  }

  c(
    standard_names,
    paste0("extra_", seq_len(n_cols - length(standard_names)))
  )
}

get_demo_bed_file <- function() {
  candidates <- c(
    app_path("data", "exampleDataFile", "ac_diff_down.bed6.bed"),
    app_path("..", "data", "exampleDataFile", "ac_diff_down.bed6.bed")
  )

  matches <- candidates[file.exists(candidates)]
  if (length(matches)) {
    return(matches[[1]])
  }

  candidates[[1]]
}

read_bed_preview <- function(file_path, n = 10, skip = 0L) {
  if (is.null(file_path) || !nzchar(file_path) || !file.exists(file_path)) {
    return(NULL)
  }

  preview_tbl <- tryCatch(
    read.delim(
      file_path,
      header = FALSE,
      sep = "\t",
      fill = TRUE,
      stringsAsFactors = FALSE,
      quote = "",
      nrows = n,
      skip = skip
    ),
    error = function(...) NULL
  )

  if (is.null(preview_tbl) || nrow(preview_tbl) == 0) {
    return(NULL)
  }

  names(preview_tbl) <- bed_column_names(ncol(preview_tbl))
  preview_tbl
}

render_preview_table <- function(preview_tbl) {
  if (is.null(preview_tbl) || nrow(preview_tbl) == 0) {
    return(tags$p(class = "upload-summary__empty", "No preview available for the saved file."))
  }

  header <- tags$tr(lapply(names(preview_tbl), function(col) tags$th(col)))
  rows <- lapply(seq_len(nrow(preview_tbl)), function(i) {
    tags$tr(lapply(preview_tbl[i, , drop = TRUE], function(value) tags$td(as.character(value))))
  })

  tags$div(
    class = "upload-preview-table-wrap",
    tags$table(
      class = "upload-preview-table",
      tags$thead(header),
      tags$tbody(rows)
    )
  )
}

mod_upload_run_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    saved_upload <- reactiveVal(NULL)
    selected_preview_index <- reactiveVal(1L)
    session_temp_files <- character()
    database_preview <- reactiveVal(NULL)
    database_preview_row_counts <- new.env(parent = emptyenv())
    database_selected <- reactiveVal(FALSE)
    database_catalog <- upload_run_database_catalog()

    track_session_temp_files <- function(saved_files) {
      paths <- vapply(saved_files, function(file) file$file_path %||% "", character(1))
      paths <- paths[nzchar(paths)]
      if (length(paths)) {
        session_temp_files <<- unique(c(session_temp_files, paths))
      }
    }

    session$onSessionEnded(function() {
      if (exists("rnameta_unlink_temp_paths", mode = "function")) {
        rnameta_unlink_temp_paths(session_temp_files)
      } else {
        unlink(session_temp_files, force = TRUE)
      }

      invisible(NULL)
    })

    publish_database_browser <- function(browser_catalog = NULL) {
      config <- list(
        database_preview_input = session$ns("database_preview_request"),
        database_preview = database_preview()
      )
      if (!is.null(browser_catalog)) {
        config$database_catalog <- browser_catalog
      }

      send_control_payload(
        session = session,
        id = session$ns("database_browser"),
        control = "upload-database-browser",
        config = config
      )
    }

    publish_upload_control <- function() {
      uploaded_files <- input$bed_file
      uploaded_count <- if (is.null(uploaded_files)) 0 else nrow(uploaded_files)
      config <- utils::modifyList(
        upload_run_select_config(id, database_catalog = database_catalog),
        list(
          data_source_value = if (is.null(input$data_source)) "Upload File" else input$data_source,
          species_value = if (is.null(input$species)) "" else input$species,
          uploaded_count = uploaded_count,
          uploaded_file_names = if (is.null(uploaded_files)) character() else as.character(uploaded_files$name)
        )
      )

      send_control_payload(
        session = session,
        id = session$ns("upload_controls"),
        control = "upload-sidebar",
        config = config
      )
      publish_database_browser()
    }

    publish_database_preview <- function() {
      publish_database_browser()
    }

    clear_saved_input_context <- function() {
      saved_upload(NULL)
      selected_preview_index(1L)
      database_preview(NULL)
    }

    session$onFlushed(function() {
      isolate(publish_upload_control())
    }, once = TRUE)

    observeEvent(input$controls_ready, {
      publish_upload_control()
    }, ignoreInit = TRUE)

    observe({
      input$data_source
      input$species
      input$bed_file
      publish_upload_control()
    })

    observeEvent(list(input$data_source, input$species), {
      clear_saved_input_context()
      publish_database_browser()
    }, ignoreInit = TRUE)

    observeEvent(input$database_context, {
      context <- input$database_context
      has_selection <- is.list(context) &&
        nzchar(as.character(context$species %||% "")) &&
        nzchar(as.character(context$source_id %||% ""))
      database_selected(has_selection)
      clear_saved_input_context()

      if (!has_selection) {
        publish_database_browser()
        return()
      }

      selected_source <- Filter(
        function(source) {
          identical(source$species, as.character(context$species)) &&
            identical(source$id, as.character(context$source_id))
        },
        database_catalog
      )
      if (length(selected_source)) {
        publish_database_browser(
          upload_run_database_control_catalog(selected_source, include_files = TRUE)
        )
      }
    }, ignoreInit = TRUE)

    # Context changes must clear first when the Demo click updates all inputs together.
    observeEvent(input$load_demo, {
      demo_file <- get_demo_bed_file()
      req(file.exists(demo_file))

      normalized_demo_file <- upload_run_normalize_bed_file(demo_file)
      preview_tbl <- read_bed_preview(normalized_demo_file, n = 10)
      saved_files <- list(
        list(
          file_name = basename(demo_file),
          file_path = normalized_demo_file
        )
      )
      track_session_temp_files(saved_files)

      saved_upload(
        upload_run_build_saved_context(
          data_source = "Upload File",
          species = "ara_TAIR10",
          saved_files = saved_files,
          preview = preview_tbl,
          is_demo = TRUE
        )
      )
      selected_preview_index(1L)
    }, priority = -10)

    observeEvent(input$save_upload, {
      req(nzchar(input$species))
      req(!is.null(input$bed_file), nrow(input$bed_file) > 0)

      saved_files <- upload_run_saved_files(input$bed_file)
      track_session_temp_files(saved_files)
      preview_tbl <- read_bed_preview(saved_files[[1]]$file_path, n = 10)

      saved_upload(
        upload_run_build_saved_context(
          data_source = if (is.null(input$data_source)) "Upload File" else input$data_source,
          species = input$species,
          saved_files = saved_files,
          preview = preview_tbl,
          is_demo = FALSE
        )
      )
      selected_preview_index(1L)
    })

    observeEvent(input$save_database_selection, {
      selection <- input$save_database_selection
      if (!is.list(selection)) {
        showNotification("Select a database BED file before saving.", type = "error")
        return()
      }

      species <- as.character(selection$species %||% "")
      source_id <- as.character(selection$source_id %||% "")
      file_keys <- selection$file_keys %||% character()
      if (!nzchar(species) || !nzchar(source_id)) {
        showNotification("Select a species and database before saving.", type = "error")
        return()
      }

      database_selection <- tryCatch(
        upload_run_database_selected_files(
          species = species,
          source_id = source_id,
          file_keys = file_keys,
          catalog = database_catalog
        ),
        error = function(error) {
          showNotification(conditionMessage(error), type = "error")
          NULL
        }
      )
      if (is.null(database_selection)) {
        return()
      }

      include_uploaded <- isTRUE(selection$include_uploaded)
      uploaded_files <- if (include_uploaded) upload_run_saved_files(input$bed_file) else list()
      if (include_uploaded && !length(uploaded_files)) {
        showNotification("Upload at least one BED file before saving the combined input.", type = "error")
        return()
      }

      if (length(uploaded_files)) {
        track_session_temp_files(uploaded_files)
      }

      saved_files <- c(uploaded_files, database_selection$files)
      preview_tbl <- read_bed_preview(saved_files[[1]]$file_path, n = 10)
      saved_upload(
        upload_run_build_saved_context(
          data_source = if (include_uploaded) "Upload + Database" else database_selection$source_label,
          species = species,
          saved_files = saved_files,
          preview = preview_tbl,
          is_demo = FALSE
        )
      )
      selected_preview_index(1L)
    })

    observeEvent(input$database_preview_request, {
      request <- input$database_preview_request
      if (!is.list(request)) {
        return()
      }

      preview <- tryCatch(
        upload_run_database_bed_preview_page(
          species = as.character(request$species %||% ""),
          source_id = as.character(request$source_id %||% ""),
          file_key = as.character(request$file_key %||% ""),
          page = request$page %||% 1L,
          page_size = 10L,
          row_count_cache = database_preview_row_counts,
          catalog = database_catalog
        ),
        error = function(error) {
          showNotification(conditionMessage(error), type = "error")
          NULL
        }
      )
      if (!is.null(preview)) {
        database_preview(preview)
        publish_database_preview()
      }
    }, ignoreInit = TRUE)

    observeEvent(input$preview_file_index, {
      saved_state <- saved_upload()
      req(!is.null(saved_state), saved_state$file_count > 0)
      next_index <- suppressWarnings(as.integer(input$preview_file_index))
      req(is.finite(next_index))
      selected_preview_index(min(max(next_index, 1L), saved_state$file_count))
    }, ignoreInit = TRUE)

    output$session_summary <- renderUI({
      saved_state <- saved_upload()
      uploaded_files <- input$bed_file
      source_label <- if (!is.null(saved_state)) {
        saved_state$data_source
      } else if (is.null(input$data_source)) {
        "Upload File"
      } else {
        input$data_source
      }
      species_label <- if (!is.null(saved_state)) {
        upload_run_display_species(saved_state$species)
      } else if (!is.null(input$species) && nzchar(input$species)) {
        upload_run_display_species(input$species)
      } else {
        "Not selected"
      }
      file_count <- if (!is.null(saved_state)) {
        saved_state$file_count
      } else if (is.null(uploaded_files)) {
        0
      } else {
        nrow(uploaded_files)
      }
      mixed_upload_count <- if (!is.null(saved_state)) {
        sum(vapply(saved_state$saved_files, function(file) identical(file$source %||% "", "Upload"), logical(1)))
      } else {
        0L
      }
      mixed_database_count <- if (!is.null(saved_state)) {
        sum(vapply(saved_state$saved_files, function(file) identical(file$source %||% "", "Database"), logical(1)))
      } else {
        0L
      }

      div(
        class = "upload-summary",
        tags$h3(class = "upload-summary__title", "Current input context"),
        div(
          class = "upload-summary__grid",
          div(
            class = "upload-summary__item",
            span(class = "upload-summary__label", "Data Source"),
            span(class = "upload-summary__value", source_label)
          ),
          div(
            class = "upload-summary__item",
            span(class = "upload-summary__label", "Species"),
            span(class = "upload-summary__value", species_label)
          ),
          div(
            class = "upload-summary__item",
            span(class = "upload-summary__label", "BED Files"),
            span(
              class = "upload-summary__value",
              if (file_count == 0) "No file uploaded" else sprintf("%d file(s)", file_count)
            )
          )
        ),
        if (identical(saved_state$data_source %||% "", "Upload + Database")) {
          tags$div(
            class = "upload-summary__source-breakdown",
            tags$span(class = "upload-summary__source-item", sprintf("Your upload: %d file(s)", mixed_upload_count)),
            tags$span(class = "upload-summary__source-item", sprintf("Database: %d file(s)", mixed_database_count))
          )
        },
        if (!is.null(saved_state)) {
          saved_file_labels <- upload_run_saved_file_labels(saved_state)
          preview_index <- min(max(selected_preview_index(), 1L), length(saved_file_labels))
          preview_file <- upload_run_selected_file(saved_state, preview_index)
          preview_tbl <- read_bed_preview(preview_file$file_path, n = 10)

          tags$div(
            class = "upload-summary__saved",
            tags$div(
              class = "upload-summary__preview-head",
              tags$h4(class = "upload-summary__preview-title", "Saved file preview (first 10 rows)"),
              tags$div(
                class = "upload-summary__saved-head",
                tags$span(class = "upload-summary__label", "Saved File"),
                tags$span(
                  class = "upload-summary__saved-name",
                  preview_file$file_name
                )
              )
            ),
            if (length(saved_file_labels) > 1) {
              tags$div(
                class = "upload-summary__saved-files",
                lapply(seq_along(saved_file_labels), function(index) {
                  tags$button(
                    type = "button",
                    class = paste(
                      "upload-summary__saved-file-chip",
                      if (index == preview_index) "is-active" else NULL
                    ),
                    onclick = sprintf(
                      "Shiny.setInputValue('%s', %d, {priority: 'event'})",
                      session$ns("preview_file_index"),
                      index
                    ),
                    saved_file_labels[[index]]
                  )
                })
              )
            },
            render_preview_table(preview_tbl)
          )
        }
      )
    })

    output$workflow_notes <- renderUI({
      if (!is.null(saved_upload()) || isTRUE(database_selected())) {
        return(NULL)
      }

      div(
        class = "upload-panel upload-panel--guide",
        div(class = "upload-panel__eyebrow", "Workflow Notes"),
        tags$h3(class = "upload-guide__title", "BED-first preparation for modular RNAmeta analysis"),
        tags$ul(
          class = "upload-guide__list",
          tags$li("Match the uploaded BED file to the correct species before running any downstream module."),
          tags$li("RNAmeta accepts BED-formatted chromosomal coordinates generated by common peak callers."),
          tags$li("Uploaded files from this page will become the shared input context for Peak Distribution, Motif, and feature statistics.")
        )
      )
    })

    list(
      upload_context = reactive(saved_upload())
    )
  })
}
