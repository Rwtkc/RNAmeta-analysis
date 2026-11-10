structure_job_worker_bootstrap <- function(app_root_path, job_id) {
  source(file.path(app_root_path, "global.R"), local = .GlobalEnv, encoding = "UTF-8")
  source(app_path("services", "cache_cleanup.R"), local = .GlobalEnv, encoding = "UTF-8")
  source(app_path("modules", "upload_run", "upload_run.shared.R"), local = .GlobalEnv, encoding = "UTF-8")
  source(app_path("modules", "upload_run", "upload_run.functions.R"), local = .GlobalEnv, encoding = "UTF-8")
  source(app_path("modules", "browser", "browser.functions.R"), local = .GlobalEnv, encoding = "UTF-8")
  source(app_path("modules", "structure", "structure.server.R"), local = .GlobalEnv, encoding = "UTF-8")

  structure_background_worker_run(job_id)
}

mod_structure_server <- function(id, upload_context = NULL, analysis_lock = NULL) {
  moduleServer(id, function(input, output, session) {
    analysis_result <- reactiveVal(NULL)
    analysis_error <- reactiveVal(NULL)
    status_message <- reactiveVal("")
    loaded_upload_context <- reactiveVal(NULL)
    is_submitting <- reactiveVal(FALSE)
    structure_job_state <- reactiveVal(NULL)
    active_job_id <- reactiveVal("")
    queue_status <- reactiveVal(structure_job_capacity_snapshot())
    mapping_table <- reactiveVal(NULL)
    mapping_status <- reactiveVal("idle")
    mapping_message <- reactiveVal("")
    mapping_progress <- reactiveVal(list(value = 8L, label = "Validating saved BED input"))
    mapping_revision <- reactiveVal(0L)
    mapping_rendered <- reactiveVal(0L)
    mapping_page <- reactiveVal(1L)
    mapping_query <- reactiveVal("")
    submission_revision <- reactiveVal(0L)
    module_lock_id <- session$ns("analysis_lock")
    last_published_job_state <- NULL
    last_published_queue_status <- NULL

    current_upload <- reactive({
      if (is.null(upload_context)) NULL else upload_context()
    })

    current_lock_owner <- reactive({
      if (is.null(analysis_lock)) {
        return(NULL)
      }

      analysis_lock()
    })

    is_analysis_locked <- reactive({
      !is.null(current_lock_owner())
    })

    visible_mapping_status <- function() {
      status <- isolate(mapping_status())
      if (
        status %in% c("ready", "error") &&
          !identical(as.integer(isolate(mapping_rendered())), as.integer(isolate(mapping_revision())))
      ) {
        return("preparing")
      }

      status
    }

    publish_controls <- function() {
      loaded_input <- isolate(loaded_upload_context())
      saved_input <- if (is.null(loaded_input)) isolate(current_upload()) else loaded_input
      active_job <- isolate(structure_job_state())
      job_lookup_value <- isolate(input$job_lookup %||% "")
      if (!nzchar(job_lookup_value) && !is.null(active_job$jobId) && nzchar(active_job$jobId)) {
        job_lookup_value <- active_job$jobId
      }
      send_control_payload(
        session = session,
        id = session$ns("controls_host"),
        control = "structure-controls",
        config = structure_controls_config(
          id = id,
          upload_context = saved_input,
          job_lookup = job_lookup_value,
          loaded = !is.null(loaded_input),
          active_job_state = active_job,
          queue_status = isolate(queue_status()),
          mapping_status = visible_mapping_status(),
          mapping_message = isolate(mapping_message()),
          analysis_locked = isolate(is_analysis_locked())
        )
      )
    }

    publish_mapping <- function() {
      send_control_payload(
        session = session,
        id = session$ns("mapping_host"),
        control = "structure-mapping",
        config = structure_mapping_config(
          id = id,
          upload_context = isolate(current_upload()),
          table = isolate(mapping_table()),
          status = isolate(mapping_status()),
          message = isolate(mapping_message()),
          revision = isolate(mapping_revision()),
          submission_revision = isolate(submission_revision()),
          active_job_state = isolate(structure_job_state()),
          result = isolate(analysis_result()),
          page = isolate(mapping_page()),
          query = isolate(mapping_query())
        )
      )
    }

    publish_progress <- function() {
      send_control_payload(
        session = session,
        id = session$ns("progress_host"),
        control = "structure-progress",
        config = list(
          workspaceStatus = visible_mapping_status(),
          progress = isolate(mapping_progress())
        )
      )
    }

    publish_export_controls <- function() {
      if (is.null(isolate(analysis_result()))) {
        clear_control_payload(session = session, id = session$ns("export_host"))
        return(invisible(NULL))
      }

      send_control_payload(
        session = session,
        id = session$ns("export_host"),
        control = "structure-export",
        config = structure_export_config(id)
      )
    }

    publish_server_snapshot <- function(force = FALSE) {
      current_state <- isolate(structure_job_state())
      current_queue <- structure_job_capacity_snapshot()
      state_changed <- isTRUE(force) || !identical(current_state, last_published_job_state)
      queue_changed <- isTRUE(force) || !identical(current_queue, last_published_queue_status)
      if (!state_changed && !queue_changed) {
        return(invisible(NULL))
      }

      queue_status(current_queue)
      last_published_job_state <<- current_state
      last_published_queue_status <<- current_queue
      publish_controls()
      publish_mapping()
      publish_progress()
      publish_export_controls()
      invisible(NULL)
    }

    refresh_job_from_disk <- function(job_id) {
      job_bundle <- structure_job_query(job_id)
      if (is.null(job_bundle) || is.null(job_bundle$state)) {
        analysis_result(NULL)
        analysis_error(sprintf("Structure job id '%s' was not found or has expired after 7 days.", job_id))
        structure_job_state(list(
          jobId = job_id,
          status = "not found",
          engineLabel = "RNAfold",
          message = sprintf("Structure job id '%s' was not found or has expired after 7 days.", job_id)
        ))
        return(invisible(NULL))
      }

      state <- job_bundle$state
      active_job_id(state$jobId %||% job_id)
      structure_job_state(state)

      if (identical(state$status, "completed")) {
        if (is.null(job_bundle$result)) {
          analysis_result(NULL)
          analysis_error(sprintf("Structure job '%s' is marked completed but no result payload was found.", job_id))
          return(invisible(NULL))
        }

        analysis_error(NULL)
        analysis_result(job_bundle$result)
        status_message(sprintf(
          "%s | %s nt | %s | %s",
          job_bundle$result$fileName %||% "Structure job",
          job_bundle$result$sequenceLength %||% nchar(job_bundle$result$sequence %||% ""),
          job_bundle$result$locus %||% "",
          job_bundle$result$method %||% ""
        ))
        return(invisible(job_bundle$result))
      }

      analysis_result(NULL)
      if (identical(state$status, "failed")) {
        analysis_error(state$errorMessage %||% state$message %||% "Structure analysis failed.")
      } else {
        analysis_error(NULL)
      }

      invisible(state)
    }

    observeEvent(input$map_structure, {
      req(!is_analysis_locked())

      saved_input <- isolate(current_upload())

      if (!structure_has_saved_upload(saved_input)) {
        loaded_upload_context(NULL)
        mapping_table(NULL)
        mapping_status("error")
        mapping_message("Save a BED file in Upload / Run first to enable Structure analysis.")
        mapping_progress(list(value = 100L, label = "Structure mapping failed"))
        mapping_revision(mapping_revision() + 1L)
        mapping_rendered(0L)
        publish_controls()
        publish_mapping()
        publish_progress()
        status_message(mapping_message())
        return(invisible(NULL))
      }

      is_submitting(TRUE)
      if (!is.null(analysis_lock)) {
        analysis_lock(module_lock_id)
      }
      on.exit({
        is_submitting(FALSE)
        if (!is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
      }, add = TRUE)

      mapping_status("preparing")
      mapping_message("Mapping saved BED intervals to transcript coordinates.")
      mapping_progress(list(value = 8L, label = "Validating saved BED input"))
      mapping_table(NULL)
      mapping_rendered(0L)
      mapping_page(1L)
      mapping_query("")
      mapping_revision(mapping_revision() + 1L)
      publish_controls()
      publish_mapping()
      publish_progress()

      mapped <- tryCatch(
        browser_build_transcript_table(
          upload_context = saved_input,
          cache_scope = sprintf("structure:%s", session$token),
          progress = function(value, detail) {
            mapping_progress(list(value = as.integer(value), label = as.character(detail)))
            publish_progress()
          }
        ),
        error = function(err) {
          mapping_message(conditionMessage(err))
          NULL
        }
      )

      if (is.null(mapped)) {
        mapping_status("error")
        mapping_progress(list(value = 100L, label = "Structure mapping failed"))
        mapping_revision(mapping_revision() + 1L)
        publish_controls()
        publish_mapping()
        publish_progress()
        return(invisible(NULL))
      }

      if (!"source_file_index" %in% names(mapped)) {
        mapped[, source_file_index := suppressWarnings(as.integer(sub("^f([0-9]+)_.*$", "\\1", row_key)))]
      }

      loaded_upload_context(saved_input)
      mapping_table(mapped)
      mapping_status("ready")
      mapping_message("")
      mapping_progress(list(
        value = 96L,
        label = sprintf("Sending %s mapped intervals to the table", nrow(mapped))
      ))
      mapping_revision(mapping_revision() + 1L)
      publish_controls()
      publish_mapping()
      publish_progress()
    }, ignoreInit = TRUE)

    observeEvent(input$mapping_rendered, {
      rendered_revision <- suppressWarnings(as.integer(input$mapping_rendered %||% NA_integer_))
      if (!is.finite(rendered_revision) || !identical(rendered_revision, as.integer(mapping_revision()))) {
        return(invisible(NULL))
      }

      mapping_rendered(rendered_revision)
      mapping_progress(list(value = 100L, label = "Transcript mapping table ready"))
      publish_controls()
      publish_progress()
    }, ignoreInit = TRUE)

    observeEvent(input$structure_submit_mapping, {
      req(!is_analysis_locked())

      row_key <- trimws(as.character(input$structure_submit_mapping %||% ""))
      table <- isolate(mapping_table())
      match_index <- if (is.null(table)) integer() else which(as.character(table$row_key) == row_key)
      if (!length(match_index)) {
        mapping_message("The selected mapping row is no longer available. Run Structure again.")
        submission_revision(submission_revision() + 1L)
        publish_mapping()
        return(invisible(NULL))
      }

      saved_input <- isolate(current_upload())
      mapping_row <- as.list(table[match_index[[1L]]])
      analysis_result(NULL)
      analysis_error(NULL)
      is_submitting(TRUE)
      if (!is.null(analysis_lock)) {
        analysis_lock(module_lock_id)
      }
      keep_job_lock <- FALSE
      on.exit({
        is_submitting(FALSE)
        if (!keep_job_lock && !is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
        submission_revision(submission_revision() + 1L)
        publish_server_snapshot(force = TRUE)
      }, add = TRUE)

      submission <- tryCatch(
        structure_submit_background_job(
          upload_context = saved_input,
          controls = list(mapping_row = mapping_row)
        ),
        error = function(err) {
          analysis_error(conditionMessage(err))
          NULL
        }
      )

      if (is.null(submission)) {
        return(invisible(NULL))
      }

      loaded_upload_context(saved_input)
      active_job_id(submission$job_id)
      structure_job_state(submission$state)
      keep_job_lock <- TRUE
      updateTextInput(session, "job_lookup", value = submission$job_id)
      status_message(sprintf("RNAfold job %s has been queued.", submission$job_id))
    }, ignoreInit = TRUE)

    observeEvent(current_upload(), {
      loaded_upload_context(NULL)
      mapping_table(NULL)
      mapping_status("idle")
      mapping_message("")
      mapping_progress(list(value = 8L, label = "Validating saved BED input"))
      mapping_revision(mapping_revision() + 1L)
      mapping_rendered(0L)
      mapping_page(1L)
      mapping_query("")
      publish_controls()
      publish_mapping()
      publish_progress()
      status_message("Saved BED context changed. Click Run Structure to load it.")
    }, ignoreInit = TRUE)

    observeEvent(input$controls_ready, {
      publish_controls()
    }, ignoreInit = TRUE)

    observeEvent(input$mapping_ready, {
      publish_mapping()
    }, ignoreInit = TRUE)

    observeEvent(input$mapping_page, {
      next_page <- suppressWarnings(as.integer(input$mapping_page %||% 1L))
      mapping_page(max(1L, next_page))
      publish_mapping()
    }, ignoreInit = TRUE)

    observeEvent(input$mapping_query, {
      mapping_query(trimws(as.character(input$mapping_query %||% "")))
      mapping_page(1L)
      publish_mapping()
    }, ignoreInit = TRUE)

    observeEvent(input$progress_ready, {
      publish_progress()
    }, ignoreInit = TRUE)

    observeEvent(input$export_ready, {
      publish_export_controls()
    }, ignoreInit = TRUE)

    export_structure <- function() {
      session$sendCustomMessage(
        "rnameta-structure-export",
        list(
          frameId = session$ns("forna_frame"),
          format = input$export_format %||% "png",
          width = suppressWarnings(as.numeric(input$export_width %||% 3000)),
          height = suppressWarnings(as.numeric(input$export_height %||% 1800)),
          dpi = suppressWarnings(as.numeric(input$export_dpi %||% 300))
        )
      )
    }

    observeEvent(input$export_figure, {
      export_structure()
    }, ignoreInit = TRUE)

    observeEvent(input$query_job, {
      requested_id <- structure_job_normalize_id(input$job_lookup)
      if (!nzchar(requested_id)) {
        analysis_result(NULL)
        analysis_error("Enter a valid Structure Job ID, for example structure_20260728084938_525831.")
        publish_controls()
        return(invisible(NULL))
      }

      refresh_job_from_disk(requested_id)
      publish_controls()
      publish_export_controls()
    }, ignoreInit = TRUE)

    observe({
      is_analysis_locked()
      publish_controls()
    })

    observe({
      invalidateLater(1000, session)
      current_job_id <- active_job_id()
      current_state <- structure_job_state()

      if (nzchar(current_job_id) && !is.null(current_state) && current_state$status %in% c("queued", "running")) {
        refresh_job_from_disk(current_job_id)
      } else if (
        nzchar(current_job_id) &&
        !is.null(current_state) &&
        !current_state$status %in% c("queued", "running") &&
        !is.null(analysis_lock) &&
        identical(analysis_lock(), module_lock_id)
      ) {
        if (!is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
      }

      publish_server_snapshot()
    })

    observe({
      analysis_result()
      publish_export_controls()
    })

    observe({
      result <- analysis_result()
      if (is.null(result)) {
        return(invisible(NULL))
      }

      session$onFlushed(function() {
        structure_send_render(session, result)
      }, once = TRUE)
    })

    observeEvent(input$status, {
      status <- input$status
      message <- status$message %||% ""
      if (!startsWith(message, "Rendered ")) {
        status_message(message)
      }
    }, ignoreInit = FALSE)

    output$status_message <- renderText({
      status_message()
    })

    output$parameter_snapshot <- renderUI({
      saved_input <- current_upload()
      file_count <- if (is.null(saved_input)) 0 else upload_run_context_file_count(saved_input)
      species_label <- if (is.null(saved_input) || is.null(saved_input$species) || !nzchar(saved_input$species)) {
        "Not saved yet"
      } else {
        upload_run_display_species(saved_input$species)
      }
      interval_label <- if (identical(mapping_status(), "ready")) "Transcript mapping table" else "Map BED rows"
      window_label <- "Full size"

      div(
        class = "structure-summary-grid",
        div(
          class = "structure-summary-item",
          span(class = "structure-summary-item__label", "Saved BED Files"),
          span(class = "structure-summary-item__value", if (is.null(saved_input) || file_count == 0) "Not saved yet" else sprintf("%d file(s)", file_count))
        ),
        div(
          class = "structure-summary-item",
          span(class = "structure-summary-item__label", "Species"),
          span(class = "structure-summary-item__value", species_label)
        ),
        div(
          class = "structure-summary-item",
          span(class = "structure-summary-item__label", "Interval Source"),
          span(class = "structure-summary-item__value", interval_label)
        ),
        div(
          class = "structure-summary-item",
          span(class = "structure-summary-item__label", "Window Size"),
          span(class = "structure-summary-item__value", window_label)
        )
      )
    })

    output$analysis_output <- renderUI({
      state <- structure_job_state()
      result <- analysis_result()
      error_message <- analysis_error()

      if (is.null(state) && is.null(result) && (is.null(error_message) || !nzchar(error_message))) {
        return(NULL)
      }

      if (is.null(result)) {
        if (is.null(error_message) || !nzchar(error_message)) {
          return(NULL)
        }

        return(structure_job_notice_card(job_state = state, error_message = error_message))
      }

      tagList(
        structure_basic_info_card(result),
        tags$iframe(
          id = session$ns("forna_frame"),
          class = "structure-forna-frame",
          src = "structure/structure-preview.html?v=20260804a",
          title = "RNA structure preview"
        )
      )
    })

    session$onFlushed(function() {
      publish_controls()
      publish_mapping()
      publish_export_controls()
      status_message("Save a BED file in Upload / Run, then click Run Structure.")
    }, once = TRUE)

    outputOptions(output, "analysis_output", suspendWhenHidden = FALSE)

    invisible(NULL)
  })
}
