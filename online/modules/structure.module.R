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
    module_lock_id <- session$ns("analysis_lock")

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

    has_structure_coordinate <- reactive({
      nzchar(trimws(as.character(input$structure_coordinate %||% "")))
    })

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
          selected_sample = isolate(input$structure_sample %||% 1L),
          coordinate = isolate(input$structure_coordinate %||% structure_controls_defaults()$coordinate),
      window_size = structure_controls_defaults()$window_size,
          job_lookup = job_lookup_value,
          loaded = !is.null(loaded_input),
          active_job_state = active_job,
          queue_status = isolate(queue_status())
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

    observeEvent(input$run_structure, {
      req(!is_analysis_locked())

      saved_input <- isolate(current_upload())
      coordinate_value <- trimws(as.character(isolate(input$structure_coordinate %||% "")))

      if (!structure_has_saved_upload(saved_input)) {
        loaded_upload_context(NULL)
        publish_controls()
        status_message("Save a BED file in Upload / Run before running Structure.")
        return(invisible(NULL))
      }

      if (!nzchar(coordinate_value)) {
        publish_controls()
        status_message("Enter a coordinate like chr1:72338-74096 before running Structure.")
        return(invisible(NULL))
      }

      analysis_result(NULL)
      analysis_error(NULL)
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

      submission <- tryCatch(
        structure_submit_background_job(
          upload_context = saved_input,
          controls = list(
            sample = input$structure_sample,
            coordinate = coordinate_value,
            window_size = structure_controls_defaults()$window_size
          )
        ),
        error = function(err) {
          analysis_error(conditionMessage(err))
          NULL
        }
      )

      if (is.null(submission)) {
        publish_controls()
        return(invisible(NULL))
      }

      loaded_upload_context(saved_input)
      active_job_id(submission$job_id)
      structure_job_state(submission$state)
      updateTextInput(session, "job_lookup", value = submission$job_id)
      publish_controls()
      publish_export_controls()
    }, ignoreInit = TRUE)

    observeEvent(input$structure_sample, {
      publish_controls()
    }, ignoreInit = TRUE)

    observeEvent(current_upload(), {
      loaded_upload_context(NULL)
      publish_controls()
      status_message("Saved BED context changed. Click Run Structure to load it.")
    }, ignoreInit = TRUE)

    observeEvent(input$controls_ready, {
      publish_controls()
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
      current_job_id <- active_job_id()
      current_state <- structure_job_state()

      if (!nzchar(current_job_id) || is.null(current_state) || !current_state$status %in% c("queued", "running")) {
        return(invisible(NULL))
      }

      invalidateLater(2000, session)
      refresh_job_from_disk(current_job_id)
      publish_controls()
      publish_export_controls()
    })

    observe({
      invalidateLater(2000, session)
      queue_status(structure_job_capacity_snapshot())
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

    observe({
      set_native_button_state(
        session = session,
        id = session$ns("run_structure"),
        disabled = isTRUE(is_submitting()) ||
          !structure_has_saved_upload(current_upload()) ||
          !has_structure_coordinate() ||
          is_analysis_locked()
      )
    })

    output$status_message <- renderText({
      status_message()
    })

    output$run_hint <- renderUI({
      saved_input <- current_upload()
      has_upload <- structure_has_saved_upload(saved_input)
      job_lookup <- trimws(as.character(input$job_lookup %||% ""))

      message <- NULL
      if (!has_upload && !nzchar(job_lookup)) {
        message <- "Enter a Structure Job ID, upload a BED file in Upload / Run, or use Demo."
      } else if (has_upload && !has_structure_coordinate()) {
        saved_files <- upload_run_context_saved_files(saved_input)
        sample_index <- suppressWarnings(as.integer(input$structure_sample %||% 1L))
        if (!is.finite(sample_index)) {
          sample_index <- 1L
        }
        sample_index <- min(max(sample_index, 1L), length(saved_files))
        placeholder <- if (length(saved_files)) structure_coordinate_placeholder(saved_files[[sample_index]]) else "chrom:start-end"
        message <- sprintf("Enter Coordinate before running Structure, for example %s.", placeholder)
      }

      if (is.null(message)) {
        return(NULL)
      }

      tags$p(class = "structure-run-hint", message)
    })

    output$parameter_snapshot <- renderUI({
      saved_input <- current_upload()
      file_count <- if (is.null(saved_input)) 0 else upload_run_context_file_count(saved_input)
      species_label <- if (is.null(saved_input) || is.null(saved_input$species) || !nzchar(saved_input$species)) {
        "Not saved yet"
      } else {
        upload_run_display_species(saved_input$species)
      }
      interval_label <- if (is.null(loaded_upload_context())) "Load on Run" else "Coordinate input"
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
        return(structure_job_notice_card(job_state = state, error_message = error_message))
      }

      tagList(
        structure_job_notice_card(job_state = state, error_message = error_message),
        div(
          class = "structure-result-summary",
          tags$strong(sprintf(
            "Job ID: %s",
            state$jobId %||% result$jobId %||% "Structure job"
          )),
          tags$span(sprintf(
            "%s | %s | %s",
            result$method %||% "RNAfold",
            state$status %||% "completed",
            result$locus %||% ""
          ))
        ),
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
      publish_export_controls()
      status_message("Save a BED file in Upload / Run, then click Run Structure.")
    }, once = TRUE)

    outputOptions(output, "run_hint", suspendWhenHidden = FALSE)
    outputOptions(output, "analysis_output", suspendWhenHidden = FALSE)

    invisible(NULL)
  })
}
