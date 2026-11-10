mod_motif_server <- function(id, upload_context = NULL, analysis_lock = NULL) {
  moduleServer(id, function(input, output, session) {
    analysis_result <- reactiveVal(NULL)
    analysis_error <- reactiveVal(NULL)
    job_state <- reactiveVal(NULL)
    queue_status <- reactiveVal(motif_job_capacity_snapshot())
    active_job_id <- reactiveVal("")
    is_submitting <- reactiveVal(FALSE)
    module_lock_id <- session$ns("analysis_lock")

    has_saved_upload <- reactive({
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      saved_files <- upload_run_context_saved_files(saved_input)
      !is.null(saved_input) &&
        length(saved_files) > 0 &&
        all(vapply(saved_files, function(file) !is.null(file$file_path) && nzchar(file$file_path) && file.exists(file$file_path), logical(1)))
    })

    sample_choices <- reactive({
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      motif_build_sample_choices(saved_input)
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

    current_controls <- reactive({
      list(
        sample_name = input$sample_name,
        engine = input$engine,
        analysis_mode = input$analysis_mode,
        sites_subset = input$sites_subset,
        motif_width = input$motif_width,
        sequence_window = input$sequence_window,
        max_motifs = input$max_motifs,
        job_lookup = input$job_lookup
      )
    })

    refresh_job_from_disk <- function(job_id) {
      job_bundle <- motif_job_query(job_id)
      if (is.null(job_bundle) || is.null(job_bundle$state)) {
        analysis_error(sprintf("Motif job id '%s' was not found.", job_id))
        return(invisible(NULL))
      }

      state <- job_bundle$state
      active_job_id(state$jobId %||% job_id)
      job_state(state)

      if (identical(state$status, "completed")) {
        if (is.null(job_bundle$result)) {
          analysis_result(NULL)
          analysis_error(sprintf("Motif job '%s' is marked completed but no result payload was found.", job_id))
          return(invisible(NULL))
        }

        analysis_error(NULL)
        analysis_result(job_bundle$result)
        return(invisible(job_bundle$result))
      }

      analysis_result(NULL)
      if (identical(state$status, "failed")) {
        analysis_error(state$errorMessage %||% state$message %||% "Motif analysis failed.")
      } else {
        analysis_error(NULL)
      }

      invisible(state)
    }

    publish_controls <- function() {
      saved_input <- if (is.null(upload_context)) NULL else upload_context()

      send_control_payload(
        session = session,
        id = session$ns("controls_host"),
        control = "motif-controls",
        config = motif_controls_config(
          id = id,
          upload_context = saved_input,
          current_values = current_controls(),
          active_job_state = job_state(),
          queue_status = queue_status()
        )
      )
    }

    publish_chart <- function() {
      result <- analysis_result()
      payload <- if (!is.null(result) && !is.null(result$chart_payload)) {
        result$chart_payload
      } else {
        motif_placeholder_chart_payload(job_state())
      }

      if (is.null(payload)) {
        clear_control_payload(session = session, id = session$ns("chart_host"))
        return(invisible(NULL))
      }

      send_control_payload(
        session = session,
        id = session$ns("chart_host"),
        control = "motif-result",
        config = payload
      )
    }

    publish_export_controls <- function() {
      result <- analysis_result()
      if (is.null(result) || is.null(result$motifs) || !length(result$motifs)) {
        clear_control_payload(session = session, id = session$ns("export_host"))
        return(invisible(NULL))
      }

      send_control_payload(
        session = session,
        id = session$ns("export_host"),
        control = "motif-export",
        config = motif_export_config(id)
      )
    }

    session$onFlushed(function() {
      isolate({
        publish_controls()
        publish_chart()
      })
    }, once = TRUE)

    observeEvent(input$controls_ready, {
      publish_controls()
    }, ignoreInit = TRUE)

    observeEvent(input$chart_ready, {
      publish_chart()
    }, ignoreInit = TRUE)

    observeEvent(input$export_ready, {
      publish_export_controls()
    }, ignoreInit = TRUE)

    observeEvent(input$engine, {
      updateNumericInput(
        session = session,
        inputId = "motif_width",
        min = motif_min_width(input$engine),
        value = motif_effective_width(input$engine, input$motif_width)
      )
    }, ignoreInit = TRUE)

    observeEvent(input$run_motif, {
      req(has_saved_upload())
      req(!is_analysis_locked())
      req(!is_submitting())

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

      submission <- tryCatch({
        controls <- motif_validate_analysis_controls(current_controls())
        motif_submit_background_job(
          upload_context = if (is.null(upload_context)) NULL else upload_context(),
          controls = controls
        )
      },
        error = function(err) {
          analysis_error(conditionMessage(err))
          NULL
        }
      )

      if (is.null(submission)) {
        return(invisible(NULL))
      }

      active_job_id(submission$job_id)
      job_state(submission$state)
      updateTextInput(session, "job_lookup", value = submission$job_id)
      publish_controls()
      publish_chart()
    }, ignoreInit = TRUE)

    observeEvent(input$query_job, {
      requested_id <- motif_job_normalize_id(input$job_lookup)
      if (!nzchar(requested_id)) {
        analysis_error("Enter a motif job id before querying.")
        return(invisible(NULL))
      }

      refresh_job_from_disk(requested_id)
      publish_controls()
      publish_chart()
    }, ignoreInit = TRUE)

    observe({
      current_job_id <- active_job_id()
      current_state <- job_state()

      if (!nzchar(current_job_id) || is.null(current_state) || !current_state$status %in% c("queued", "running")) {
        return(invisible(NULL))
      }

      invalidateLater(2000, session)
      refresh_job_from_disk(current_job_id)
    })

    observe({
      invalidateLater(2000, session)
      queue_status(motif_job_capacity_snapshot())
    })

    observeEvent(input$export_plot, {
      result <- analysis_result()
      req(!is.null(result))
      req(!is.null(result$motifs))
      req(length(result$motifs) > 0)

      settings <- motif_export_settings(input)
      session$sendCustomMessage(
        "motif-export",
        list(
          hostId = session$ns("chart_host"),
          triggerId = session$ns("export_plot"),
          format = settings$format,
          width = settings$width,
          height = settings$height,
          dpi = settings$dpi,
          filename = motif_export_filename(result$engine, result$sample_name, settings$format),
          analysis = result$analysis_metadata,
          motifs = result$motifs
        )
      )
    }, ignoreInit = TRUE)

    if (!is.null(upload_context)) {
      observeEvent(upload_context(), {
        analysis_result(NULL)
        analysis_error(NULL)
        job_state(NULL)
        active_job_id("")
        choices <- sample_choices()
        updateSelectInput(session, "sample_name", choices = choices, selected = if (length(choices)) unname(choices[[1]]) else "")
        updateTextInput(session, "job_lookup", value = "")
        publish_controls()
        publish_chart()
      }, ignoreInit = TRUE)
    }

    observe({
      analysis_result()
      job_state()
      publish_chart()
      publish_export_controls()
    })

    observe({
      job_state()
      publish_controls()
    })

    observe({
      queue_status()
      publish_controls()
    })

    observeEvent(current_controls(), {
      publish_controls()
    }, ignoreInit = TRUE)

    observe({
      set_native_button_state(
        session = session,
        id = session$ns("run_motif"),
        disabled = !has_saved_upload() || is_analysis_locked() || is_submitting()
      )
    })

    output$run_hint <- renderUI({
      if (has_saved_upload()) {
        return(NULL)
      }

      tags$p(
        class = "site-profile-sidebar__hint",
        "Save a BED file in Upload / Run first to enable motif job submission."
      )
    })

    output$parameter_snapshot <- renderUI({
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      file_count <- if (is.null(saved_input)) 0 else upload_run_context_file_count(saved_input)

      div(
        class = "site-profile-summary-grid",
        div(
          class = "site-profile-summary-item",
          span(class = "site-profile-summary-item__label", "Saved BED Files"),
          span(class = "site-profile-summary-item__value", if (is.null(saved_input) || file_count == 0) "Not saved yet" else sprintf("%d file(s)", file_count))
        ),
        div(
          class = "site-profile-summary-item",
          span(class = "site-profile-summary-item__label", "Species"),
          span(class = "site-profile-summary-item__value", if (is.null(saved_input)) "Not saved yet" else saved_input$species)
        )
      )
    })

    output$analysis_notice <- renderUI({
      motif_job_notice_card(job_state = job_state(), error_message = analysis_error())
    })

    output$export_panel <- renderUI({
      result <- analysis_result()
      req(!is.null(result))
      req(!is.null(result$motifs))
      req(length(result$motifs) > 0)

      react_control_host(
        id = session$ns("export_host"),
        control_type = "motif-export",
        class = "motif-export-root",
        ready_input_id = session$ns("export_ready")
      )
    })

    outputOptions(output, "run_hint", suspendWhenHidden = FALSE)
  })
}
