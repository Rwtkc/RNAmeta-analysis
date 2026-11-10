mod_transcription_server <- function(id, upload_context = NULL, analysis_lock = NULL) {
  moduleServer(id, function(input, output, session) {
    analysis_result <- reactiveVal(NULL)
    analysis_error <- reactiveVal(NULL)
    module_lock_id <- session$ns("analysis_lock")

    has_saved_upload <- reactive({
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      saved_files <- upload_run_context_saved_files(saved_input)
      !is.null(saved_input) &&
        length(saved_files) > 0 &&
        all(vapply(saved_files, function(file) !is.null(file$file_path) && nzchar(file$file_path) && file.exists(file$file_path), logical(1)))
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

    publish_chart <- function() {
      result <- analysis_result()

      if (is.null(result) || is.null(result$chart_payload)) {
        clear_control_payload(session = session, id = session$ns("chart_host"))
        return(invisible(NULL))
      }

      send_control_payload(
        session = session,
        id = session$ns("chart_host"),
        control = "transcription-chart",
        config = result$chart_payload
      )
    }

    publish_export_controls <- function() {
      send_control_payload(
        session = session,
        id = session$ns("export_host"),
        control = "transcription-export",
        config = transcription_export_config(id)
      )
    }

    observeEvent(input$run_transcription, {
      req(has_saved_upload())
      req(!is_analysis_locked())

      analysis_error(NULL)
      if (!is.null(analysis_lock)) {
        analysis_lock(module_lock_id)
      }
      progress <- transcription_progress_start(session)
      on.exit({
        transcription_progress_close(progress)
        if (!is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
      }, add = TRUE)

      result <- tryCatch(
        run_transcription_analysis(
          upload_context = if (is.null(upload_context)) NULL else upload_context(),
          progress = progress
        ),
        error = function(err) {
          analysis_error(conditionMessage(err))
          NULL
        }
      )

      analysis_result(result)
    })

    observeEvent(input$export_plot, {
      result <- analysis_result()
      req(!is.null(result), !is.null(result$chart_payload))

      settings <- transcription_export_settings(input)

      session$sendCustomMessage(
        "transcription-export",
        list(
          hostId = session$ns("chart_host"),
          format = settings$format,
          filename = transcription_export_filename(settings$format),
          width = settings$width,
          height = settings$height,
          dpi = settings$dpi
        )
      )
    })

    observeEvent(input$data_export, {
      result <- analysis_result()
      req(!is.null(result))

      settings <- rnameta_data_export_settings(input$data_export_format)

      session$sendCustomMessage(
        "rnameta-data-export",
        list(
          filename = transcription_data_export_filename(settings$extension),
          mimeType = rnameta_data_export_mime_type(settings$format),
          content = rnameta_data_export_content(transcription_data_export_table(result), delimiter = settings$delimiter)
        )
      )
    })

    observeEvent(input$chart_ready, {
      publish_chart()
    }, ignoreInit = TRUE)

    observeEvent(input$export_ready, {
      publish_export_controls()
    }, ignoreInit = TRUE)

    if (!is.null(upload_context)) {
      observeEvent(upload_context(), {
        analysis_result(NULL)
        analysis_error(NULL)
      }, ignoreInit = TRUE)
    }

    observe({
      analysis_result()
      publish_chart()
    })

    observe({
      set_native_button_state(
        session = session,
        id = session$ns("run_transcription"),
        disabled = !has_saved_upload() || is_analysis_locked()
      )
    })

    output$run_hint <- renderUI({
      if (has_saved_upload()) {
        return(NULL)
      }

      tags$p(
        class = "site-profile-sidebar__hint",
        "Save a BED file in Upload / Run first to enable Transcription analysis."
      )
    })

    output$export_panel <- renderUI({
      result <- analysis_result()
      req(!is.null(result), isTRUE(result$has_signal))

      react_control_host(
        id = session$ns("export_host"),
        control_type = "transcription-export",
        class = "site-profile-export-root",
        ready_input_id = session$ns("export_ready")
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
      if (!is.null(analysis_error())) {
        return(
          div(
            class = "site-profile-result-card site-profile-result-card--error",
            tags$p(class = "site-profile-result-card__copy", analysis_error())
          )
        )
      }

      result <- analysis_result()
      if (!is.null(result) && !isTRUE(result$has_signal) && !is.null(result$empty_state_message)) {
        return(
          div(
            class = "site-profile-result-card site-profile-result-card--info",
            tags$p(class = "site-profile-result-card__copy", result$empty_state_message)
          )
        )
      }

      if (!is.null(result) && isTRUE(result$has_signal) && !isTRUE(result$has_heatmap_signal) && !is.null(result$heatmap_notice_message)) {
        return(
          div(
            class = "site-profile-result-card site-profile-result-card--info",
            tags$p(class = "site-profile-result-card__copy", result$heatmap_notice_message)
          )
        )
      }

      NULL
    })

    output$analysis_summary <- renderUI({
      result <- analysis_result()
      req(!is.null(result))

      div(
        class = "site-profile-summary-grid site-profile-summary-grid--result",
        div(
          class = "site-profile-summary-item",
          span(class = "site-profile-summary-item__label", "Species"),
          span(class = "site-profile-summary-item__value", result$species)
        ),
        div(
          class = "site-profile-summary-item",
          span(class = "site-profile-summary-item__label", "Samples"),
          span(class = "site-profile-summary-item__value", result$sample_count)
        ),
        div(
          class = "site-profile-summary-item",
          span(class = "site-profile-summary-item__label", "Intervals"),
          span(class = "site-profile-summary-item__value", result$interval_count)
        ),
        div(
          class = "site-profile-summary-item",
          span(class = "site-profile-summary-item__label", "Heatmap Sample"),
          span(
            class = "site-profile-summary-item__value",
            if (!isTRUE(result$has_signal)) {
              "No boundary hits"
            } else if (!isTRUE(result$has_heatmap_signal)) {
              "No heatmap signal"
            } else {
              result$heatmap_sample
            }
          )
        )
      )
    })

    outputOptions(output, "run_hint", suspendWhenHidden = FALSE)
  })
}
