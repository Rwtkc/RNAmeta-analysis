mod_peak_distribution_server <- function(id, upload_context = NULL, analysis_lock = NULL) {
  moduleServer(id, function(input, output, session) {
    analysis_result <- reactiveVal(NULL)
    analysis_error <- reactiveVal(NULL)
    detail_output_dir <- reactiveVal(NULL)
    detail_output_dir_state <- new.env(parent = emptyenv())
    detail_output_dir_state$path <- NULL
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

    current_controls <- reactive({
      selected_features <- input$selected_features
      if (is.null(selected_features) || length(selected_features) == 0) {
        selected_features <- peak_distribution_feature_order()
      }

      list(
        selected_features = selected_features,
        detail_output_dir = detail_output_dir()
      )
    })

    build_detail_output_dir <- function() {
      old_dir <- detail_output_dir_state$path
      if (!is.null(old_dir) && dir.exists(old_dir)) {
        unlink(old_dir, recursive = TRUE, force = TRUE)
      }

      new_dir <- tempfile(pattern = "peak_distribution_", tmpdir = tempdir())
      dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)
      detail_output_dir_state$path <- new_dir
      detail_output_dir(new_dir)
      new_dir
    }

    publish_controls <- function() {
      send_control_payload(
        session = session,
        id = session$ns("controls_host"),
        control = "peak-distribution-controls",
        config = peak_distribution_controls_config(id, current_controls())
      )
    }

    publish_chart <- function() {
      result <- analysis_result()

      if (is.null(result) || is.null(result$chart_payload)) {
        clear_control_payload(session = session, id = session$ns("chart_host"))
        return(invisible(NULL))
      }

      send_control_payload(
        session = session,
        id = session$ns("chart_host"),
        control = "peak-distribution-chart",
        config = result$chart_payload
      )
    }

    publish_export_controls <- function() {
      send_control_payload(
        session = session,
        id = session$ns("export_host"),
        control = "peak-distribution-export",
        config = peak_distribution_export_config(id)
      )
    }

    session$onSessionEnded(function() {
      output_dir <- detail_output_dir_state$path
      if (!is.null(output_dir) && dir.exists(output_dir)) {
        unlink(output_dir, recursive = TRUE, force = TRUE)
      }
    })

    session$onFlushed(function() {
      isolate(publish_controls())
    }, once = TRUE)

    observeEvent(input$run_peak_distribution, {
      req(has_saved_upload())
      req(!is_analysis_locked())

      analysis_error(NULL)
      build_detail_output_dir()
      if (!is.null(analysis_lock)) {
        analysis_lock(module_lock_id)
      }
      progress <- peak_distribution_progress_start(session)
      on.exit({
        peak_distribution_progress_close(progress)
        if (!is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
      }, add = TRUE)

      result <- tryCatch(
        run_peak_distribution_analysis(
          upload_context = if (is.null(upload_context)) NULL else upload_context(),
          controls = current_controls(),
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

      settings <- peak_distribution_export_settings(input)

      session$sendCustomMessage(
        "peak-distribution-export",
        list(
          hostId = session$ns("chart_host"),
          format = settings$format,
          filename = peak_distribution_export_filename(settings$format),
          width = settings$width,
          height = settings$height,
          dpi = settings$dpi
        )
      )
    })

    observeEvent(input$data_export, {
      result <- analysis_result()
      req(!is.null(result), !is.null(result$feature_counts))

      settings <- peak_distribution_data_export_settings(input)

      session$sendCustomMessage(
        "rnameta-data-export",
        list(
          filename = peak_distribution_data_export_filename(settings$extension),
          mimeType = peak_distribution_data_export_mime_type(settings$format),
          content = peak_distribution_data_export_content(result, delimiter = settings$delimiter)
        )
      )
    })

    observeEvent(input$controls_ready, {
      publish_controls()
    }, ignoreInit = TRUE)

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
      current_controls()
      publish_controls()
    })

    observe({
      analysis_result()
      publish_chart()
    })

    observe({
      set_native_button_state(
        session = session,
        id = session$ns("run_peak_distribution"),
        disabled = !has_saved_upload() || is_analysis_locked()
      )
    })

    output$run_hint <- renderUI({
      if (has_saved_upload()) {
        return(NULL)
      }

      tags$p(
        class = "peak-distribution-sidebar__hint",
        "Save a BED file in Upload / Run first to enable Peak Distribution analysis."
      )
    })

    output$export_panel <- renderUI({
      req(!is.null(analysis_result()))

      react_control_host(
        id = session$ns("export_host"),
        control_type = "peak-distribution-export",
        class = "peak-distribution-export-root",
        ready_input_id = session$ns("export_ready")
      )
    })

    output$parameter_snapshot <- renderUI({
      controls <- current_controls()
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      file_count <- if (is.null(saved_input)) 0 else upload_run_context_file_count(saved_input)
      feature_label <- sprintf("%d selected", length(controls$selected_features))

      div(
        class = "peak-distribution-summary-grid",
        div(
          class = "peak-distribution-summary-item",
          span(class = "peak-distribution-summary-item__label", "Visible Features"),
          span(class = "peak-distribution-summary-item__value", feature_label)
        ),
        div(
          class = "peak-distribution-summary-item",
          span(class = "peak-distribution-summary-item__label", "Saved BED Files"),
          span(class = "peak-distribution-summary-item__value", if (is.null(saved_input) || file_count == 0) "Not saved yet" else sprintf("%d file(s)", file_count))
        )
      )
    })

    output$analysis_notice <- renderUI({
      if (!is.null(analysis_error())) {
        return(
          div(
            class = "peak-distribution-result-card peak-distribution-result-card--error",
            tags$p(class = "peak-distribution-result-card__copy", analysis_error())
          )
        )
      }

      NULL
    })

    output$analysis_summary <- renderUI({
      result <- analysis_result()
      req(!is.null(result))

      div(
        class = "peak-distribution-summary-grid peak-distribution-summary-grid--result",
        div(
          class = "peak-distribution-summary-item",
          span(class = "peak-distribution-summary-item__label", "Species"),
          span(class = "peak-distribution-summary-item__value", result$species)
        ),
        div(
          class = "peak-distribution-summary-item",
          span(class = "peak-distribution-summary-item__label", "Samples"),
          span(class = "peak-distribution-summary-item__value", result$sample_count)
        ),
        div(
          class = "peak-distribution-summary-item",
          span(class = "peak-distribution-summary-item__label", "Intervals"),
          span(class = "peak-distribution-summary-item__value", result$interval_count)
        ),
        div(
          class = "peak-distribution-summary-item",
          span(class = "peak-distribution-summary-item__label", "Feature Rows"),
          span(class = "peak-distribution-summary-item__value", nrow(result$feature_counts))
        )
      )
    })

    output$chart_empty_state <- renderUI({
      NULL
    })

    outputOptions(output, "run_hint", suspendWhenHidden = FALSE)
  })
}
