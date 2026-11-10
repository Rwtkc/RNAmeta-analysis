mod_meta_plot_server <- function(id, upload_context = NULL, analysis_lock = NULL) {
  moduleServer(id, function(input, output, session) {
    defaults <- meta_plot_controls_config(id)$defaults
    analysis_result <- reactiveVal(NULL)
    analysis_error <- reactiveVal(NULL)
    module_lock_id <- session$ns("analysis_lock")
    session_cache_dirs <- character()
    session_temp_files <- character()

    session$onSessionEnded(function() {
      paths <- unique(session_cache_dirs)
      if (length(paths) && exists("rnameta_unlink_cache_paths", mode = "function")) {
        rnameta_unlink_cache_paths(paths, recursive = TRUE)
      } else if (length(paths)) {
        unlink(paths, recursive = TRUE, force = TRUE)
      }

      temp_paths <- unique(session_temp_files)
      if (length(temp_paths) && exists("rnameta_unlink_temp_paths", mode = "function")) {
        rnameta_unlink_temp_paths(temp_paths)
      } else if (length(temp_paths)) {
        unlink(temp_paths, force = TRUE)
      }

      invisible(NULL)
    })

    apply_chart_tuning <- function(payload) {
      if (is.null(payload)) {
        return(payload)
      }

      payload$componentOverrides <- list(
        promoter = list(scale = 1.0, fill = "#111111", alpha = 1),
        tail = list(scale = 1.0, fill = "#111111", alpha = 1)
      )

      if (!is.null(payload$yDomain) && length(payload$yDomain) >= 2) {
        payload$yDomain[[1]] <- -0.6
      }

      payload
    }

    current_controls <- reactive({
      meta_plot_control_values(input, defaults)
    })

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

    publish_controls <- function() {
      send_control_payload(
        session = session,
        id = session$ns("controls_host"),
        control = "meta-plot-controls",
        config = meta_plot_controls_config(id)
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
        control = "meta-plot-chart",
        config = apply_chart_tuning(result$chart_payload)
      )
    }

    publish_export_controls <- function() {
      send_control_payload(
        session = session,
        id = session$ns("export_host"),
        control = "meta-plot-export",
        config = meta_plot_export_config(id)
      )
    }

    session$onFlushed(function() {
      isolate(publish_controls())
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

    observeEvent(input$run_plot, {
      req(has_saved_upload())
      req(!is_analysis_locked())

      analysis_error(NULL)
      analysis_id <- meta_plot_new_analysis_id()
      if (!is.null(analysis_lock)) {
        analysis_lock(module_lock_id)
      }
      progress <- meta_plot_progress_start(session)
      on.exit({
        meta_plot_progress_close(progress)
        if (!is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
      }, add = TRUE)

      result <- tryCatch(
        run_meta_plot_analysis(
          upload_context = if (is.null(upload_context)) NULL else upload_context(),
          controls = current_controls(),
          progress = progress
        ),
        error = function(err) {
          analysis_error(conditionMessage(err))
          NULL
        }
      )

      if (!is.null(result) && !is.null(result$chart_payload)) {
        result$chart_payload <- apply_chart_tuning(result$chart_payload)
        cache_write <- write_meta_plot_chart_cache(result$chart_payload, analysis_id = analysis_id)
        result$analysis_id <- cache_write$analysis_id
        result$chart_cache_path <- cache_write$path
        session_cache_dirs <<- unique(c(session_cache_dirs, dirname(cache_write$path)))
      }

      if (!is.null(result) && length(result$prepared_bed_files %||% character())) {
        session_temp_files <<- unique(c(session_temp_files, result$prepared_bed_files))
      }

      analysis_result(result)
    })

    observeEvent(input$export_plot, {
      result <- analysis_result()
      req(!is.null(result), !is.null(result$chart_payload))

      settings <- meta_plot_export_settings(input)

      session$sendCustomMessage(
        "meta-plot-export",
        list(
          hostId = session$ns("chart_host"),
          format = settings$format,
          filename = meta_plot_export_filename(settings$format),
          width = settings$width,
          height = settings$height,
          dpi = settings$dpi
        )
      )
    })

    observeEvent(input$data_export, {
      result <- analysis_result()
      req(!is.null(result), !is.null(result$chart_payload))

      settings <- rnameta_data_export_settings(input$data_export_format)

      session$sendCustomMessage(
        "rnameta-data-export",
        list(
          filename = meta_plot_data_export_filename(settings$extension),
          mimeType = rnameta_data_export_mime_type(settings$format),
          content = rnameta_data_export_content(meta_plot_data_export_table(result), delimiter = settings$delimiter)
        )
      )
    })

    observe({
      set_native_button_state(
        session = session,
        id = session$ns("run_plot"),
        disabled = !has_saved_upload() || is_analysis_locked()
      )
    })

    output$run_hint <- renderUI({
      if (has_saved_upload()) {
        return(NULL)
      }

      tags$p(
        class = "meta-plot-sidebar__hint",
        "Save a BED file in Upload / Run first to enable Meta Plot analysis."
      )
    })

    observe({
      analysis_result()
      publish_chart()
    })

    if (!is.null(upload_context)) {
      observeEvent(upload_context(), {
        analysis_result(NULL)
        analysis_error(NULL)
      }, ignoreInit = TRUE)
    }

    output$parameter_snapshot <- renderUI({
      controls <- current_controls()
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      file_count <- if (is.null(saved_input)) 0 else upload_run_context_file_count(saved_input)
      transcript_label <- switch(
        controls$pltTxType,
        tx = "exon",
        mrna = "mRNA",
        ncrna = "lncRNA",
        controls$pltTxType
      )

      div(
        class = "meta-plot-summary-grid",
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Transcript Mode"),
          span(class = "meta-plot-summary-item__value", transcript_label)
        ),
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Promoter / Tail"),
          span(class = "meta-plot-summary-item__value", controls$headOrtail)
        ),
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Confidence Interval"),
          span(class = "meta-plot-summary-item__value", controls$enableCI)
        ),
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Saved BED Files"),
          span(
            class = "meta-plot-summary-item__value",
            if (is.null(saved_input) || file_count == 0) "Not saved yet" else sprintf("%d file(s)", file_count)
          )
        )
      )
    })

    output$analysis_notice <- renderUI({
      if (!is.null(analysis_error())) {
        return(
          div(
            class = "meta-plot-result-card meta-plot-result-card--error",
            tags$p(class = "meta-plot-result-card__copy", analysis_error())
          )
        )
      }

      NULL
    })

    output$export_panel <- renderUI({
      req(!is.null(analysis_result()))

      react_control_host(
        id = session$ns("export_host"),
        control_type = "meta-plot-export",
        class = "meta-plot-export-root",
        ready_input_id = session$ns("export_ready")
      )
    })

    output$analysis_summary <- renderUI({
      result <- analysis_result()
      req(!is.null(result))

      div(
        class = "meta-plot-summary-grid meta-plot-summary-grid--result",
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Species"),
          span(class = "meta-plot-summary-item__value", result$species)
        ),
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Samples"),
          span(class = "meta-plot-summary-item__value", result$sample_count)
        ),
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Intervals"),
          span(class = "meta-plot-summary-item__value", result$interval_count)
        ),
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Sampled Points"),
          span(class = "meta-plot-summary-item__value", result$sampled_point_count)
        ),
        div(
          class = "meta-plot-summary-item",
          span(class = "meta-plot-summary-item__label", "Transcript Hits"),
          span(class = "meta-plot-summary-item__value", result$overlap_count)
        )
      )
    })

    outputOptions(output, "run_hint", suspendWhenHidden = FALSE)

  })
}
