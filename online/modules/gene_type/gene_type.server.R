mod_gene_type_server <- function(id, upload_context = NULL, analysis_lock = NULL) {
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
        control = "gene-type-chart",
        config = result$chart_payload
      )
    }

    publish_export_controls <- function() {
      send_control_payload(
        session = session,
        id = session$ns("export_host"),
        control = "gene-type-export",
        config = gene_type_export_config(id)
      )
    }

    observeEvent(input$run_gene_type, {
      req(has_saved_upload())
      req(!is_analysis_locked())

      analysis_error(NULL)
      if (!is.null(analysis_lock)) {
        analysis_lock(module_lock_id)
      }
      progress <- gene_type_progress_start(session)
      on.exit({
        gene_type_progress_close(progress)
        if (!is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
      }, add = TRUE)

      result <- tryCatch(
        run_gene_type_analysis(
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

      settings <- gene_type_export_settings(input)

      session$sendCustomMessage(
        "gene-type-export",
        list(
          hostId = session$ns("chart_host"),
          format = settings$format,
          filename = gene_type_export_filename(settings$format),
          width = settings$width,
          height = settings$height,
          dpi = settings$dpi
        )
      )
    })

    observeEvent(input$data_export, {
      result <- analysis_result()
      req(!is.null(result), !is.null(result$gene_type_counts))

      settings <- rnameta_data_export_settings(input$data_export_format)

      session$sendCustomMessage(
        "rnameta-data-export",
        list(
          filename = gene_type_data_export_filename(settings$extension),
          mimeType = rnameta_data_export_mime_type(settings$format),
          content = rnameta_data_export_content(gene_type_data_export_table(result), delimiter = settings$delimiter)
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
        id = session$ns("run_gene_type"),
        disabled = !has_saved_upload() || is_analysis_locked()
      )
    })

    output$run_hint <- renderUI({
      if (has_saved_upload()) {
        return(NULL)
      }

      tags$p(
        class = "gene-type-sidebar__hint",
        "Save a BED file in Upload / Run first to enable Gene Type analysis."
      )
    })

    output$export_panel <- renderUI({
      req(!is.null(analysis_result()))

      react_control_host(
        id = session$ns("export_host"),
        control_type = "gene-type-export",
        class = "gene-type-export-root",
        ready_input_id = session$ns("export_ready")
      )
    })

    output$parameter_snapshot <- renderUI({
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      file_count <- if (is.null(saved_input)) 0 else upload_run_context_file_count(saved_input)

      div(
        class = "gene-type-summary-grid",
        div(
          class = "gene-type-summary-item",
          span(class = "gene-type-summary-item__label", "Saved BED Files"),
          span(class = "gene-type-summary-item__value", if (is.null(saved_input) || file_count == 0) "Not saved yet" else sprintf("%d file(s)", file_count))
        ),
        div(
          class = "gene-type-summary-item",
          span(class = "gene-type-summary-item__label", "Species"),
          span(class = "gene-type-summary-item__value", if (is.null(saved_input)) "Not saved yet" else saved_input$species)
        )
      )
    })

    output$analysis_notice <- renderUI({
      if (!is.null(analysis_error())) {
        return(
          div(
            class = "gene-type-result-card gene-type-result-card--error",
            tags$p(class = "gene-type-result-card__copy", analysis_error())
          )
        )
      }

      NULL
    })

    output$analysis_summary <- renderUI({
      result <- analysis_result()
      req(!is.null(result))

      div(
        class = "gene-type-summary-grid gene-type-summary-grid--result",
        div(
          class = "gene-type-summary-item",
          span(class = "gene-type-summary-item__label", "Species"),
          span(class = "gene-type-summary-item__value", result$species)
        ),
        div(
          class = "gene-type-summary-item",
          span(class = "gene-type-summary-item__label", "Samples"),
          span(class = "gene-type-summary-item__value", result$sample_count)
        ),
        div(
          class = "gene-type-summary-item",
          span(class = "gene-type-summary-item__label", "Intervals"),
          span(class = "gene-type-summary-item__value", result$interval_count)
        ),
        div(
          class = "gene-type-summary-item",
          span(class = "gene-type-summary-item__label", "Biotype Hits"),
          span(class = "gene-type-summary-item__value", result$overlap_count)
        )
      )
    })

    outputOptions(output, "run_hint", suspendWhenHidden = FALSE)
  })
}
