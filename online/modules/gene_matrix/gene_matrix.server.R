mod_gene_matrix_server <- function(id, upload_context = NULL, analysis_lock = NULL) {
  moduleServer(id, function(input, output, session) {
    analysis_result <- reactiveVal(NULL)
    analysis_error <- reactiveVal(NULL)
    module_lock_id <- session$ns("analysis_lock")

    saved_files <- reactive({
      context <- if (is.null(upload_context)) NULL else upload_context()
      gene_matrix_saved_files(context)
    })

    has_saved_upload <- reactive({
      files <- saved_files()
      gene_matrix_sample_count_ok(length(files)) &&
        all(vapply(files, function(file) !is.null(file$file_path) && nzchar(file$file_path) && file.exists(file$file_path), logical(1)))
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
        control = "gene-matrix-chart",
        config = result$chart_payload
      )
    }

    publish_export_controls <- function() {
      send_control_payload(
        session = session,
        id = session$ns("export_host"),
        control = "gene-matrix-export",
        config = gene_matrix_export_config(id)
      )
    }

    observeEvent(input$run_gene_matrix, {
      req(has_saved_upload())
      req(!is_analysis_locked())

      analysis_error(NULL)
      if (!is.null(analysis_lock)) {
        analysis_lock(module_lock_id)
      }
      progress <- gene_matrix_progress_start(session)
      on.exit({
        gene_matrix_progress_close(progress)
        if (!is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
      }, add = TRUE)

      result <- tryCatch(
        run_gene_matrix_analysis(
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

      settings <- gene_matrix_export_settings(input)
      session$sendCustomMessage(
        "gene-matrix-export",
        list(
          hostId = session$ns("chart_host"),
          format = settings$format,
          filename = gene_matrix_export_filename(settings$format),
          width = settings$width,
          height = settings$height,
          dpi = settings$dpi
        )
      )
    })

    observeEvent(input$data_export, {
      result <- analysis_result()
      req(!is.null(result), !is.null(result$gene_sets))

      settings <- rnameta_data_export_settings(input$data_export_format)

      session$sendCustomMessage(
        "rnameta-data-export",
        list(
          filename = gene_matrix_data_export_filename(settings$extension),
          mimeType = rnameta_data_export_mime_type(settings$format),
          content = rnameta_data_export_content(gene_matrix_data_export_table(result), delimiter = settings$delimiter)
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
        id = session$ns("run_gene_matrix"),
        disabled = !has_saved_upload() || is_analysis_locked()
      )
    })

    output$run_hint <- renderUI({
      file_count <- length(saved_files())

      if (has_saved_upload()) {
        return(NULL)
      }

      tags$p(
        class = "gene-matrix-sidebar__hint",
        if (!file_count) {
          "Save 2 to 5 BED files in Upload / Run first to enable Gene Matrix."
        } else {
          sprintf("Gene Matrix needs 2 to 5 saved BED files. Current saved file count: %d.", file_count)
        }
      )
    })

    output$export_panel <- renderUI({
      req(!is.null(analysis_result()))

      react_control_host(
        id = session$ns("export_host"),
        control_type = "gene-matrix-export",
        class = "gene-matrix-export-root",
        ready_input_id = session$ns("export_ready")
      )
    })

    output$parameter_snapshot <- renderUI({
      context <- if (is.null(upload_context)) NULL else upload_context()

      div(
        class = "gene-matrix-summary-grid",
        div(
          class = "gene-matrix-summary-item",
          span(class = "gene-matrix-summary-item__label", "Saved BED Files"),
          span(class = "gene-matrix-summary-item__value", if (is.null(context)) "Not saved yet" else context$file_count)
        ),
        div(
          class = "gene-matrix-summary-item",
          span(class = "gene-matrix-summary-item__label", "Species"),
          span(class = "gene-matrix-summary-item__value", if (is.null(context)) "Not saved yet" else context$species)
        )
      )
    })

    output$analysis_notice <- renderUI({
      if (!is.null(analysis_error())) {
        return(
          div(
            class = "gene-matrix-result-card gene-matrix-result-card--error",
            tags$p(class = "gene-matrix-result-card__copy", analysis_error())
          )
        )
      }

      NULL
    })

    output$analysis_summary <- renderUI({
      result <- analysis_result()
      req(!is.null(result))

      sample_gene_counts <- vapply(result$sample_results, function(item) length(item$gene_set), integer(1))

      div(
        class = "gene-matrix-summary-grid gene-matrix-summary-grid--result",
        div(
          class = "gene-matrix-summary-item",
          span(class = "gene-matrix-summary-item__label", "Species"),
          span(class = "gene-matrix-summary-item__value", result$species)
        ),
        div(
          class = "gene-matrix-summary-item",
          span(class = "gene-matrix-summary-item__label", "Samples"),
          span(class = "gene-matrix-summary-item__value", result$sample_count)
        ),
        div(
          class = "gene-matrix-summary-item",
          span(class = "gene-matrix-summary-item__label", "Union Genes"),
          span(class = "gene-matrix-summary-item__value", result$union_gene_count)
        ),
        div(
          class = "gene-matrix-summary-item",
          span(class = "gene-matrix-summary-item__label", "Per-sample Genes"),
          span(class = "gene-matrix-summary-item__value", paste(sample_gene_counts, collapse = " / "))
        )
      )
    })

    outputOptions(output, "run_hint", suspendWhenHidden = FALSE)
  })
}
