mod_browser_server <- function(id, upload_context = NULL, analysis_lock = NULL) {
  moduleServer(id, function(input, output, session) {
    module_lock_id <- session$ns("analysis_lock")
    browser_cache_scope_id <- browser_cache_scope(paste(
      "browser",
      session$token %||% paste(sample(c(letters, 0:9), size = 12L, replace = TRUE), collapse = ""),
      sep = "_"
    ))
    refresh_nonce <- reactiveVal(0L)
    browser_requested <- reactiveVal(FALSE)
    browser_preparing <- reactiveVal(FALSE)
    browser_prepare_error <- reactiveVal(NULL)
    loaded_upload_context <- reactiveVal(NULL)
    browser_locus <- reactiveVal(NULL)
    browser_table_data <- reactiveVal(NULL)
    browser_table_revision <- reactiveVal(0L)
    browser_last_published_table_revision <- reactiveVal(NULL)
    selected_browser_row_key <- reactiveVal(NULL)
    last_browser_table_click <- reactiveVal(NULL)

    current_lock_owner <- reactive({
      if (is.null(analysis_lock)) {
        return(NULL)
      }

      analysis_lock()
    })

    is_analysis_locked <- reactive({
      !is.null(current_lock_owner())
    })

    increment_table_revision <- function() {
      browser_table_revision(isolate(browser_table_revision()) + 1L)
    }

    selected_browser_row <- reactive({
      row_key <- trimws(as.character(selected_browser_row_key() %||% ""))
      table <- browser_table_data()
      if (!nzchar(row_key) || is.null(table) || !nrow(table)) {
        return(NULL)
      }

      match_index <- which(as.character(table[["row_key"]]) == row_key)
      if (!length(match_index)) {
        return(NULL)
      }

      table[match_index[[1L]]]
    })

    upload_signature <- function(saved_input) {
      files <- upload_run_context_saved_files(saved_input)
      if (!length(files)) {
        return("")
      }

      paste(vapply(files, function(file) {
        path <- as.character(file$file_path %||% "")
        stamp <- if (nzchar(path) && file.exists(path)) as.character(file.info(path)$mtime) else ""
        paste(as.character(file$name %||% ""), path, stamp, sep = "::")
      }, character(1)), collapse = "||")
    }
    loaded_upload_signature <- reactiveVal("")

    has_saved_upload <- reactive({
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      saved_files <- upload_run_context_saved_files(saved_input)
      !is.null(saved_input) &&
        length(saved_files) > 0 &&
        all(vapply(saved_files, function(file) {
          !is.null(file$file_path) && nzchar(file$file_path) && file.exists(file$file_path)
        }, logical(1)))
    })

    reset_browser_state <- function() {
      browser_requested(FALSE)
      browser_preparing(FALSE)
      browser_prepare_error(NULL)
      loaded_upload_context(NULL)
      browser_locus(NULL)
      browser_table_data(NULL)
      selected_browser_row_key(NULL)
      last_browser_table_click(NULL)
      increment_table_revision()
    }

    session$onSessionEnded(function() {
      browser_cleanup_cache_scope(browser_cache_scope_id)
    })

    build_browser_config <- function(saved_input, use_filtered_cache = isTRUE(browser_requested()), selected_transcript = NULL) {
      config <- browser_build_config(
        saved_input,
        "transcript",
        use_filtered_cache = use_filtered_cache,
        cache_scope = browser_cache_scope_id,
        selected_transcript = selected_transcript
      )
      config$assetKey <- paste(config$assetKey %||% config$status, refresh_nonce(), sep = "::refresh::")
      config
    }

    current_config <- reactive({
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      build_browser_config(saved_input, selected_transcript = selected_browser_row())
    })

    loaded_config <- reactive({
      saved_input <- loaded_upload_context()
      if (is.null(saved_input)) {
        return(NULL)
      }
      build_browser_config(
        saved_input,
        use_filtered_cache = TRUE,
        selected_transcript = selected_browser_row()
      )
    })

    observeEvent(if (is.null(upload_context)) NULL else upload_context(), {
      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      next_signature <- upload_signature(saved_input)
      if (!identical(next_signature, loaded_upload_signature())) {
        reset_browser_state()
        loaded_upload_signature(next_signature)
      }
    }, ignoreInit = TRUE)

    observe({
      if (!has_saved_upload()) {
        reset_browser_state()
        loaded_upload_signature("")
      }
    })

    observeEvent(input$refresh_browser, {
      req(has_saved_upload())
      req(!is_analysis_locked())

      if (!is.null(analysis_lock)) {
        analysis_lock(module_lock_id)
      }
      on.exit({
        if (!is.null(analysis_lock) && identical(analysis_lock(), module_lock_id)) {
          analysis_lock(NULL)
        }
      }, add = TRUE)

      saved_input <- if (is.null(upload_context)) NULL else upload_context()
      browser_preparing(TRUE)
      browser_requested(FALSE)
      browser_prepare_error(NULL)
      browser_table_data(NULL)
      selected_browser_row_key(NULL)
      browser_locus(NULL)
      last_browser_table_click(NULL)
      increment_table_revision()

      is_ready <- tryCatch(
        browser_build_transcript_table(
          upload_context = saved_input,
          cache_scope = browser_cache_scope_id,
          progress = NULL
        ),
        error = function(err) {
          browser_prepare_error(conditionMessage(err))
          NULL
        }
      )

      browser_preparing(FALSE)
      if (is.null(is_ready)) {
        return()
      }

      refresh_nonce(refresh_nonce() + 1L)
      loaded_upload_context(saved_input)
      loaded_upload_signature(upload_signature(saved_input))
      browser_table_data(is_ready)
      increment_table_revision()
      first_mapped <- is_ready[nzchar(transcript_id %||% "") & nzchar(igv_locus %||% "")]
      if (nrow(first_mapped)) {
        selected_browser_row_key(first_mapped$row_key[[1]])
        browser_locus(first_mapped$igv_locus[[1]])
      }
      browser_requested(TRUE)
    }, ignoreInit = TRUE)

    select_browser_table_row <- function(row_key) {
      row_key <- trimws(as.character(row_key %||% ""))
      table <- browser_table_data()
      if (!nzchar(row_key) || is.null(table) || !nrow(table)) {
        return(invisible(NULL))
      }

      match_index <- which(as.character(table[["row_key"]]) == row_key)
      if (!length(match_index)) {
        return(invisible(NULL))
      }

      match <- table[match_index[[1L]]]
      selected_browser_row_key(match$row_key[[1]])
      browser_locus(match$igv_locus[[1]])
      refresh_nonce(refresh_nonce() + 1L)
      invisible(NULL)
    }

    observeEvent(input$browser_table_click, {
      row_key <- trimws(as.character(input$browser_table_click %||% ""))
      if (nzchar(row_key) && !identical(row_key, last_browser_table_click())) {
        last_browser_table_click(row_key)
        select_browser_table_row(row_key)
      }
    }, ignoreInit = TRUE)

    browser_status_label <- function(status) {
      switch(
        status,
        ready = "Ready",
        reference_only = "Reference only",
        missing_reference = "Missing reference",
        missing_annotation = "Missing annotation",
        unsupported_species = "Waiting for saved BED input",
        "Pending"
      )
    }

    browser_react_config <- reactive({
      requested <- isTRUE(browser_requested())
      config <- if (requested && !is.null(loaded_upload_context())) loaded_config() else current_config()
      selected <- selected_browser_row()
      workspace_status <- if (isTRUE(browser_preparing())) {
        "preparing"
      } else if (!is.null(browser_prepare_error())) {
        "error"
      } else if (requested) {
        "ready"
      } else {
        "idle"
      }
      message <- if (identical(workspace_status, "error")) {
        browser_prepare_error()
      } else if (identical(workspace_status, "idle") && !has_saved_upload()) {
        "Save a BED file in Upload / Run first to enable Genome Browser."
      } else if (identical(workspace_status, "idle")) {
        "Refresh Browser to map the saved BED file to transcript coordinates."
      } else if (is.null(selected)) {
        "No mapped transcript is available for IGV."
      } else {
        ""
      }
      tracks <- if (requested && !is.null(config) && !is.null(selected)) config$tracks %||% list() else list()

      list(
        workspaceStatus = workspace_status,
        message = as.character(message %||% ""),
        hasSavedUpload = isTRUE(has_saved_upload()),
        refreshInput = session$ns("refresh_browser"),
        rowClickInput = session$ns("browser_table_click"),
        tableRevision = browser_table_revision(),
        tablePageSize = 10L,
        selectedRowKey = as.character(selected_browser_row_key() %||% ""),
        summary = list(
          reference = config$speciesLabel %||% config$species %||% "Not selected",
          view = "Transcript table",
          locus = browser_locus() %||% config$locus %||% "Not selected",
          tracks = as.character(length(tracks)),
          status = browser_status_label(config$status %||% "pending")
        ),
        reference = if (requested && !is.null(selected)) config$reference else NULL,
        locus = if (requested && !is.null(selected)) config$locus else NULL,
        tracks = tracks,
        assetKey = if (requested && !is.null(selected)) config$assetKey else NULL
      )
    })

    publish_react_browser <- function(force_table = FALSE) {
      config <- browser_react_config()
      table_revision <- browser_table_revision()
      include_table <- isTRUE(force_table) || !identical(table_revision, isolate(browser_last_published_table_revision()))
      if (include_table) {
        config$table <- browser_react_transcript_table_payload(browser_table_data())
        browser_last_published_table_revision(table_revision)
      }

      send_control_payload(
        session = session,
        id = session$ns("browser_react_view"),
        control = "genome-browser",
        config = config
      )
    }

    observeEvent(input$browser_react_ready, {
      browser_last_published_table_revision(NULL)
      publish_react_browser(force_table = TRUE)
    }, ignoreInit = TRUE)

    observe({
      browser_react_config()
      publish_react_browser()
    })
  })
}
