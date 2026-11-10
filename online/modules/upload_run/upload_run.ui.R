mod_upload_run_ui <- function(id) {
  ns <- NS(id)

  tagList(
    div(
      class = "page-shell upload-shell",
      div(
        class = "upload-workspace",
        tags$aside(
          class = "upload-sidebar",
          div(
            class = "upload-panel upload-panel--sidebar",
            div(class = "upload-panel__eyebrow", "Choose Data Source"),
            tags$h2(class = "upload-panel__title", "Upload / Run"),
            tags$p(
              class = "upload-panel__copy",
              "Choose how RNAmeta should load your data, then set species context and provide BED files for downstream analysis."
            ),
            div(
              class = "upload-form",
              react_control_host(
                id = ns("upload_controls"),
                control_type = "upload-sidebar",
                class = "upload-controls-root",
                ready_input_id = ns("controls_ready")
              ),
              div(
                class = "upload-native-state",
                textInput(
                  inputId = ns("data_source"),
                  label = NULL,
                  value = "Upload File"
                ),
                textInput(
                  inputId = ns("species"),
                  label = NULL,
                  value = ""
                ),
                div(
                  class = "upload-native-save",
                  actionButton(
                    inputId = ns("save_upload"),
                    label = "Save"
                  ),
                  actionButton(
                    inputId = ns("load_demo"),
                    label = "Demo"
                  )
                )
              ),
              div(
                class = "upload-native-file",
                fileInput(
                  inputId = ns("bed_file"),
                  label = NULL,
                  multiple = TRUE,
                  accept = c(".bed", ".BED")
                )
              )
            )
          )
        ),
        tags$section(
          class = "upload-main",
          div(
            class = "upload-panel upload-panel--summary",
            div(class = "upload-panel__eyebrow", "Session Summary"),
            uiOutput(ns("session_summary"))
          ),
          react_control_host(
            id = ns("database_browser"),
            control_type = "upload-database-browser",
            class = "upload-database-host"
          ),
          uiOutput(ns("workflow_notes"))
        )
      )
    )
  )
}
