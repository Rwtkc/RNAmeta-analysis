mod_transcription_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell site-profile-shell",
    div(
      class = "site-profile-workspace",
      tags$aside(
        class = "site-profile-sidebar",
        div(
          class = "site-profile-panel site-profile-panel--sidebar",
          div(class = "site-profile-panel__eyebrow", "Transcription Controls"),
          tags$h2(class = "site-profile-panel__title", "Transcription"),
          tags$p(
            class = "site-profile-panel__copy",
            "Profile peak density around transcription start and end boundaries while preserving the legacy first-sample heatmap view inside the refactored workflow."
          ),
          div(
            class = "site-profile-sidebar__action",
            actionButton(
              inputId = ns("run_transcription"),
              label = "Run Transcription",
              class = "meta-btn meta-btn--primary",
              `data-rnameta-analysis-trigger` = "true",
              `data-rnameta-analysis-owner` = ns("progress_slot")
            )
          ),
          uiOutput(ns("run_hint")),
          uiOutput(ns("export_panel"))
        )
      ),
      tags$section(
        class = "site-profile-main",
        div(
          class = "site-profile-panel site-profile-panel--intro",
          div(class = "site-profile-panel__eyebrow", "Transcription Workspace"),
          tags$h3(class = "site-profile-main__title", "Boundary Density and Heatmap View"),
          tags$p(
            class = "site-profile-main__copy",
            "This module keeps the legacy transcription boundary workflow, including the first-sample heatmap behavior, while routing the output through the current React-Shiny bridge."
          ),
          uiOutput(ns("parameter_snapshot"))
        ),
        div(
          class = "site-profile-panel site-profile-panel--canvas",
          div(class = "site-profile-panel__eyebrow", "Analysis Output"),
          div(class = "site-profile-progress-slot", `data-rnameta-progress-slot` = ns("progress_slot")),
          uiOutput(ns("analysis_notice")),
          uiOutput(ns("analysis_summary")),
          react_control_host(
            id = ns("chart_host"),
            control_type = "transcription-chart",
            class = "site-profile-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
