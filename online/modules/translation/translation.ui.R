mod_translation_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell site-profile-shell",
    div(
      class = "site-profile-workspace",
      tags$aside(
        class = "site-profile-sidebar",
        div(
          class = "site-profile-panel site-profile-panel--sidebar",
          div(class = "site-profile-panel__eyebrow", "Translation Controls"),
          tags$h2(class = "site-profile-panel__title", "Translation"),
          tags$p(
            class = "site-profile-panel__copy",
            "Summarize peak density around translation start and end boundaries and keep the legacy first-sample heatmap behavior inside the refactored workflow."
          ),
          div(
            class = "site-profile-sidebar__action",
            actionButton(
              inputId = ns("run_translation"),
              label = "Run Translation",
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
          div(class = "site-profile-panel__eyebrow", "Translation Workspace"),
          tags$h3(class = "site-profile-main__title", "Boundary Density and Heatmap View"),
          tags$p(
            class = "site-profile-main__copy",
            "This module reuses the saved Upload / Run context and preserves the old translation start and end visualization structure while shifting rendering to the bridge-driven frontend."
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
            control_type = "translation-chart",
            class = "site-profile-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
