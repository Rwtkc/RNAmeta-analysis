mod_splicesite_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell site-profile-shell",
    div(
      class = "site-profile-workspace",
      tags$aside(
        class = "site-profile-sidebar",
        div(
          class = "site-profile-panel site-profile-panel--sidebar",
          div(class = "site-profile-panel__eyebrow", "Splicesite Controls"),
          tags$h2(class = "site-profile-panel__title", "Splicesite"),
          tags$p(
            class = "site-profile-panel__copy",
            "Measure peak density around 5' and 3' splice-site neighborhoods while keeping the refactored workflow tied to the shared Upload / Run context."
          ),
          div(
            class = "site-profile-sidebar__action",
            actionButton(
              inputId = ns("run_splicesite"),
              label = "Run Splicesite",
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
          div(class = "site-profile-panel__eyebrow", "Splicesite Workspace"),
          tags$h3(class = "site-profile-main__title", "Splice-Site Density View"),
          tags$p(
            class = "site-profile-main__copy",
            "This module lifts the legacy splice-site neighborhood density plots into the same bridge-driven export and multi-file workflow used by the stable modules."
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
            control_type = "splicesite-chart",
            class = "site-profile-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
