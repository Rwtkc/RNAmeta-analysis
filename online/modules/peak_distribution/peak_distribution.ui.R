mod_peak_distribution_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell peak-distribution-shell",
    div(
      class = "peak-distribution-workspace",
      tags$aside(
        class = "peak-distribution-sidebar",
        div(
          class = "peak-distribution-panel peak-distribution-panel--sidebar",
          div(class = "peak-distribution-panel__eyebrow", "Peak Distribution Controls"),
          tags$h2(class = "peak-distribution-panel__title", "Peak Distribution"),
          tags$p(
            class = "peak-distribution-panel__copy",
            "Run feature-level peak annotation against the saved Upload / Run BED context and compare how intervals distribute across promoter, UTR, CDS, intron, and intergenic regions."
          ),
          div(
            class = "peak-distribution-sidebar__action",
            actionButton(
              inputId = ns("run_peak_distribution"),
              label = "Run Peak Distribution",
              class = "meta-btn meta-btn--primary",
              `data-rnameta-analysis-trigger` = "true",
              `data-rnameta-analysis-owner` = ns("progress_slot")
            )
          ),
          uiOutput(ns("run_hint")),
          uiOutput(ns("export_panel")),
          div(
            class = "peak-distribution-form",
            div(
              class = "peak-distribution-card",
              tags$h3(class = "peak-distribution-card__title", "Display Settings"),
              react_control_host(
                id = ns("controls_host"),
                control_type = "peak-distribution-controls",
                class = "peak-distribution-controls-root",
                ready_input_id = ns("controls_ready")
              ),
              div(
                class = "peak-distribution-native-state",
                checkboxGroupInput(
                  inputId = ns("selected_features"),
                  label = "Visible Features",
                  choices = peak_distribution_feature_order(),
                  selected = peak_distribution_feature_order()
                )
              )
            )
          )
        )
      ),
      tags$section(
        class = "peak-distribution-main",
        div(
          class = "peak-distribution-panel peak-distribution-panel--intro",
          div(class = "peak-distribution-panel__eyebrow", "Peak Distribution Workspace"),
          tags$h3(class = "peak-distribution-main__title", "Feature-Level Peak Annotation"),
          tags$p(
            class = "peak-distribution-main__copy",
            "This refactored module keeps the legacy Peak Distribution annotation logic while isolating it from unrelated analysis branches, so accuracy stays aligned and the execution path stays readable."
          ),
          uiOutput(ns("parameter_snapshot"))
        ),
        div(
          class = "peak-distribution-panel peak-distribution-panel--canvas",
          div(class = "peak-distribution-panel__eyebrow", "Analysis Output"),
          div(class = "peak-distribution-progress-slot", `data-rnameta-progress-slot` = ns("progress_slot")),
          uiOutput(ns("analysis_notice")),
          uiOutput(ns("analysis_summary")),
          uiOutput(ns("chart_empty_state")),
          react_control_host(
            id = ns("chart_host"),
            control_type = "peak-distribution-chart",
            class = "peak-distribution-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
