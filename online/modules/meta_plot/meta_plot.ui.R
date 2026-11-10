mod_meta_plot_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell meta-plot-shell",
    div(
      class = "meta-plot-workspace",
      tags$aside(
        class = "meta-plot-sidebar",
        div(
          class = "meta-plot-panel meta-plot-panel--sidebar",
          div(class = "meta-plot-panel__eyebrow", "Meta Plot Controls"),
          tags$h2(class = "meta-plot-panel__title", "Meta Plot"),
          tags$p(
            class = "meta-plot-panel__copy",
            "Configure transcript filters, sampling rules, and smoothing options before rendering transcript-scale metagene distributions."
          ),
          div(
            class = "meta-plot-sidebar__action",
            actionButton(
              inputId = ns("run_plot"),
              label = "Run Meta Plot Analysis",
              class = "meta-btn meta-btn--primary",
              `data-rnameta-analysis-trigger` = "true",
              `data-rnameta-analysis-owner` = ns("progress_slot")
            )
          ),
          uiOutput(ns("run_hint")),
          div(
            class = "meta-plot-sidebar__export",
            uiOutput(ns("export_panel"))
          ),
          div(
            class = "meta-plot-form",
            react_control_host(
              id = ns("controls_host"),
              control_type = "meta-plot-controls",
              class = "meta-plot-controls-root",
              ready_input_id = ns("controls_ready")
            )
          )
        )
      ),
      tags$section(
        class = "meta-plot-main",
        div(
          class = "meta-plot-panel meta-plot-panel--intro",
          div(class = "meta-plot-panel__eyebrow", "Meta Plot Workspace"),
          tags$h3(class = "meta-plot-main__title", "Meta Gene Stat"),
          tags$p(
            class = "meta-plot-main__copy",
            "The refactored Meta Plot module preserves the original parameter surface while reorganizing it into a cleaner workspace for transcript-scale distribution analysis."
          ),
          uiOutput(ns("parameter_snapshot"))
        ),
        div(
          class = "meta-plot-panel meta-plot-panel--canvas",
          div(class = "meta-plot-panel__eyebrow", "Analysis Flow"),
          div(
            class = "meta-plot-progress-slot",
            `data-meta-plot-progress-slot` = "true",
            `data-rnameta-progress-slot` = ns("progress_slot")
          ),
          uiOutput(ns("analysis_notice")),
          uiOutput(ns("analysis_summary")),
          react_control_host(
            id = ns("chart_host"),
            control_type = "meta-plot-chart",
            class = "meta-plot-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
