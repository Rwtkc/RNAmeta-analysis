mod_peak_gene_size_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell peak-gene-size-shell",
    div(
      class = "peak-gene-size-workspace",
      tags$aside(
        class = "peak-gene-size-sidebar",
        div(
          class = "peak-gene-size-panel peak-gene-size-panel--sidebar",
          div(class = "peak-gene-size-panel__eyebrow", "Peak Gene Size Controls"),
          tags$h2(class = "peak-gene-size-panel__title", "Peak Gene Size"),
          tags$p(
            class = "peak-gene-size-panel__copy",
            "Summarize transcript lengths associated with the saved Upload / Run BED context and reproduce the legacy Peak Gene Size distribution path in an isolated module."
          ),
          div(
            class = "peak-gene-size-sidebar__action",
            actionButton(
              inputId = ns("run_peak_gene_size"),
              label = "Run Peak Gene Size",
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
        class = "peak-gene-size-main",
        div(
          class = "peak-gene-size-panel peak-gene-size-panel--intro",
          div(class = "peak-gene-size-panel__eyebrow", "Peak Gene Size Workspace"),
          tags$h3(class = "peak-gene-size-main__title", "Transcript Length Distribution"),
          tags$p(
            class = "peak-gene-size-main__copy",
            "This module keeps the legacy transcript-length extraction logic while lifting Peak Gene Size into the refactored React-Shiny workflow."
          ),
          uiOutput(ns("parameter_snapshot"))
        ),
        div(
          class = "peak-gene-size-panel peak-gene-size-panel--canvas",
          div(class = "peak-gene-size-panel__eyebrow", "Analysis Output"),
          div(class = "peak-gene-size-progress-slot", `data-rnameta-progress-slot` = ns("progress_slot")),
          uiOutput(ns("analysis_notice")),
          uiOutput(ns("analysis_summary")),
          react_control_host(
            id = ns("chart_host"),
            control_type = "peak-gene-size-chart",
            class = "peak-gene-size-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
