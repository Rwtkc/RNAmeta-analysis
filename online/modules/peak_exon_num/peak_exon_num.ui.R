mod_peak_exon_num_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell peak-gene-size-shell",
    div(
      class = "peak-gene-size-workspace",
      tags$aside(
        class = "peak-gene-size-sidebar",
        div(
          class = "peak-gene-size-panel peak-gene-size-panel--sidebar",
          div(class = "peak-gene-size-panel__eyebrow", "Peak Exon Num Controls"),
          tags$h2(class = "peak-gene-size-panel__title", "Peak Exon Num"),
          tags$p(
            class = "peak-gene-size-panel__copy",
            "Summarize exon-count distributions for transcripts connected to the saved Upload / Run BED context."
          ),
          div(
            class = "peak-gene-size-sidebar__action",
            actionButton(
              inputId = ns("run_peak_exon_num"),
              label = "Run Peak Exon Num",
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
          div(class = "peak-gene-size-panel__eyebrow", "Peak Exon Num Workspace"),
          tags$h3(class = "peak-gene-size-main__title", "Exon Count Distribution"),
          tags$p(
            class = "peak-gene-size-main__copy",
            "This module keeps the legacy exon-count summary path while exposing the result inside the refactored React-Shiny shell."
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
            control_type = "peak-exon-num-chart",
            class = "peak-exon-num-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
