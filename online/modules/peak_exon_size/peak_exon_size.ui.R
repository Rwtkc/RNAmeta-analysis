mod_peak_exon_size_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell peak-gene-size-shell",
    div(
      class = "peak-gene-size-workspace",
      tags$aside(
        class = "peak-gene-size-sidebar",
        div(
          class = "peak-gene-size-panel peak-gene-size-panel--sidebar",
          div(class = "peak-gene-size-panel__eyebrow", "Peak Exon Size Controls"),
          tags$h2(class = "peak-gene-size-panel__title", "Peak Exon Size"),
          tags$p(
            class = "peak-gene-size-panel__copy",
            "Summarize exon-length distributions linked to the saved Upload / Run BED context while preserving the legacy first, middle, and last exon split."
          ),
          div(
            class = "peak-gene-size-sidebar__action",
            actionButton(
              inputId = ns("run_peak_exon_size"),
              label = "Run Peak Exon Size",
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
          div(class = "peak-gene-size-panel__eyebrow", "Peak Exon Size Workspace"),
          tags$h3(class = "peak-gene-size-main__title", "Exon Length Distribution"),
          tags$p(
            class = "peak-gene-size-main__copy",
            "This module keeps the legacy exon-length extraction path while lifting Peak Exon Size into the refactored React-Shiny workflow."
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
            control_type = "peak-exon-size-chart",
            class = "peak-exon-size-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
