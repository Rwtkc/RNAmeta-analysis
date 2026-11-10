mod_peak_exon_type_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell gene-type-shell",
    div(
      class = "gene-type-workspace",
      tags$aside(
        class = "gene-type-sidebar",
        div(
          class = "gene-type-panel gene-type-panel--sidebar",
          div(class = "gene-type-panel__eyebrow", "Peak Exon Type Controls"),
          tags$h2(class = "gene-type-panel__title", "Peak Exon Type"),
          tags$p(
            class = "gene-type-panel__copy",
            "Compare the relative share of first, middle, and last exon overlaps across the saved Upload / Run BED context."
          ),
          div(
            class = "gene-type-sidebar__action",
            actionButton(
              inputId = ns("run_peak_exon_type"),
              label = "Run Peak Exon Type",
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
        class = "gene-type-main",
        div(
          class = "gene-type-panel gene-type-panel--intro",
          div(class = "gene-type-panel__eyebrow", "Peak Exon Type Workspace"),
          tags$h3(class = "gene-type-main__title", "Exon Position Composition"),
          tags$p(
            class = "gene-type-main__copy",
            "This module preserves the legacy exon-position percentage summary while moving the visible chart into the refactored React-Shiny shell."
          ),
          uiOutput(ns("parameter_snapshot"))
        ),
        div(
          class = "gene-type-panel gene-type-panel--canvas",
          div(class = "gene-type-panel__eyebrow", "Analysis Output"),
          div(class = "gene-type-progress-slot", `data-rnameta-progress-slot` = ns("progress_slot")),
          uiOutput(ns("analysis_notice")),
          uiOutput(ns("analysis_summary")),
          react_control_host(
            id = ns("chart_host"),
            control_type = "peak-exon-type-chart",
            class = "peak-exon-type-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
