mod_gene_type_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell gene-type-shell",
    div(
      class = "gene-type-workspace",
      tags$aside(
        class = "gene-type-sidebar",
        div(
          class = "gene-type-panel gene-type-panel--sidebar",
          div(class = "gene-type-panel__eyebrow", "Gene Type Controls"),
          tags$h2(class = "gene-type-panel__title", "Gene Type"),
          tags$p(
            class = "gene-type-panel__copy",
            "Summarize transcript biotypes by reusing the saved Upload / Run BED context while keeping the legacy overlap-counting path isolated from unrelated feature statistics."
          ),
          div(
            class = "gene-type-sidebar__action",
            actionButton(
              inputId = ns("run_gene_type"),
              label = "Run Gene Type",
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
          div(class = "gene-type-panel__eyebrow", "Gene Type Workspace"),
          tags$h3(class = "gene-type-main__title", "Transcript Biotype Summary"),
          tags$p(
            class = "gene-type-main__copy",
            "This module extracts the legacy Gene Type execution path into its own analysis flow so the result stays aligned while the code becomes easier to maintain."
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
            control_type = "gene-type-chart",
            class = "gene-type-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
