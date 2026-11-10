mod_gene_matrix_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell gene-matrix-shell",
    div(
      class = "gene-matrix-workspace",
      tags$aside(
        class = "gene-matrix-sidebar",
        div(
          class = "gene-matrix-panel gene-matrix-panel--sidebar",
          div(class = "gene-matrix-panel__eyebrow", "Gene Matrix Controls"),
          tags$h2(class = "gene-matrix-panel__title", "Gene Matrix"),
          tags$p(
            class = "gene-matrix-panel__copy",
            "Compare shared gene coverage across 2 to 5 saved BED files with the current intersection matrix view inside the refactored workflow."
          ),
          div(
            class = "gene-matrix-sidebar__action",
            actionButton(
              inputId = ns("run_gene_matrix"),
              label = "Run Gene Matrix",
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
        class = "gene-matrix-main",
        div(
          class = "gene-matrix-panel gene-matrix-panel--intro",
          div(class = "gene-matrix-panel__eyebrow", "Gene Matrix Workspace"),
          tags$h3(class = "gene-matrix-main__title", "Shared Gene Overlap Summary"),
          tags$p(
            class = "gene-matrix-main__copy",
            "This module isolates the legacy Gene Matrix execution chain so multi-sample overlap views can be reproduced without reusing the old monolithic pipeline."
          ),
          uiOutput(ns("parameter_snapshot"))
        ),
        div(
          class = "gene-matrix-panel gene-matrix-panel--canvas",
          div(class = "gene-matrix-panel__eyebrow", "Analysis Output"),
          div(class = "gene-matrix-progress-slot", `data-rnameta-progress-slot` = ns("progress_slot")),
          uiOutput(ns("analysis_notice")),
          uiOutput(ns("analysis_summary")),
          react_control_host(
            id = ns("chart_host"),
            control_type = "gene-matrix-chart",
            class = "gene-matrix-chart-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
