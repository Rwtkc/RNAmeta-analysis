mod_motif_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell site-profile-shell motif-shell",
    div(
      class = "site-profile-workspace motif-workspace",
      tags$aside(
        class = "site-profile-sidebar",
        div(
          class = "site-profile-panel site-profile-panel--sidebar motif-panel",
          div(class = "site-profile-panel__eyebrow", "Motif Controls"),
          tags$h2(class = "site-profile-panel__title", "Motif"),
          tags$p(
            class = "site-profile-panel__copy",
            "Submit HOMER, MEME, or STREME motif jobs to the server, keep them running in the background, and look up completed results later with a job id."
          ),
          div(
            class = "motif-retention-tip",
            tags$strong("Result retention"),
            tags$span("Motif job results are retained for 7 days. Save the Job ID to query the result after leaving this page.")
          ),
          div(
            class = "site-profile-sidebar__export",
            uiOutput(ns("export_panel"))
          ),
          div(
            class = "motif-form",
            react_control_host(
              id = ns("controls_host"),
              control_type = "motif-controls",
              class = "motif-controls-root",
              ready_input_id = ns("controls_ready")
            ),
            div(
              class = "motif-native-state",
              selectInput(ns("sample_name"), NULL, choices = character(), selected = ""),
              selectInput(ns("engine"), NULL, choices = motif_engine_choices(), selected = motif_default_controls()$engine),
              selectInput(ns("analysis_mode"), NULL, choices = motif_analysis_mode_choices(), selected = motif_default_controls()$analysis_mode),
              selectInput(ns("sites_subset"), NULL, choices = motif_sites_subset_choices(), selected = motif_default_controls()$sites_subset),
              numericInput(ns("motif_width"), NULL, value = motif_default_controls()$motif_width, min = 1, step = 1),
              numericInput(ns("sequence_window"), NULL, value = motif_default_controls()$sequence_window, min = 1, step = 1),
              numericInput(ns("max_motifs"), NULL, value = motif_default_controls()$max_motifs, min = 1, step = 1),
              textInput(ns("job_lookup"), NULL, value = "")
            ),
            div(
              class = "motif-native-file",
              fileInput(ns("result_bundle"), NULL, accept = c(".json"))
            )
          ),
          div(
            class = "site-profile-sidebar__action",
            actionButton(
              inputId = ns("run_motif"),
              label = "Submit Motif Job",
              class = "meta-btn meta-btn--primary",
              `data-rnameta-analysis-trigger` = "true"
            )
          ),
          uiOutput(ns("run_hint"))
        )
      ),
      tags$section(
        class = "site-profile-main motif-main",
        div(
          class = "site-profile-panel site-profile-panel--intro motif-panel",
          div(class = "site-profile-panel__eyebrow", "Motif Workspace"),
          tags$h3(class = "site-profile-main__title", "Background Motif Jobs"),
          tags$p(
            class = "site-profile-main__copy",
            "Submitted jobs continue running on the server even if you leave this page. Enter a saved job id in the sidebar at any time to reopen the matching result."
          ),
          div(
            class = "motif-retention-tip motif-retention-tip--main",
            tags$strong("Results are available for 7 days."),
            tags$span("Keep the Motif Job ID if you want to reopen completed HOMER, MEME, or STREME results later.")
          )
        ),
        div(
          class = "site-profile-panel site-profile-panel--canvas motif-panel",
          div(class = "site-profile-panel__eyebrow", "Motif Output"),
          uiOutput(ns("analysis_notice")),
          react_control_host(
            id = ns("chart_host"),
            control_type = "motif-result",
            class = "motif-result-root",
            ready_input_id = ns("chart_ready")
          )
        )
      )
    )
  )
}
