mod_welcome_ui <- function(id) {
  ns <- NS(id)
  shell_config <- welcome_shell_json(welcome_shell_config())

  tagList(
    tags$section(
      class = "welcome-module",
      div(
        id = ns("app_shell_root"),
        class = "app-shell-root",
        `data-shell-config` = shell_config
      )
    )
  )
}
