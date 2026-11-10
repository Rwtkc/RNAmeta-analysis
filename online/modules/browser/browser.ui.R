mod_browser_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell browser-shell",
    react_control_host(
      id = ns("browser_react_view"),
      control_type = "genome-browser",
      class = "browser-react-view",
      ready_input_id = ns("browser_react_ready")
    )
  )
}
