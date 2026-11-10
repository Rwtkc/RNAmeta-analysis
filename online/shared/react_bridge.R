react_control_host <- function(id, control_type, class = NULL, ready_input_id = NULL) {
  tags$div(
    id = id,
    class = paste("rnameta-control-host", class),
    `data-rnameta-control` = control_type,
    `data-rnameta-ready-input` = ready_input_id
  )
}

send_control_payload <- function(session, id, control, config = list()) {
  session$sendCustomMessage(
    "rnameta:update-control",
    list(
      id = id,
      control = control,
      config = config
    )
  )
}

clear_control_payload <- function(session, id = NULL) {
  session$sendCustomMessage(
    "rnameta:clear-control",
    list(id = id)
  )
}

set_native_button_state <- function(session, id, disabled = FALSE) {
  session$sendCustomMessage(
    "rnameta:set-button-state",
    list(
      id = id,
      disabled = isTRUE(disabled)
    )
  )
}
