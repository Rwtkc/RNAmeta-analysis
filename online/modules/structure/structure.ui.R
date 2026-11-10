mod_structure_ui <- function(id) {
  ns <- NS(id)

  div(
    class = "page-shell structure-shell",
    tags$script(HTML(sprintf(
      "(function(){
        var nsPrefix = %s;
        if (window.__rnametaStructureBridgeAttached) {
          return;
        }
        window.__rnametaStructureBridgeAttached = true;
        window.__rnametaStructureLatestPayload = null;
        function postStructurePayload(payload) {
          if (!payload || !payload.frameId) {
            return;
          }
          var frame = document.getElementById(payload.frameId);
          if (!frame || !frame.contentWindow) {
            return;
          }
          frame.contentWindow.postMessage({
            type: 'rnameta-structure-render',
            sequence: payload.sequence,
            structure: payload.structure,
            highlights: payload.highlights,
            highlightStart: payload.highlightStart,
            highlightEnd: payload.highlightEnd
          }, '*');
        }
        function scheduleStructureRender(payload) {
          postStructurePayload(payload);
          [150, 500, 1200].forEach(function(delay) {
            window.setTimeout(function() {
              postStructurePayload(payload);
            }, delay);
          });
        }
        Shiny.addCustomMessageHandler('rnameta-structure-render', function(message) {
          window.__rnametaStructureLatestPayload = message;
          scheduleStructureRender(message);
        });
        Shiny.addCustomMessageHandler('rnameta-structure-export', function(message) {
          var frame = document.getElementById(message.frameId);
          if (!frame || !frame.contentWindow) {
            return;
          }
          frame.contentWindow.postMessage({
            type: 'rnameta-structure-export',
            format: message.format,
            width: message.width,
            height: message.height,
            dpi: message.dpi
          }, '*');
        });
        window.addEventListener('message', function(event) {
          var data = event.data || {};
          if (data.type !== 'rnameta-structure-status') {
            return;
          }
          if (window.Shiny && typeof window.Shiny.setInputValue === 'function') {
            window.Shiny.setInputValue(nsPrefix + 'status', data, { priority: 'event' });
          }
          if (data.message === 'Ready.' && window.__rnametaStructureLatestPayload) {
            scheduleStructureRender(window.__rnametaStructureLatestPayload);
          }
        });
      })();",
      jsonlite::toJSON(ns(""), auto_unbox = TRUE)
    ))),
    div(
      class = "structure-workspace",
      tags$aside(
        class = "structure-sidebar",
        div(
          class = "structure-panel structure-panel--sidebar",
          div(class = "structure-panel__eyebrow", "Structure Test"),
          tags$h2(class = "structure-panel__title", "RNA Structure"),
          tags$p(
            class = "structure-panel__copy",
            "Render RNA secondary structures directly from the saved Upload / Run BED context."
          ),
          react_control_host(
            id = ns("controls_host"),
            control_type = "structure-controls",
            class = "structure-controls-root",
            ready_input_id = ns("controls_ready")
          ),
          react_control_host(
            id = ns("export_host"),
            control_type = "structure-export",
            class = "structure-export-root",
            ready_input_id = ns("export_ready")
          )
        )
      ),
      tags$section(
        class = "structure-main",
        div(
          class = "structure-panel structure-panel--intro",
          div(class = "structure-panel__eyebrow", "Structure Workspace"),
          tags$h3(class = "structure-main__title", "BED Structure View"),
          tags$p(
            class = "structure-main__copy",
            "The Structure module uses the saved Upload / Run BED context to extract reference sequence and render RNA secondary structure without manual sequence input."
          ),
          uiOutput(ns("parameter_snapshot"))
        ),
        react_control_host(
          id = ns("mapping_host"),
          control_type = "structure-mapping",
          class = "structure-mapping-root",
          ready_input_id = ns("mapping_ready")
        ),
        div(
          class = "structure-panel structure-panel--canvas",
          div(class = "structure-panel__eyebrow", "Analysis Output"),
          react_control_host(
            id = ns("progress_host"),
            control_type = "structure-progress",
            class = "structure-progress-root",
            ready_input_id = ns("progress_ready")
          ),
          uiOutput(ns("analysis_output"))
        )
      )
    )
  )
}
