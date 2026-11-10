source(app_path("modules", "help", "help.content.R"), local = environment(), encoding = "UTF-8")

help_list_card_ui <- function(list_definition) {
  div(
    class = "help-list-card",
    tags$strong(class = "help-list-card__title", list_definition$label),
    tags$ul(
      class = "help-list-card__list",
      lapply(list_definition$items, function(item) tags$li(item))
    )
  )
}

help_block_ui <- function(block, block_id, figure_slots = list()) {
  tags$section(
    id = block_id,
    class = "help-doc-block",
    tags$h3(class = "help-section-card__title", block$title),
    tagList(lapply(block$paragraphs, function(paragraph) tags$p(class = "help-section-card__copy", paragraph))),
    if (length(figure_slots)) tagList(lapply(figure_slots, help_figure_slot_ui)),
    if (length(block$lists)) {
      div(
        class = "help-list-grid",
        lapply(block$lists, help_list_card_ui)
      )
    }
  )
}

help_nav_group_ui <- function(ns, section) {
  topic_id <- ns(section$id)
  div(
    class = "help-nav__group",
    `data-help-nav-group` = topic_id,
    `data-help-search` = paste(unlist(section, use.names = FALSE), collapse = " "),
    tags$button(
      type = "button",
      class = "help-nav__toggle",
      `data-help-topic-trigger` = topic_id,
      `aria-controls` = paste0(topic_id, "-children"),
      `aria-expanded` = "false",
      tags$span(section$label),
      tags$span(class = "help-nav__chevron", `aria-hidden` = TRUE)
    ),
    div(
      id = paste0(topic_id, "-children"),
      class = "help-nav__children",
      lapply(seq_along(section$blocks), function(index) {
        tags$a(
          href = paste0("#", ns(paste0(section$id, "-block-", index))),
          `data-help-topic-target` = topic_id,
          section$blocks[[index]]$title
        )
      })
    )
  )
}

help_figure_slots <- function() {
  list(
    "upload-run" = list(
      "1" = list(list("upload-run-controls", "Choose the data source and species", "Data-source selector, species selector, and BED-file picker.")),
      "2" = list(list("upload-run-multiple-files", "Save multiple BED files", "Saved-file chips or selectors showing more than one input file.")),
      "3" = list(list("upload-run-preview", "Inspect the saved-file preview", "Session Summary and the first rows of the saved BED file."))
    ),
    "peak-distribution" = list(
      "2" = list(
        list("peak-distribution-controls", "Configure visible feature classes", "Sidebar display controls before running the analysis."),
        list("peak-distribution-result", "Inspect and export the distribution", "Feature-category result, active legend, and export panel.")
      )
    ),
    "genome-browser" = list(
      "1" = list(list("genome-browser-view", "Inspect a genomic locus", "Loaded IGV-Web view with the saved BED track and reference annotation.")),
      "2" = list(list("genome-browser-refresh", "Refresh after changing input", "Refresh Browser action and the active session snapshot."))
    ),
    "structure" = list(
      "1" = list(list("structure-mapping-table", "Map BED intervals to transcripts", "Transcript Mapping Table after Run Structure.")),
      "2" = list(list("structure-rnafold-submit", "Submit a row to RNAfold", "Selected mapping row and RNAfold confirmation dialog.")),
      "3" = list(list("structure-result", "Review a completed structure job", "Job status, RNA secondary structure, transcript track, and export controls."))
    ),
    "gene-statistics" = list(
      "1" = list(list("gene-statistics-gene-type", "Gene Type result", "Gene Type controls, result, and export panel.")),
      "2" = list(list("gene-statistics-peak-gene-size", "Peak Gene Size result", "Peak Gene Size controls, result, and export panel.")),
      "3" = list(list("gene-statistics-gene-matrix", "Gene Matrix result", "Gene Matrix input selection, overlap result, and export panel."))
    ),
    "exon-statistics" = list(
      "1" = list(list("exon-statistics-size", "Peak Exon Size result", "Peak Exon Size controls, result, and export panel.")),
      "2" = list(list("exon-statistics-type", "Peak Exon Type result", "Peak Exon Type controls, result, and export panel.")),
      "3" = list(list("exon-statistics-num", "Peak Exon Num result", "Peak Exon Num controls, result, and export panel."))
    ),
    "site" = list(
      "1" = list(list("site-transcription", "Transcription profile", "Transcription controls and completed profile.")),
      "2" = list(list("site-translation", "Translation profile", "Translation controls and completed profile.")),
      "3" = list(list("site-splicesite", "Splicesite profile", "Splicesite controls and completed profile."))
    ),
    "motif" = list(
      "1" = list(list("motif-submit", "Submit a motif job", "Motif discovery settings and job submission controls.")),
      "2" = list(list("motif-status", "Track motif job status", "Queue status and Motif Job ID lookup.")),
      "3" = list(list("motif-result", "Inspect and export motifs", "Completed motif result and export panel."))
    )
  )
}

help_figure_slot_ui <- function(slot) {
  tags$figure(
    id = paste0("figure-", slot[[1]]),
    class = "help-figure-slot",
    `data-help-figure-slot` = slot[[1]],
    div(class = "help-figure-slot__label", "Screenshot placeholder"),
    div(class = "help-figure-slot__frame", tags$span("Add screenshot")),
    tags$figcaption(
      tags$strong(slot[[2]]),
      tags$span(slot[[3]])
    )
  )
}

help_topic_ui <- function(ns, section) {
  figure_slots <- help_figure_slots()[[section$id]]
  topic_id <- ns(section$id)
  tags$section(
    id = topic_id,
    class = "help-topic",
    `data-help-topic` = TRUE,
    `data-help-search` = paste(unlist(section, use.names = FALSE), collapse = " "),
    div(
      class = "help-topic__head",
      tags$h2(class = "help-topic__title", section$label),
      tags$p(class = "help-topic__summary", section$summary)
    ),
    div(
      class = "help-topic__stack",
      lapply(seq_along(section$blocks), function(index) {
        help_block_ui(
          section$blocks[[index]],
          ns(paste0(section$id, "-block-", index)),
          figure_slots[[as.character(index)]]
        )
      })
    )
  )
}

mod_help_ui <- function(id) {
  ns <- NS(id)
  sections <- help_sections()
  search_id <- ns("doc_search")

  div(
    class = "page-shell help-shell",
    div(
      class = "help-workspace",
      tags$aside(
        class = "help-sidebar",
        div(
          class = "help-sidebar__brand",
          tags$span("RNAmeta Online"),
          tags$strong("Documentation")
        ),
        tags$label(
          class = "help-search",
          `for` = search_id,
          tags$span("Search docs"),
          tags$input(
            id = search_id,
            class = "help-search__input",
            type = "search",
            placeholder = "Search documentation",
            autocomplete = "off"
          )
        ),
        tags$nav(
          class = "help-nav",
          `aria-label` = "Documentation navigation",
          lapply(sections, function(section) help_nav_group_ui(ns, section))
        )
      ),
      tags$section(
        class = "help-main",
        div(
          class = "help-doc-header",
          tags$div(class = "help-breadcrumb", tags$span("Docs"), tags$span("RNAmeta Online")),
          tags$h1(class = "help-hero__title", "RNAmeta Online documentation"),
          tags$p(
            class = "help-hero__copy",
            "A practical reference for preparing BED input, running RNAmeta Online modules, interpreting results, and exporting analysis outputs."
          )
        ),
        lapply(sections, function(section) help_topic_ui(ns, section)),
        tags$script(src = "js/help-search.js?v=help-search-20260822a")
      ),
      tags$script(HTML(sprintf("(function () { const input = document.getElementById('%s'); if (!input) return; const root = input.closest('.help-shell'); const topics = Array.from(root.querySelectorAll('[data-help-topic]')); const groups = Array.from(root.querySelectorAll('[data-help-nav-group]')); const normalize = (value) => String(value || '').toLowerCase(); let activeId = topics[0] ? topics[0].id : ''; const setActive = (topicId, scroll) => { if (!topicId) return; activeId = topicId; topics.forEach((topic) => { topic.hidden = topic.id !== topicId; }); groups.forEach((group) => { const active = group.dataset.helpNavGroup === topicId; group.classList.toggle('is-active', active); const toggle = group.querySelector('[data-help-topic-trigger]'); if (toggle) toggle.setAttribute('aria-expanded', String(active)); }); if (scroll) document.getElementById(topicId)?.scrollIntoView({ behavior: 'smooth', block: 'start' }); }; root.querySelectorAll('[data-help-topic-trigger]').forEach((trigger) => { trigger.addEventListener('click', () => setActive(trigger.dataset.helpTopicTrigger, true)); }); root.querySelectorAll('[data-help-topic-target]').forEach((link) => { link.addEventListener('click', () => setActive(link.dataset.helpTopicTarget, false)); }); const filter = () => { const query = normalize(input.value).trim(); const matches = groups.filter((group) => !query || normalize(group.dataset.helpSearch).includes(query)); groups.forEach((group) => { group.hidden = !matches.includes(group); }); if (matches.length && !matches.some((group) => group.dataset.helpNavGroup === activeId)) setActive(matches[0].dataset.helpNavGroup, false); }; input.addEventListener('input', filter); setActive(activeId, false); })();", search_id)))
    )
  )
}
