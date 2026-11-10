source(app_path("shared", "react_bridge.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "welcome", "welcome.adapter.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "welcome", "welcome.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "help", "help.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "upload_run", "upload_run.shared.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "upload_run", "upload_run.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "browser", "browser.functions.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "browser", "browser.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "structure", "structure.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "gene_type", "gene_type.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "gene_matrix", "gene_matrix.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "peak_gene_size", "peak_gene_size.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "peak_exon_size", "peak_exon_size.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "peak_exon_type", "peak_exon_type.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "peak_exon_num", "peak_exon_num.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "peak_distribution", "peak_distribution.shared.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "peak_distribution", "peak_distribution.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "transcription", "transcription.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "translation", "translation.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "splicesite", "splicesite.ui.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "motif", "motif.functions.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "motif", "motif.controls.R"), local = TRUE, encoding = "UTF-8")
source(app_path("modules", "motif", "motif.ui.R"), local = TRUE, encoding = "UTF-8")

if (!"app-static" %in% names(shiny::resourcePaths()) && file.exists(app_path("logo.webp"))) {
  addResourcePath("app-static", app_root)
}

browser_register_resource_paths()

react_shell_asset_tags <- function() {
  asset_dir <- app_path("www", "react", "app_shell", "assets")
  css_file <- file.path(asset_dir, "app-shell.css")
  js_file <- file.path(asset_dir, "app-shell.js")
  asset_tags <- list()
  asset_version <- function(path) {
    if (!file.exists(path)) {
      return(NULL)
    }

    unname(tools::md5sum(path))
  }

  if (file.exists(css_file)) {
    css_version <- asset_version(css_file)
    css_href <- "react/app_shell/assets/app-shell.css"

    if (!is.null(css_version)) {
      css_href <- sprintf("%s?v=%s", css_href, css_version)
    }

    asset_tags <- c(
      asset_tags,
      list(tags$link(rel = "stylesheet", href = css_href))
    )
  }

  if (file.exists(js_file)) {
    js_version <- asset_version(js_file)
    js_src <- "react/app_shell/assets/app-shell.js"

    if (!is.null(js_version)) {
      js_src <- sprintf("%s?v=%s", js_src, js_version)
    }

    asset_tags <- c(
      asset_tags,
      list(tags$script(type = "module", src = js_src))
    )
  }

  tagList(asset_tags)
}

app_head_resources <- function() {
  tags$head(
    tags$title("RNAmeta Online"),
    tags$meta(charset = "utf-8"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$link(rel = "stylesheet", href = "css/app.css?v=20260920a"),
    react_shell_asset_tags()
  )
}

ui <- shinyUI(
  tagList(
    navbarPage(
      title = "RNAmeta Online",
      id = "main_navbar",
      selected = "welcome",
      header = app_head_resources(),
      tabPanel("Welcome", value = "welcome", mod_welcome_ui("welcome")),
      tabPanel(
        "Upload / Run",
        value = "upload_run",
        mod_upload_run_ui("upload_run")
      ),
      tabPanel(
        "Peak Distribution",
        value = "peak_distribution",
        mod_peak_distribution_ui("peak_distribution")
      ),
      tabPanel(
        "Genome Browser",
        value = "genome_browser",
        mod_browser_ui("genome_browser")
      ),
      tabPanel(
        "Structure",
        value = "structure",
        mod_structure_ui("structure")
      ),
      navbarMenu(
        "Gene Statistics",
        tabPanel(
          "Gene Type",
          value = "gene_type",
          mod_gene_type_ui("gene_type")
        ),
        tabPanel(
          "Peak Gene Size",
          value = "peak_gene_size",
          mod_peak_gene_size_ui("peak_gene_size")
        ),
        tabPanel(
          "Gene Matrix",
          value = "gene_matrix",
          mod_gene_matrix_ui("gene_matrix")
        )
      ),
      navbarMenu(
        "Exon Statistics",
        tabPanel(
          "Peak Exon Size",
          value = "peak_exon_size",
          mod_peak_exon_size_ui("peak_exon_size")
        ),
        tabPanel(
          "Peak Exon Type",
          value = "peak_exon_type",
          mod_peak_exon_type_ui("peak_exon_type")
        ),
        tabPanel(
          "Peak Exon Num",
          value = "peak_exon_num",
          mod_peak_exon_num_ui("peak_exon_num")
        )
      ),
      navbarMenu(
        "Site",
        tabPanel(
          "Transcription",
          value = "transcription",
          mod_transcription_ui("transcription")
        ),
        tabPanel(
          "Translation",
          value = "translation",
          mod_translation_ui("translation")
        ),
        tabPanel(
          "Splicesite",
          value = "splicesite",
          mod_splicesite_ui("splicesite")
        )
      ),
      tabPanel(
        "Motif",
        value = "motif",
        mod_motif_ui("motif")
      ),
      tabPanel(
        "Help",
        value = "help",
        mod_help_ui("help")
      )
    ),
    tags$footer(
      class = "rm-footer",
      "RNAmeta Online | Rice Research Institute, Guangdong Academy of Agricultural Sciences | Guangzhou, Guangdong, China | 2026"
    )
  )
)
