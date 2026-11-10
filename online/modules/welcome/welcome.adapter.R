welcome_shell_config <- function() {
  list(
    hero = list(
      eyebrow = "mRNA Modification Analysis Workspace",
      title = "RNAmeta Online",
      description = paste(
        "RNAmeta Online is an advanced workspace for meta-analysis and functional",
        "annotation of mRNA modification signals. It accepts BED-formatted",
        "chromosomal coordinates from widely used peak callers and reorganizes",
        "downstream interpretation into clearer module-oriented views."
      ),
      logo_src = "app-static/logo.webp",
      supporting = paste(
        "The platform supports BED coordinates compatible with tools such as",
        "MACS, exomePeak, and MetaTeak, and helps users move from genomic",
        "positions to transcript-aware interpretation more efficiently."
      ),
      badges = c(
        "BED-compatible input",
        "Meta-analysis and annotation",
        "Functional interpretation"
      )
    ),
    modules = list(
      list(
        id = "upload_run",
        title = "Upload / Run",
        hint = "Input and session setup"
      ),
      list(
        id = "peak_distribution",
        title = "Peak Distribution",
        hint = "Feature-level distribution"
      ),
      list(
        id = "motif",
        title = "Motif",
        hint = "Sequence context"
      ),
      list(
        id = "feature_statistics",
        title = "Feature Statistics",
        hint = "Gene and exon summaries"
      )
    )
  )
}

welcome_shell_json <- function(config = welcome_shell_config()) {
  jsonlite::toJSON(
    config,
    auto_unbox = TRUE,
    pretty = FALSE,
    null = "null"
  )
}
