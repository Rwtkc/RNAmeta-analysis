browser_supported_references <- function() {
  list(
    ara_TAIR10 = list(
      id = "ara_TAIR10",
      name = "Arabidopsis thaliana (TAIR10)",
      fasta = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
      index = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai",
      asset_dirs = c("ara_TAIR10", "igv_ara_TAIR10"),
      alias = "chrom_alias.tsv",
      chrom_sizes = "chrom.sizes",
      transcripts = "transcripts.bed12.gz",
      transcript_index = "transcripts.bed12.gz.tbi",
      annotation_views = list(
        transcript = list(
          label = "Transcript",
          track_name = "TAIR10 transcripts",
          files = c("transcripts.bed12.gz", "transcripts.bed12"),
          index = "transcripts.bed12.gz.tbi",
          color = "#65795b",
          height = 180
        ),
        five_utr = list(
          label = "5'UTR",
          track_name = "TAIR10 5'UTR",
          files = c("regions_5utr.bed.gz", "regions_5utr.bed"),
          index = "regions_5utr.bed.gz.tbi",
          color = "#80a876",
          height = 120
        ),
        cds = list(
          label = "CDS",
          track_name = "TAIR10 CDS",
          files = c("regions_cds.bed.gz", "regions_cds.bed"),
          index = "regions_cds.bed.gz.tbi",
          color = "#4d7848",
          height = 130
        ),
        three_utr = list(
          label = "3'UTR",
          track_name = "TAIR10 3'UTR",
          files = c("regions_3utr.bed.gz", "regions_3utr.bed"),
          index = "regions_3utr.bed.gz.tbi",
          color = "#9cb581",
          height = 120
        ),
        noncoding = list(
          label = "Non-coding",
          track_name = "TAIR10 non-coding transcripts",
          files = c("transcripts_noncoding.bed12.gz", "transcripts_noncoding.bed12"),
          index = "transcripts_noncoding.bed12.gz.tbi",
          color = "#80897a",
          height = 180
        )
      ),
      chromosomes = c("1", "2", "3", "4", "5", "Mt", "Pt"),
      default_locus = "1:1-100000"
    ),
    osa_IRGSP_1 = list(
      id = "osa_IRGSP_1",
      name = "Oryza sativa (IRGSP 1.0)",
      fasta = "Oryza_sativa.IRGSP-1.0.dna.toplevel.fa",
      index = "Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.fai",
      asset_dirs = c("osa_IRGSP_1", "igv_osa_IRGSP_1"),
      alias = "chrom_alias.tsv",
      chrom_sizes = "chrom.sizes",
      transcripts = "transcripts.bed12.gz",
      transcript_index = "transcripts.bed12.gz.tbi",
      annotation_views = list(
        transcript = list(
          label = "Transcript",
          track_name = "IRGSP 1.0 transcripts",
          files = c("transcripts.bed12.gz", "transcripts.bed12"),
          index = "transcripts.bed12.gz.tbi",
          color = "#65795b",
          height = 180
        ),
        five_utr = list(
          label = "5'UTR",
          track_name = "IRGSP 1.0 5'UTR",
          files = c("regions_5utr.bed.gz", "regions_5utr.bed"),
          index = "regions_5utr.bed.gz.tbi",
          color = "#80a876",
          height = 120
        ),
        cds = list(
          label = "CDS",
          track_name = "IRGSP 1.0 CDS",
          files = c("regions_cds.bed.gz", "regions_cds.bed"),
          index = "regions_cds.bed.gz.tbi",
          color = "#4d7848",
          height = 130
        ),
        three_utr = list(
          label = "3'UTR",
          track_name = "IRGSP 1.0 3'UTR",
          files = c("regions_3utr.bed.gz", "regions_3utr.bed"),
          index = "regions_3utr.bed.gz.tbi",
          color = "#9cb581",
          height = 120
        ),
        noncoding = list(
          label = "Non-coding",
          track_name = "IRGSP 1.0 non-coding transcripts",
          files = c("transcripts_noncoding.bed12.gz", "transcripts_noncoding.bed12"),
          index = "transcripts_noncoding.bed12.gz.tbi",
          color = "#80897a",
          height = 180
        )
      ),
      chromosomes = c(as.character(1:12), "Mt", "Pt"),
      default_locus = "1:1-100000"
    )
  )
}

browser_reference_spec <- function(species) {
  references <- browser_supported_references()
  species <- as.character(species %||% "")

  if (!species %in% names(references)) {
    return(NULL)
  }

  references[[species]]
}

browser_reference_files_ready <- function(spec) {
  if (is.null(spec)) {
    return(FALSE)
  }

  file.exists(app_path("data", "fa", spec$fasta)) &&
    file.exists(app_path("data", "fa", spec$index))
}

browser_igv_asset_dir <- function(spec) {
  candidates <- file.path(app_path("data", "igv"), spec$asset_dirs)
  matches <- candidates[dir.exists(candidates)]

  if (!length(matches)) {
    return(NULL)
  }

  matches[[1]]
}

browser_igv_asset_file <- function(spec, filename) {
  asset_dir <- browser_igv_asset_dir(spec)

  if (is.null(asset_dir)) {
    return(NULL)
  }

  path <- file.path(asset_dir, filename)
  if (!file.exists(path)) {
    return(NULL)
  }

  path
}

browser_igv_asset_first <- function(spec, filenames) {
  filenames <- as.character(filenames %||% character())
  for (filename in filenames) {
    path <- browser_igv_asset_file(spec, filename)
    if (!is.null(path)) {
      return(path)
    }
  }

  NULL
}

browser_igv_asset_url <- function(path) {
  root <- normalizePath(app_path("data", "igv"), winslash = "/", mustWork = TRUE)
  relative <- sub(
    paste0("^", gsub("\\\\", "/", root), "/?"),
    "",
    normalizePath(path, winslash = "/", mustWork = TRUE)
  )

  paste0("rnameta-igv-assets/", relative)
}

browser_sanitize_filename <- function(value) {
  value <- basename(as.character(value %||% "track.bed"))
  gsub("[^A-Za-z0-9_.-]", "_", value)
}

browser_public_upload_path <- function(file_entry, cache_scope = NULL) {
  source_path <- file_entry$file_path
  if (is.null(source_path) || !nzchar(source_path) || !file.exists(source_path)) {
    return(NULL)
  }

  upload_dir <- browser_public_subdir("uploads", cache_scope)

  fingerprint <- paste(
    tools::md5sum(source_path),
    file.info(source_path)$size,
    sep = "_"
  )
  target_name <- sprintf(
    "%s_%s",
    substr(gsub("[^A-Za-z0-9]", "", fingerprint), 1, 16),
    browser_sanitize_filename(file_entry$file_name)
  )
  target_path <- file.path(upload_dir, target_name)

  if (!file.exists(target_path)) {
    file.copy(source_path, target_path, overwrite = TRUE)
  }

  target_path
}
