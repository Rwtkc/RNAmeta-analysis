browser_supported_references <- function() {
  references <- list(
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
    hg38.gencode_v46 = list(
      id = "hg38.gencode_v46",
      name = "Homo sapiens (GRCh38, GENCODE v46)",
      fasta = "hg38.transcripts.fa",
      index = "hg38.transcripts.fa.fai",
      asset_dirs = c("hg38.gencode_v46", "igv_hg38.gencode_v46"),
      alias = "chrom_alias.tsv",
      chrom_sizes = "chrom.sizes",
      transcripts = "transcripts.bed12.gz",
      transcript_index = "transcripts.bed12.gz.tbi",
      annotation_views = list(
        transcript = list(
          label = "Transcript",
          track_name = "GENCODE v46 transcripts",
          files = c("transcripts.bed12.gz", "transcripts.bed12"),
          index = "transcripts.bed12.gz.tbi",
          color = "#65795b",
          height = 180
        ),
        five_utr = list(
          label = "5'UTR",
          track_name = "GENCODE v46 5'UTR",
          files = c("regions_5utr.bed.gz", "regions_5utr.bed"),
          index = "regions_5utr.bed.gz.tbi",
          color = "#80a876",
          height = 120
        ),
        cds = list(
          label = "CDS",
          track_name = "GENCODE v46 CDS",
          files = c("regions_cds.bed.gz", "regions_cds.bed"),
          index = "regions_cds.bed.gz.tbi",
          color = "#4d7848",
          height = 130
        ),
        three_utr = list(
          label = "3'UTR",
          track_name = "GENCODE v46 3'UTR",
          files = c("regions_3utr.bed.gz", "regions_3utr.bed"),
          index = "regions_3utr.bed.gz.tbi",
          color = "#9cb581",
          height = 120
        ),
        noncoding = list(
          label = "Non-coding",
          track_name = "GENCODE v46 non-coding transcripts",
          files = c("transcripts_noncoding.bed12.gz", "transcripts_noncoding.bed12"),
          index = "transcripts_noncoding.bed12.gz.tbi",
          color = "#80897a",
          height = 180
        )
      ),
      chromosomes = c(as.character(1:22), "X", "Y", "Mt"),
      default_locus = "1:1-100000"
    ),
    hg19.gencode_v19 = list(
      id = "hg19",
      name = "Homo sapiens (GRCh37, GENCODE v19)",
      fasta = "hg19.transcripts.fa",
      index = "hg19.transcripts.fa.fai",
      asset_dirs = c("hg19", "igv_hg19"),
      alias = "chrom_alias.tsv",
      chrom_sizes = "chrom.sizes",
      transcripts = "transcripts.bed12.gz",
      transcript_index = "transcripts.bed12.gz.tbi",
      annotation_views = list(
        transcript = list(
          label = "Transcript",
          track_name = "GENCODE v19 transcripts",
          files = c("transcripts.bed12.gz", "transcripts.bed12"),
          index = "transcripts.bed12.gz.tbi",
          color = "#65795b",
          height = 180
        ),
        five_utr = list(
          label = "5'UTR",
          track_name = "GENCODE v19 5'UTR",
          files = c("regions_5utr.bed.gz", "regions_5utr.bed"),
          index = "regions_5utr.bed.gz.tbi",
          color = "#80a876",
          height = 120
        ),
        cds = list(
          label = "CDS",
          track_name = "GENCODE v19 CDS",
          files = c("regions_cds.bed.gz", "regions_cds.bed"),
          index = "regions_cds.bed.gz.tbi",
          color = "#4d7848",
          height = 130
        ),
        three_utr = list(
          label = "3'UTR",
          track_name = "GENCODE v19 3'UTR",
          files = c("regions_3utr.bed.gz", "regions_3utr.bed"),
          index = "regions_3utr.bed.gz.tbi",
          color = "#9cb581",
          height = 120
        ),
        noncoding = list(
          label = "Non-coding",
          track_name = "GENCODE v19 non-coding transcripts",
          files = c("transcripts_noncoding.bed12.gz", "transcripts_noncoding.bed12"),
          index = "transcripts_noncoding.bed12.gz.tbi",
          color = "#80897a",
          height = 180
        )
      ),
      chromosomes = c(as.character(1:22), "X", "Y", "Mt"),
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
  references$hg38 <- utils::modifyList(references$hg38.gencode_v46, list(
    id = "hg38",
    name = "Homo sapiens (hg38, GENCODE v46)"
  ))
  references$hg19 <- references$hg19.gencode_v19
  references
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
