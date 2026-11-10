meta_plot_load_rda_object <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required annotation file not found: %s", path))
  }

  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)

  if (length(object_names) == 0) {
    stop(sprintf("No object found in R data file: %s", path))
  }

  env[[object_names[[1]]]]
}

resolve_meta_plot_annotation_bundle <- function(species) {
  peak_distribution_resolve_annotation_bundle(species)
}

load_meta_plot_annotation_bundle <- function(species) {
  meta_plot_cached_annotation_bundle(
    species = species,
    loader = function() {
      bundle <- resolve_meta_plot_annotation_bundle(species)
      txdb <- AnnotationDbi::loadDb(bundle$txdb_path)
      GenomeInfoDb::seqlevels(txdb) <- ifelse(startsWith(GenomeInfoDb::seqlevels(txdb), "chr"), GenomeInfoDb::seqlevels(txdb), paste0("chr", GenomeInfoDb::seqlevels(txdb)))

      list(
        species = species,
        txdb = txdb,
        txlens = meta_plot_load_rda_object(bundle$txlens_path),
        gff = meta_plot_load_rda_object(bundle$gff_path),
        txdb_path = bundle$txdb_path,
        txlens_path = bundle$txlens_path,
        gff_path = bundle$gff_path
      )
    }
  )
}

meta_plot_detect_delimiter <- function(file_path) {
  first_line <- readLines(file_path, n = 1, warn = FALSE)

  if (length(first_line) == 0 || !nzchar(first_line)) {
    return("\t")
  }

  if (grepl("\t", first_line, fixed = TRUE)) {
    return("\t")
  }

  if (grepl(",", first_line, fixed = TRUE)) {
    return(",")
  }

  if (grepl(";", first_line, fixed = TRUE)) {
    return(";")
  }

  ""
}

meta_plot_prepare_bed_file <- function(file_path, txdb) {
  if (is.null(file_path) || !nzchar(file_path) || !file.exists(file_path)) {
    stop("Saved BED file was not found. Please save the upload again.")
  }

  delimiter <- meta_plot_detect_delimiter(file_path)
  data <- data.table::fread(file_path, sep = delimiter, header = FALSE, data.table = TRUE)

  if (ncol(data) < 3) {
    stop("BED input must contain at least 3 columns.")
  }

  if (ncol(data) == 3) {
    data[, V4 := paste0("site", .I)]
    data[, V5 := "."]
    data[, V6 := "*"]
  } else if (ncol(data) == 4) {
    data[, V5 := "."]
    data[, V6 := "*"]
  } else if (ncol(data) == 5) {
    data[, V6 := "*"]
  }

  data <- data[, 1:6]
  data[, V1 := trimws(as.character(V1))]
  data[, V1 := ifelse(startsWith(V1, "chr"), V1, paste0("chr", V1))]
  data[, V2 := suppressWarnings(as.numeric(V2))]
  data[, V3 := suppressWarnings(as.numeric(V3))]
  data <- data[!is.na(V2) & !is.na(V3)]
  data <- data[V3 > V2]

  valid_chromosomes <- GenomeInfoDb::seqlevels(txdb)
  data <- data[V1 %in% valid_chromosomes]

  if (nrow(data) == 0) {
    stop("No valid BED intervals remained after transcriptome chromosome filtering.")
  }

  output_file <- tempfile(pattern = "meta_plot_", fileext = ".bed")
  data.table::fwrite(data, file = output_file, sep = "\t", col.names = FALSE, scipen = 999L)
  output_file
}

meta_plot_group_name <- function(file_info) {
  if (!is.null(file_info$file_name) && nzchar(file_info$file_name)) {
    return(tools::file_path_sans_ext(basename(file_info$file_name)))
  }

  "Group1"
}

meta_plot_build_series_metadata <- function(saved_files) {
  original_names <- vapply(saved_files, meta_plot_group_name, character(1))
  display_names <- sprintf("Sample %d", seq_along(original_names))

  lapply(seq_along(saved_files), function(index) {
    list(
      display_name = display_names[[index]],
      original_name = original_names[[index]]
    )
  })
}

meta_plot_require_upload <- function(upload_context) {
  if (is.null(upload_context)) {
    stop("Save a BED file in Upload / Run before running Meta Plot.")
  }

  if (is.null(upload_context$species) || !nzchar(upload_context$species)) {
    stop("A species must be saved in Upload / Run before running Meta Plot.")
  }

  saved_files <- upload_run_context_saved_files(upload_context)
  if (!length(saved_files)) {
    stop("At least one BED file must be saved in Upload / Run before running Meta Plot.")
  }

  upload_context
}
