peak_distribution_detect_delimiter <- function(filepath) {
  lines <- readLines(filepath, n = 5, warn = FALSE)
  delimiters <- c(",", "\t", " ")
  counts <- vapply(
    delimiters,
    function(delimiter) {
      sum(vapply(lines, function(line) length(strsplit(line, delimiter, fixed = TRUE)[[1]]), integer(1)))
    },
    integer(1)
  )
  delimiters[[which.max(counts)]]
}

peak_distribution_load_rda_object <- function(path) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  env[[object_names[[1]]]]
}

peak_distribution_annotation_data_root <- function() {
  app_path("data")
}

peak_distribution_species_code <- function(species) {
  species_value <- as.character(species)[[1]]
  code_match <- regexec("\\(([^)]+)\\)", species_value)
  code_parts <- regmatches(species_value, code_match)[[1]]

  if (length(code_parts) >= 2L && nzchar(code_parts[[2]])) {
    return(code_parts[[2]])
  }

  gsub("[^A-Za-z0-9._-]+", "_", trimws(species_value))
}

peak_distribution_species_tokens <- function(species) {
  species_value <- as.character(species)[[1]]
  species_name <- trimws(sub("\\s*\\(.*\\)$", "", species_value))
  tokens <- unique(strsplit(tolower(species_name), "[^a-z0-9]+")[[1]])
  tokens[nzchar(tokens) & nchar(tokens) >= 4L]
}

peak_distribution_find_annotation_file <- function(data_dir, species, extension, species_code = NULL) {
  if (!dir.exists(data_dir)) {
    stop(sprintf("Required annotation directory not found: %s", data_dir))
  }

  if (is.null(species_code) || !nzchar(species_code)) {
    species_code <- peak_distribution_species_code(species)
  }

  exact_path <- file.path(data_dir, paste0(species_code, extension))
  if (file.exists(exact_path)) {
    return(exact_path)
  }

  files <- list.files(data_dir, pattern = paste0("\\Q", extension, "\\E$"), full.names = TRUE)
  basenames <- tools::file_path_sans_ext(basename(files))
  lower_basenames <- tolower(basenames)

  code_hits <- files[lower_basenames == tolower(species_code)]
  if (length(code_hits) == 1L) {
    return(code_hits[[1]])
  }

  tokens <- peak_distribution_species_tokens(species)
  token_hits <- files[vapply(lower_basenames, function(name) {
    any(vapply(tokens, function(token) grepl(token, name, fixed = TRUE), logical(1)))
  }, logical(1))]

  if (length(token_hits) == 1L) {
    return(token_hits[[1]])
  }

  stop(sprintf(
    paste(
      "Annotation file for species '%s' was not found in %s.",
      "Expected either %s or a uniquely matching file in that directory."
    ),
    as.character(species)[[1]],
    data_dir,
    basename(exact_path)
  ))
}

peak_distribution_resolve_annotation_bundle <- function(species) {
  data_root <- peak_distribution_annotation_data_root()
  species_abb <- peak_distribution_species_code(species)
  list(
    species = species,
    species_abb = species_abb,
    txdb_path = peak_distribution_find_annotation_file(
      file.path(data_root, "txdb"),
      species = species,
      extension = ".txdb.sqlite",
      species_code = species_abb
    ),
    txlens_path = peak_distribution_find_annotation_file(
      file.path(data_root, "txlens"),
      species = species,
      extension = ".txlens.rda",
      species_code = species_abb
    ),
    gff_path = peak_distribution_find_annotation_file(
      file.path(data_root, "gff"),
      species = species,
      extension = ".gff.rda",
      species_code = species_abb
    )
  )
}

peak_distribution_load_annotation_bundle <- function(species) {
  bundle <- peak_distribution_resolve_annotation_bundle(species)
  txdb <- AnnotationDbi::loadDb(bundle$txdb_path)
  GenomeInfoDb::seqlevels(txdb) <- ifelse(startsWith(GenomeInfoDb::seqlevels(txdb), "chr"), GenomeInfoDb::seqlevels(txdb), paste0("chr", GenomeInfoDb::seqlevels(txdb)))

  txlens <- peak_distribution_load_rda_object(bundle$txlens_path)
  gff <- peak_distribution_load_rda_object(bundle$gff_path)

  if ("gene_biotype" %in% names(S4Vectors::mcols(gff))) {
    names(S4Vectors::mcols(gff))[names(S4Vectors::mcols(gff)) == "gene_biotype"] <- "gene_type"
  }

  if ("transcript_biotype" %in% names(S4Vectors::mcols(gff))) {
    names(S4Vectors::mcols(gff))[names(S4Vectors::mcols(gff)) == "transcript_biotype"] <- "transcript_type"
  }

  GenomeInfoDb::seqlevels(gff) <- ifelse(startsWith(GenomeInfoDb::seqlevels(gff), "chr"), GenomeInfoDb::seqlevels(gff), paste0("chr", GenomeInfoDb::seqlevels(gff)))

  list(
    species = species,
    species_abb = bundle$species_abb,
    txdb = txdb,
    txlens = txlens,
    gff = gff,
    txdb_path = bundle$txdb_path,
    txlens_path = bundle$txlens_path,
    gff_path = bundle$gff_path
  )
}

peak_distribution_valid_chromosomes <- function(txdb) {
  chroms <- GenomeInfoDb::seqlevels(txdb)
  unique(ifelse(startsWith(chroms, "chr"), chroms, paste0("chr", chroms)))
}

peak_distribution_read_peaks <- function(source_file, txdb, sample_name = "Group1", subset_for_motif = "All") {
  delimiter <- peak_distribution_detect_delimiter(source_file)
  peak_df <- data.table::fread(source_file, header = FALSE, sep = delimiter)

  if (ncol(peak_df) == 3) {
    peak_df[, site := paste0("site", .I)]
    peak_df[, dot_column := "."]
    peak_df[, star_column := "*"]
  }

  if (ncol(peak_df) == 4) {
    peak_df[, dot_column := "."]
    peak_df[, star_column := "*"]
  }

  if (ncol(peak_df) == 5) {
    peak_df[, star_column := "*"]
  }

  peak_df <- peak_df[, 1:6]

  if (sum(!is.na(as.numeric(peak_df[[5]]))) > sum(vapply(peak_df[, -5, with = FALSE], function(col) sum(is.na(as.numeric(col))), integer(1)))) {
    peak_df <- peak_df[!is.na(as.numeric(V5)), ]
  }

  peak_df <- peak_df[peak_df[[3]] > peak_df[[2]], ]

  bed_headers <- c("chr", "start", "end", "name", "score", "strand")
  data.table::setnames(peak_df, bed_headers[seq_len(ncol(peak_df))])

  valid_chromosomes <- peak_distribution_valid_chromosomes(txdb)
  peak_df[, chr := trimws(as.character(chr))]
  peak_df[, chr := ifelse(startsWith(chr, "chr"), chr, paste0("chr", chr))]
  peak_df <- peak_df[chr %in% valid_chromosomes]

  peak_df[!strand %in% c("+", "-", "*"), strand := "*"]
  peak_df[, score := as.numeric(score)]

  if (ncol(peak_df) >= 5 && !identical(subset_for_motif, "All")) {
    data.table::setorder(peak_df, -score)
  }

  peak_df[, start := start + 1]
  GenomeInfoDb::seqlevels(txdb) <- ifelse(startsWith(GenomeInfoDb::seqlevels(txdb), "chr"), GenomeInfoDb::seqlevels(txdb), paste0("chr", GenomeInfoDb::seqlevels(txdb)))
  peak_df <- peak_df[chr %in% GenomeInfoDb::seqlevels(txdb)]

  peak_gr <- GenomicRanges::GRanges(
    seqnames = peak_df$chr,
    ranges = IRanges::IRanges(peak_df$start, peak_df$end)
  )

  if (ncol(peak_df) >= 6) {
    GenomicRanges::strand(peak_gr) <- peak_df$strand
  }

  if (ncol(peak_df) > 3) {
    for (column_name in colnames(peak_df)[4:ncol(peak_df)]) {
      if (!identical(column_name, "strand")) {
        S4Vectors::mcols(peak_gr)[[column_name]] <- peak_df[[column_name]]
      }
    }
  }

  peak_gr$groupname <- sample_name
  peak_gr$peakname <- paste0("Locus_", seq_len(nrow(peak_df)))
  unique(peak_gr)
}
