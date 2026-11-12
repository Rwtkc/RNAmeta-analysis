peak_gene_size_compute_lengths <- function(peaks, transcript_ranges, txlens, sample_name) {
  transcript_ranges_unstranded <- transcript_ranges
  GenomicRanges::strand(transcript_ranges_unstranded) <- "*"
  overlaps <- GenomicRanges::findOverlaps(peaks, transcript_ranges_unstranded)

  if (!length(overlaps)) {
    return(data.table::data.table(
      tx_name = character(),
      tx_len = numeric(),
      Sample = character()
    ))
  }

  hit_index <- which(!duplicated(S4Vectors::queryHits(overlaps)))
  subject_index <- S4Vectors::subjectHits(overlaps[hit_index])
  transcript_hits <- unique(as.character(transcript_ranges$tx_name[subject_index]))

  txlens_dt <- data.table::as.data.table(txlens)
  lengths_dt <- txlens_dt[tx_name %in% transcript_hits, .(tx_name, tx_len)]
  lengths_dt <- unique(lengths_dt, by = "tx_name")
  lengths_dt <- lengths_dt[!is.na(tx_len) & tx_len > 0]
  lengths_dt[, Sample := sample_name]
  lengths_dt[order(Sample, tx_len)]
}

peak_gene_size_detect_delimiter <- function(filepath) {
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

peak_gene_size_load_rda_object <- function(path) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  env[[object_names[[1]]]]
}

peak_gene_size_resolve_annotation_bundle <- function(species) {
  bundle <- peak_distribution_resolve_annotation_bundle(species)
  list(
    species = bundle$species,
    species_abb = bundle$species_abb,
    txdb_path = bundle$txdb_path,
    txlens_path = bundle$txlens_path
  )
}

peak_gene_size_load_annotation_bundle <- function(species) {
  bundle <- peak_gene_size_resolve_annotation_bundle(species)
  txdb <- AnnotationDbi::loadDb(bundle$txdb_path)
  GenomeInfoDb::seqlevels(txdb) <- sub("^([^chr])", "chr\\1", GenomeInfoDb::seqlevels(txdb))
  txlens <- peak_gene_size_load_rda_object(bundle$txlens_path)

  list(
    species = species,
    species_abb = bundle$species_abb,
    txdb = txdb,
    txlens = txlens,
    txdb_path = bundle$txdb_path,
    txlens_path = bundle$txlens_path
  )
}

peak_gene_size_valid_chromosomes <- function(txdb) {
  chroms <- GenomeInfoDb::seqlevels(txdb)
  unique(ifelse(startsWith(chroms, "chr"), chroms, paste0("chr", chroms)))
}

peak_gene_size_prepare_transcript_ranges <- function(txdb) {
  transcripts <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"))
  GenomeInfoDb::seqlevels(transcripts) <- sub("^([^chr])", "chr\\1", GenomeInfoDb::seqlevels(transcripts))
  transcripts
}

peak_gene_size_read_peaks <- function(source_file, txdb, sample_name = "Group1") {
  delimiter <- peak_gene_size_detect_delimiter(source_file)
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

  valid_chromosomes <- peak_gene_size_valid_chromosomes(txdb)
  peak_df[, chr := trimws(as.character(chr))]
  peak_df[, chr := ifelse(startsWith(chr, "chr"), chr, paste0("chr", chr))]
  peak_df <- peak_df[chr %in% valid_chromosomes]
  peak_df[!strand %in% c("+", "-", "*"), strand := "*"]
  peak_df[, score := as.numeric(score)]
  peak_df[, start := start + 1]
  peak_df <- peak_df[chr %in% GenomeInfoDb::seqlevels(txdb)]

  peak_gr <- GenomicRanges::GRanges(
    seqnames = peak_df$chr,
    ranges = IRanges::IRanges(peak_df$start, peak_df$end)
  )

  if (ncol(peak_df) >= 6) {
    GenomicRanges::strand(peak_gr) <- peak_df$strand
  }

  peak_gr$groupname <- sample_name
  peak_gr$peakname <- paste0("Locus_", seq_len(nrow(peak_df)))
  unique(peak_gr)
}

peak_gene_size_default_controls <- function() {
  list()
}
