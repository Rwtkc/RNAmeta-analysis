structure_extract_sequence <- function(upload_context, selection, window_size = "fullsize") {
  spec <- structure_transcript_reference_spec(upload_context$species)
  if (is.null(spec) || !structure_transcript_reference_ready(spec)) {
    stop("The selected species does not have ready cDNA FASTA and GFF3 references for Structure.")
  }

  interval <- selection$interval
  transcript <- structure_map_bed_to_transcript(selection, spec)
  fasta_path <- app_path("data", "fa", spec$cdna)
  seqname <- structure_match_transcript_seqname(
    transcript_id = transcript$transcript_id,
    transcript_name = transcript$transcript_name,
    fasta_path = fasta_path
  )
  fasta_index <- Rsamtools::scanFaIndex(fasta_path)
  seq_widths <- IRanges::width(fasta_index)
  seq_names <- as.character(GenomeInfoDb::seqnames(fasta_index))
  seq_length <- as.integer(seq_widths[match(seqname, seq_names)])

  range <- GenomicRanges::GRanges(seqnames = seqname, ranges = IRanges::IRanges(start = 1L, end = seq_length))
  dna <- Rsamtools::scanFa(Rsamtools::FaFile(fasta_path), param = range)[[1]]
  if (length(dna) == 0) {
    stop("No transcript sequence was extracted for the selected BED interval.")
  }

  sequence <- toupper(as.character(dna))
  sequence <- chartr("T", "U", sequence)

  list(
    sequence = sequence,
    locus = sprintf("%s:%s-%s", transcript$transcript_id, 1L, nchar(sequence)),
    bed_locus = sprintf("%s:%s-%s", interval$chrom, transcript$bed_start_1based, transcript$bed_end_1based),
    strand = transcript$strand %||% ".",
    window_size = "fullsize",
    interval_length = as.integer(interval$end - interval$start),
    used_bed_interval = FALSE,
    matched_bed_row = isTRUE(selection$matched_bed_row),
    transcript = transcript,
    highlight_ranges = transcript$highlight_ranges,
    highlight_start = transcript$highlight_start,
    highlight_end = transcript$highlight_end,
    transcript_length = transcript$transcript_length,
    transcript_seqname = seqname
  )
}

structure_rnafold_path <- function() {
  candidates <- c(
    Sys.getenv("RNAMETA_RNAFOLD", ""),
    "/home/shiny/miniconda3/envs/RNAhalf_RNAfold/bin/RNAfold",
    Sys.which("RNAfold")
  )
  candidates <- candidates[nzchar(candidates)]
  matches <- candidates[file.exists(path.expand(candidates))]

  if (length(matches)) {
    return(path.expand(matches[[1]]))
  }

  ""
}

structure_fold_with_rnafold <- function(sequence) {
  rnafold <- structure_rnafold_path()
  if (!nzchar(rnafold)) {
    return(NULL)
  }

  input <- tempfile(pattern = "rnameta_structure_", fileext = ".fa")
  writeLines(c(">rnameta", sequence), input)
  on.exit(unlink(input, force = TRUE), add = TRUE)

  output <- tryCatch(
    system2(rnafold, args = "--noPS", stdin = input, stdout = TRUE, stderr = TRUE),
    error = function(...) character()
  )
  structure_line <- output[grepl("^[().]+\\s", output)]
  if (!length(structure_line)) {
    details <- paste(utils::tail(output, 6L), collapse = " | ")
    folded <- NA_character_
    attr(folded, "error") <- if (nzchar(details)) details else "RNAfold produced no parseable output."
    return(folded)
  }

  sub("\\s.*$", "", structure_line[[1]])
}

structure_can_pair <- function(left, right) {
  paste0(left, right) %in% c("AU", "UA", "GC", "CG", "GU", "UG")
}

structure_fold_preview <- function(sequence) {
  bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
  n <- length(bases)
  structure <- rep(".", n)
  left <- 1L
  right <- n

  while (left + 3L < right) {
    if (structure_can_pair(bases[[left]], bases[[right]])) {
      structure[[left]] <- "("
      structure[[right]] <- ")"
      left <- left + 1L
      right <- right - 1L
    } else {
      right_match <- which(vapply(seq.int(right, left + 4L), function(candidate) {
        structure_can_pair(bases[[left]], bases[[candidate]])
      }, logical(1)))
      if (length(right_match)) {
        candidate <- seq.int(right, left + 4L)[[right_match[[1]]]]
        structure[[left]] <- "("
        structure[[candidate]] <- ")"
        left <- left + 1L
        right <- candidate - 1L
      } else {
        left <- left + 1L
      }
    }
  }

  paste(structure, collapse = "")
}

structure_fold_sequence <- function(sequence) {
  folded <- structure_fold_with_rnafold(sequence)
  if (!is.null(folded) && nchar(folded) == nchar(sequence)) {
    return(list(structure = folded, method = "RNAfold"))
  }

  rnafold_error <- attr(folded, "error") %||% ""
  stop(sprintf(
    "RNAfold did not return a valid dot-bracket structure. %s",
    rnafold_error
  ))
}

structure_progress_start <- function(session) {
  progress <- shiny::Progress$new(session = session, min = 0, max = 100)
  progress$set(
    value = 0,
    message = "Running Structure",
    detail = "0% | Initializing structure workflow"
  )
  progress
}

structure_progress_update <- function(progress, value, detail) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  safe_value <- as.integer(round(max(0, min(100, value))))
  progress$set(
    value = safe_value,
    message = "Running Structure",
    detail = sprintf("%s%% | %s", safe_value, detail)
  )

  invisible(NULL)
}

structure_progress_close <- function(progress) {
  if (is.null(progress)) {
    return(invisible(NULL))
  }

  Sys.sleep(0.6)
  progress$close()
  invisible(NULL)
}

STRUCTURE_JOB_ACTIVE_LIMIT <- 3L
STRUCTURE_JOB_ACTIVE_STATUSES <- c("queued", "running")
