MOTIF_ENGINE_THREADS <- list(
  homer = 32L,
  meme = 32L
)

motif_engine_output_dir_name <- function(engine) {
  engine_key <- tolower(as.character(engine)[[1]])

  if (identical(engine_key, "homer")) {
    return("homer_out")
  }

  if (identical(engine_key, "meme")) {
    return("meme_out")
  }

  "streme_out"
}

motif_engine_thread_count <- function(engine) {
  engine_key <- tolower(as.character(engine)[[1]])
  thread_count <- MOTIF_ENGINE_THREADS[[engine_key]]

  if (is.null(thread_count)) {
    return(NA_integer_)
  }

  motif_safe_integer(thread_count, NA_integer_)
}

motif_engine_command_args <- function(engine, controls) {
  engine_key <- tolower(as.character(engine)[[1]])
  motif_width <- motif_safe_integer(controls$motif_width, 12L)
  max_motifs <- motif_safe_integer(controls$max_motifs, 6L)
  output_dir <- motif_engine_output_dir_name(engine_key)
  thread_count <- motif_engine_thread_count(engine_key)
  matched <- identical(controls$analysis_mode, "matched_full")

  if (matched) {
    if (identical(engine_key, "homer")) {
      return(c("query.fa", "fasta", output_dir, "-fastaBg", "control.fa",
        "-nogo", "-noknown", "-rna", "-nlen", "1", "-S", as.character(max_motifs),
        "-p", as.character(thread_count), "-len", as.character(motif_width)))
    }
    if (identical(engine_key, "meme")) {
      return(c("query.fa", "-neg", "control.fa", "-objfun", "de", "-mod", "zoops",
        "-dna", "-bfile", "shared_background.txt", "-nmotifs", as.character(max_motifs),
        "-minw", as.character(motif_width), "-maxw", as.character(motif_width),
        "-seed", "0", "-p", as.character(thread_count), "-oc", output_dir))
    }
    if (identical(engine_key, "streme")) {
      return(c("--p", "query.fa", "--n", "control.fa", "--rna", "--order", "1",
        "--nmotifs", as.character(max_motifs), "--minw", as.character(motif_width),
        "--maxw", as.character(motif_width), "--seed", "0", "--oc", output_dir))
    }
    stop("Unsupported motif engine for matched full-peak analysis.")
  }

  if (identical(engine_key, "homer")) {
    return(c(
      "query.fa",
      "fasta",
      output_dir,
      "-fasta",
      "control.fa",
      "-nogo",
      "-noknown",
      "-rna",
      "-p",
      as.character(thread_count),
      "-len",
      as.character(motif_width)
    ))
  }

  if (identical(engine_key, "meme")) {
    return(c(
      "query.fa",
      "-mod",
      "zoops",
      "-dna",
      "-nmotifs",
      as.character(max_motifs),
      "-minw",
      as.character(motif_width),
      "-maxw",
      as.character(motif_width),
      "-revcomp",
      "-p",
      as.character(thread_count),
      "-oc",
      output_dir
    ))
  }

  c(
    "--p",
    "query.fa",
    "--dna",
    "--nmotifs",
    as.character(max_motifs),
    "--minw",
    as.character(motif_width),
    "--maxw",
    as.character(motif_width),
    "--oc",
    output_dir
  )
}
