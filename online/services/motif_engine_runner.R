MOTIF_ENGINE_BINARIES <- list(
  homer = "/home/shiny/miniconda3/envs/rnameta_motif/bin/findMotifs.pl",
  meme = "/home/shiny/miniconda3/envs/rnameta_motif/bin/meme",
  streme = "/home/shiny/miniconda3/envs/rnameta_motif/bin/streme"
)

MOTIF_ENGINE_THREADS <- list(
  homer = 32L,
  meme = 32L
)

motif_engine_result_directory <- function(engine, sample_name) {
  engine_key <- tolower(as.character(engine)[[1]])
  sample_key <- gsub("[^A-Za-z0-9._-]+", "-", as.character(sample_name)[[1]])
  file.path("motif_results", engine_key, sample_key)
}

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

motif_engine_binary_path <- function(engine) {
  engine_key <- tolower(as.character(engine)[[1]])
  binary_path <- MOTIF_ENGINE_BINARIES[[engine_key]]

  if (is.null(binary_path) || !nzchar(binary_path)) {
    stop(sprintf("Motif engine binary is not configured for '%s'.", engine_key))
  }

  binary_path
}

motif_engine_thread_count <- function(engine) {
  engine_key <- tolower(as.character(engine)[[1]])
  thread_count <- MOTIF_ENGINE_THREADS[[engine_key]]

  if (is.null(thread_count)) {
    return(NA_integer_)
  }

  motif_safe_integer(thread_count, NA_integer_)
}

motif_engine_run_script_name <- function(engine) {
  sprintf("run_%s.sh", tolower(as.character(engine)[[1]]))
}

motif_engine_expected_outputs <- function(engine) {
  engine_key <- tolower(as.character(engine)[[1]])

  if (identical(engine_key, "homer")) {
    return(c(
      "query.fa",
      "control.fa",
      "knownResults.txt",
      "homerResults.html",
      "homerMotifs.all.motifs"
    ))
  }

  if (identical(engine_key, "meme")) {
    return(c(
      "query.fa",
      "meme.txt",
      "meme.xml",
      "meme.html"
    ))
  }

  c(
    "query.fa",
    "streme.txt",
    "streme.xml",
    "streme.html"
  )
}

motif_engine_command_preview <- function(engine, sample_name, controls) {
  command <- c(motif_engine_binary_path(engine), motif_engine_command_args(engine, controls))
  paste(vapply(command, shQuote, character(1), type = "sh"), collapse = " ")
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

motif_engine_workspace_output_dir <- function(engine, package_dir) {
  file.path(package_dir, motif_engine_output_dir_name(engine))
}

motif_engine_env_vars <- function(binary_path) {
  binary_dir <- dirname(binary_path)
  current_path <- Sys.getenv("PATH", unset = "")
  updated_path <- if (nzchar(current_path)) {
    paste(binary_dir, current_path, sep = .Platform$path.sep)
  } else {
    binary_dir
  }

  c(
    sprintf("PATH=%s", updated_path)
  )
}

motif_run_engine <- function(engine, package_dir, controls, progress = NULL) {
  engine_key <- tolower(as.character(engine)[[1]])
  binary_path <- motif_engine_binary_path(engine_key)
  output_dir <- motif_engine_workspace_output_dir(engine_key, package_dir)
  stdout_log <- file.path(package_dir, sprintf("%s.stdout.log", engine_key))
  stderr_log <- file.path(package_dir, sprintf("%s.stderr.log", engine_key))

  if (!file.exists(binary_path)) {
    stop(sprintf("Configured motif engine binary does not exist: %s", binary_path))
  }

  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }

  if (!dir.exists(package_dir)) {
    dir.create(package_dir, recursive = TRUE, showWarnings = FALSE)
  }

  args <- motif_engine_command_args(engine_key, controls)
  timeout_seconds <- motif_engine_timeout_seconds(engine_key)
  motif_progress_update(
    progress,
    90,
    sprintf(
      "Running %s on the server workspace (timeout: %d min)",
      motif_engine_label(engine_key),
      timeout_seconds / 60L
    )
  )

  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  setwd(package_dir)

  unlink(stdout_log, force = TRUE)
  unlink(stderr_log, force = TRUE)

  exit_status <- suppressWarnings(
    system2(
      command = binary_path,
      args = args,
      env = motif_engine_env_vars(binary_path),
      stdout = stdout_log,
      stderr = stderr_log,
      timeout = timeout_seconds
    )
  )

  if (is.null(exit_status)) {
    exit_status <- 0L
  }

  if (!identical(as.integer(exit_status), 0L)) {
    if (identical(as.integer(exit_status), 124L)) {
      stop(sprintf(
        "%s exceeded the %d minute runtime limit.",
        motif_engine_label(engine_key),
        timeout_seconds / 60L
      ))
    }

    error_lines <- character()
    if (file.exists(stderr_log)) {
      error_lines <- readLines(stderr_log, warn = FALSE, encoding = "UTF-8")
    }
    if (!length(error_lines) && file.exists(stdout_log)) {
      error_lines <- readLines(stdout_log, warn = FALSE, encoding = "UTF-8")
    }
    error_lines <- error_lines[nzchar(trimws(error_lines))]
    error_message <- if (length(error_lines)) {
      tail(error_lines, 8L)
    } else {
      sprintf("%s exited with status %s.", motif_engine_label(engine_key), exit_status)
    }

    stop(paste(error_message, collapse = "\n"))
  }

  list(
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    stdout_log = normalizePath(stdout_log, winslash = "/", mustWork = FALSE),
    stderr_log = normalizePath(stderr_log, winslash = "/", mustWork = FALSE)
  )
}

motif_bundle_schema_notice <- function() {
  c(
    "Local development currently prepares motif jobs and can replay uploaded JSON result bundles.",
    "Recommended result bundle fields: engine, sampleName, sourceFile, motifs[].name, motifs[].consensus, motifs[].pwm."
  )
}
