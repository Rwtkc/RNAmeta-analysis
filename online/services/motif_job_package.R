motif_build_job_manifest <- function(result, package_files, command_preview, script_name) {
  list(
    engine = result$engine,
    analysisMode = result$analysis_mode,
    engineLabel = result$engine_label,
    jobId = if (is.null(result$job_id)) "" else result$job_id,
    sampleName = result$sample_name,
    sourceFile = result$source_file,
    sourceInputMd5 = result$source_input_md5,
    species = result$species,
    subset = result$sites_subset,
    motifWidth = result$motif_width,
    maxMotifs = result$max_motifs,
    peakCount = result$peak_count,
    commandPreview = command_preview,
    runScript = script_name,
    packageFiles = package_files,
    expectedOutputs = result$expected_outputs,
    createdAt = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  )
}

motif_prepare_job_package <- function(result, motif_regions, annotation_bundle = NULL, exclusion = motif_regions, progress = NULL) {
  package_dir <- motif_prepare_job_directory(result$engine, result$sample_name, result$job_id)
  matched <- identical(result$analysis_mode, "matched_full")
  matching_summary <- NULL
  extra_files <- character()
  shared_background <- NULL
  if (matched) {
    if (is.null(annotation_bundle)) annotation_bundle <- peak_distribution_load_annotation_bundle(result$species)
    prepared <- motif_prepare_matched_input(motif_regions, result, annotation_bundle, exclusion, package_dir, progress)
    query_regions <- prepared$query_regions
    query_sequences <- prepared$query_sequences
    control_regions <- prepared$control_regions
    control_sequences <- prepared$control_sequences
    matching_summary <- prepared$matching_summary
    extra_files <- prepared$extra_files
    if (identical(result$engine, "meme")) {
      shared_background <- motif_write_shared_background(query_sequences, file.path(package_dir, "shared_background.txt"))
      extra_files <- c(extra_files, "shared_background.txt")
    }
  } else {
    query_regions <- motif_center_regions(motif_regions, width = result$sequence_window)
    query_input <- motif_extract_sequences(query_regions, result$species, prefix = "query")
    query_regions <- query_input$regions
    query_sequences <- query_input$sequences
    control_regions <- motif_shuffle_regions(query_regions, result$species)
    control_input <- motif_extract_sequences(control_regions, result$species, prefix = "control")
    control_regions <- control_input$regions
    control_sequences <- control_input$sequences
    query_filtered <- motif_filter_sequences_by_width(query_regions, query_sequences, result$motif_width, "query", result$engine)
    control_filtered <- motif_filter_sequences_by_width(control_regions, control_sequences, result$motif_width, "control", result$engine)
    query_regions <- query_filtered$regions
    control_regions <- control_filtered$regions
    query_sequences <- query_filtered$sequences
    control_sequences <- control_filtered$sequences
  }

  query_fasta <- file.path(package_dir, "query.fa")
  control_fasta <- file.path(package_dir, "control.fa")
  query_bed <- file.path(package_dir, "query_regions.bed")
  control_bed <- file.path(package_dir, "control_regions.bed")
  script_name <- motif_engine_run_script_name(result$engine)
  run_script <- file.path(package_dir, script_name)
  manifest_path <- file.path(package_dir, "manifest.json")
  output_dir_name <- motif_engine_output_dir_name(result$engine)
  command_preview <- motif_engine_command_preview(result$engine, result$sample_name, list(
    motif_width = result$motif_width,
    max_motifs = result$max_motifs,
    analysis_mode = result$analysis_mode
  ))

  Biostrings::writeXStringSet(query_sequences, filepath = query_fasta, format = "fasta")
  Biostrings::writeXStringSet(control_sequences, filepath = control_fasta, format = "fasta")
  motif_write_bed(query_regions, query_bed)
  motif_write_bed(control_regions, control_bed)

  motif_write_unix_lines(
    run_script,
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "cd -- \"$(dirname -- \"$0\")\"",
      "",
      sprintf("mkdir -p %s", output_dir_name),
      command_preview
    )
  )

  package_files <- c(
    basename(query_fasta),
    basename(control_fasta),
    basename(query_bed),
    basename(control_bed),
    basename(run_script),
    basename(manifest_path),
    extra_files
  )
  manifest <- motif_build_job_manifest(result, package_files, command_preview, script_name)
  manifest$sequenceWindowWidth <- result$sequence_window
  manifest$queryCount <- length(query_regions)
  manifest$controlCount <- length(control_regions)
  manifest$queryLengthRange <- range(Biostrings::width(query_sequences))
  manifest$rVersion <- R.version.string
  manifest$matching <- matching_summary
  manifest$sharedBackground <- shared_background
  manifest$backgroundStrategy <- if (matched) "Same chromosome, RNA strand, gene, exact width, centre feature, GC and disjoint feature fractions; first accepted control per peak" else switch(result$engine,
    meme = "MEME default: background estimated from query sequences",
    streme = "STREME default: shuffled query sequences",
    homer = "Chromosome-matched random genomic control.fa")
  manifest$orientation <- if (matched) "BED RNA strand, forward search only" else "Engine default (MEME/STREME: both strands; HOMER: RNA strand)"
  manifest$inputMd5 <- as.list(tools::md5sum(file.path(package_dir, c("query.fa", "control.fa", extra_files))))
  jsonlite::write_json(manifest, path = manifest_path, auto_unbox = TRUE, pretty = TRUE)

  list(
    package_dir = normalizePath(package_dir, winslash = "/", mustWork = FALSE),
    package_files = package_files,
    script_name = script_name,
    manifest_name = basename(manifest_path),
    query_count = length(query_regions),
    control_count = length(control_regions),
    command_preview = command_preview,
    manifest = manifest
  )
}
