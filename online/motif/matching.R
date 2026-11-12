motif_matching_config <- function() {
  list(seed = 20261026L, controls_per_peak = 5L, gc_tolerance = 0.03,
    feature_fraction_tolerance = 0.10, rounds = 10L, candidates_per_round = 30L,
    engine_control = "First accepted background per matched peak; exact 1:1 pairing",
    feature_priority = c("CDS", "UTR3", "UTR5", "exon_other", "intron", "intergenic"))
}

motif_match_background <- function(query, sequences, annotation, genome, exclusion, progress = NULL) {
  config <- motif_matching_config()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(config$seed)
  qa <- motif_matching_features(query, annotation)
  gene_index <- motif_matching_genes(query, annotation)
  composition <- motif_matching_composition(sequences)
  eligible <- !is.na(gene_index) & qa$label != "intergenic" & composition$unambiguous
  domains <- motif_matching_domains(annotation, unique(gene_index[eligible]))
  chosen <- vector("list", length(query)); selected <- integer(length(query))
  excluded_index <- GenomicRanges::GNCList(exclusion)
  chromosome_lengths <- motif_reference_lengths(genome)
  for (round in seq_len(config$rounds)) {
    active <- which(eligible & selected < config$controls_per_peak)
    if (!length(active)) break
    sampled <- lapply(active, function(i) {
      domain <- domains[[paste(gene_index[[i]], qa$label[[i]], sep = ":")]]
      if (is.null(domain) || !nrow(domain)) return(NULL)
      widths <- domain$end - domain$start + 1L
      draws <- sample.int(sum(widths), config$candidates_per_round, replace = TRUE)
      offsets <- c(0, cumsum(widths))
      segment <- findInterval(draws - 1L, offsets)
      data.table::data.table(query = i, center = domain$start[segment] + draws - 1L - offsets[segment])
    })
    candidates <- unique(data.table::rbindlist(sampled))
    if (!nrow(candidates)) next
    qi <- candidates$query
    starts <- candidates$center - (IRanges::width(query)[qi] - 1L) %/% 2L
    regions <- GenomicRanges::GRanges(GenomicRanges::seqnames(query)[qi],
      IRanges::IRanges(starts, width = IRanges::width(query)[qi]), strand = GenomicRanges::strand(query)[qi])
    valid <- starts >= 1L & IRanges::end(regions) <= chromosome_lengths[as.character(GenomicRanges::seqnames(regions))]
    regions <- regions[valid]; candidates <- candidates[valid]
    valid <- !IRanges::overlapsAny(regions, excluded_index, ignore.strand = TRUE)
    regions <- regions[valid]; candidates <- candidates[valid]
    if (!length(regions)) next
    ca <- motif_matching_features(regions, annotation)
    difference <- abs(ca$fractions - qa$fractions[candidates$query, , drop = FALSE])
    valid <- ca$label == qa$label[candidates$query] & apply(difference, 1L, max) <= config$feature_fraction_tolerance
    regions <- regions[valid]; candidates <- candidates[valid]
    difference <- difference[valid, , drop = FALSE]
    if (!length(regions)) next
    candidate_sequences <- motif_reference_sequences(genome, regions)
    cc <- motif_matching_composition(candidate_sequences)
    valid <- cc$unambiguous & abs(cc$gc - composition$gc[candidates$query]) <= config$gc_tolerance
    regions <- regions[valid]; candidates <- candidates[valid]
    candidate_sequences <- candidate_sequences[valid]; cc <- cc[valid, , drop = FALSE]
    difference <- difference[valid, , drop = FALSE]
    if (!length(regions)) next
    candidates$start <- IRanges::start(regions)
    candidates$end <- IRanges::end(regions)
    candidates$gc <- cc$gc
    candidates$max_feature_difference <- apply(difference, 1L, max)
    candidates$round <- round
    for (i in unique(candidates$query)) {
      options <- which(candidates$query == i)
      prior <- if (is.null(chosen[[i]])) integer() else chosen[[i]]$start
      options <- options[!candidates$start[options] %in% prior]
      options <- head(options, config$controls_per_peak - selected[[i]])
      if (!length(options)) next
      rows <- data.table::copy(candidates[options])
      rows$sequence <- as.character(candidate_sequences[options])
      chosen[[i]] <- data.table::rbindlist(list(chosen[[i]], rows))
      selected[[i]] <- selected[[i]] + length(options)
    }
    motif_progress_update(progress, 82 + round / config$rounds * 6,
      sprintf("Matching background: %d / %d peaks", sum(selected > 0L), length(query)))
  }
  controls <- data.table::rbindlist(chosen)
  query_table <- data.frame(query = seq_along(query), name = names(query),
    chrom = as.character(GenomicRanges::seqnames(query)), start = IRanges::start(query) - 1L,
    end = IRanges::end(query), strand = as.character(GenomicRanges::strand(query)),
    width = IRanges::width(query), gc = composition$gc, feature = qa$label,
    gene_key = annotation$genes$gene_key[gene_index], controls = selected,
    reason = ifelse(is.na(gene_index), "no_gene_at_center",
      ifelse(!composition$unambiguous, "ambiguous_query_sequence",
        ifelse(selected == 0L, "no_match_within_fixed_tolerances", "matched"))))
  if (nrow(controls)) {
    controls$gene_key <- annotation$genes$gene_key[gene_index[controls$query]]
    controls$feature <- qa$label[controls$query]
  }
  list(queries = query_table, controls = controls, config = config)
}

motif_write_shared_background <- function(sequences, path) {
  counts <- colSums(Biostrings::alphabetFrequency(sequences, baseOnly = TRUE)[, c("A", "C", "G", "T"), drop = FALSE])
  frequencies <- (counts + 1) / (sum(counts) + 4)
  motif_write_unix_lines(path, c("# Shared zero-order query background; add-one pseudocount per base",
    paste(names(frequencies), sprintf("%.17g", frequencies))))
  list(source = "query.fa", order = 0L, pseudocountPerBase = 1L,
    sharedByQueryAndControl = TRUE, frequencies = as.list(frequencies), file = basename(path))
}

motif_prepare_matched_input <- function(regions, result, annotation_bundle, exclusion, package_dir, progress = NULL) {
  genome <- motif_load_seqdb(result$species)
  chromosome_lengths <- motif_reference_lengths(genome)
  chrom <- as.character(GenomicRanges::seqnames(regions))
  valid_genome <- chrom %in% names(chromosome_lengths) & IRanges::start(regions) >= 1L &
    IRanges::end(regions) <= chromosome_lengths[chrom]
  valid_genome[is.na(valid_genome)] <- FALSE
  known_strand <- as.character(GenomicRanges::strand(regions)) %in% c("+", "-")
  wide_enough <- IRanges::width(regions) >= result$motif_width
  audit <- data.frame(input_index = seq_along(regions), chrom = chrom,
    start = IRanges::start(regions) - 1L, end = IRanges::end(regions),
    strand = as.character(GenomicRanges::strand(regions)),
    reason = ifelse(!known_strand, "unknown_RNA_strand",
      ifelse(!valid_genome, "outside_selected_genome",
        ifelse(!wide_enough, "shorter_than_motif_width", "pending_match"))))
  write_audit <- function() utils::write.table(audit, file.path(package_dir, "input_audit.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  keep <- which(audit$reason == "pending_match")
  if (!length(keep)) {
    write_audit()
    stop("No full peaks are eligible. Matched analysis requires BED RNA strands (+/-), valid genome coordinates, and peaks at least as wide as the motif.")
  }
  query <- regions[keep]
  names(query) <- sprintf("query_%d", keep)
  sequences <- motif_reference_sequences(genome, query)
  names(sequences) <- names(query)
  annotation <- motif_matching_annotation(annotation_bundle$txdb, names(chromosome_lengths))
  exclusion <- exclusion[as.character(GenomicRanges::seqnames(exclusion)) %in% names(chromosome_lengths)]
  exclusion <- GenomicRanges::reduce(exclusion, ignore.strand = TRUE)
  matched <- motif_match_background(query, sequences, annotation, genome, exclusion, progress)
  audit$reason[keep] <- matched$queries$reason
  write_audit()
  utils::write.table(matched$queries, file.path(package_dir, "matching_queries.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  controls <- matched$controls
  if (!nrow(controls)) stop("No non-peak background matched the fixed gene, feature and GC criteria. See input_audit.tsv for excluded peaks.")
  controls$chrom <- as.character(GenomicRanges::seqnames(query))[controls$query]
  controls$strand <- as.character(GenomicRanges::strand(query))[controls$query]
  controls$control_id <- sprintf("control_%d", seq_len(nrow(controls)))
  pair_table <- controls[, setdiff(names(controls), "sequence"), with = FALSE]
  utils::write.table(pair_table, file.path(package_dir, "matching_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  first <- controls[!duplicated(controls$query)]
  control_regions <- GenomicRanges::GRanges(first$chrom, IRanges::IRanges(first$start, first$end), strand = first$strand)
  names(control_regions) <- first$control_id
  control_sequences <- Biostrings::DNAStringSet(first$sequence)
  names(control_sequences) <- names(control_regions)
  query <- query[first$query]; sequences <- sequences[first$query]
  if (!all(IRanges::width(query) == Biostrings::width(control_sequences))) stop("Matched motif input lost length pairing.")
  summary <- list(inputCount = length(regions), eligibleCount = length(keep),
    matchedCount = length(query), matchedFraction = length(query) / length(regions),
    excludedReasons = as.list(table(audit$reason[audit$reason != "matched"])),
    config = matched$config, exclusionScope = "Union of all saved BED files, before subset selection; either strand",
    annotationMd5 = as.list(tools::md5sum(annotation_bundle$txdb_path)),
    selection = "No motif counts used; only peaks with at least one accepted background enter discovery",
    orientation = "RNA strand from BED; unknown strands excluded; no reverse-complement search")
  list(query_regions = query, query_sequences = sequences, control_regions = control_regions,
    control_sequences = control_sequences, matching_summary = summary,
    extra_files = c("input_audit.tsv", "matching_queries.tsv", "matching_pairs.tsv"))
}

motif_matching_annotation <- function(txdb, valid_chromosomes) {
  trim <- function(gr) {
    S4Vectors::mcols(gr) <- NULL
    gr <- gr[as.character(GenomicRanges::seqnames(gr)) %in% valid_chromosomes &
      as.character(GenomicRanges::strand(gr)) %in% c("+", "-")]
    GenomicRanges::reduce(gr, ignore.strand = FALSE)
  }
  flatten <- function(x) trim(unlist(x, use.names = FALSE))
  genes <- GenomicFeatures::genes(txdb)
  genes <- genes[as.character(GenomicRanges::seqnames(genes)) %in% valid_chromosomes]
  genes$gene_key <- names(genes)
  if (!length(genes)) stop("Matched background requires annotated genes on the selected genome.")
  raw <- list(CDS = flatten(GenomicFeatures::cdsBy(txdb, by = "tx")),
    UTR3 = flatten(GenomicFeatures::threeUTRsByTranscript(txdb)),
    UTR5 = flatten(GenomicFeatures::fiveUTRsByTranscript(txdb)),
    exon_other = trim(GenomicFeatures::exons(txdb)), intron = trim(genes))
  used <- raw[[1]][FALSE]
  features <- lapply(names(raw), function(label) {
    part <- GenomicRanges::setdiff(raw[[label]], used, ignore.strand = FALSE)
    used <<- GenomicRanges::reduce(c(used, raw[[label]]), ignore.strand = FALSE)
    part$feature <- rep(label, length(part))
    part
  })
  masks <- do.call(c, unname(features))
  list(genes = genes, masks = masks, index = GenomicRanges::GNCList(masks),
    types = c(names(raw), "intergenic"))
}

motif_matching_centers <- function(gr) {
  GenomicRanges::GRanges(GenomicRanges::seqnames(gr),
    IRanges::IRanges(IRanges::start(gr) + (IRanges::width(gr) - 1L) %/% 2L, width = 1L),
    strand = GenomicRanges::strand(gr))
}

motif_matching_features <- function(gr, annotation) {
  types <- annotation$types
  h <- GenomicRanges::findOverlaps(motif_matching_centers(gr), annotation$index, ignore.strand = FALSE)
  if (anyDuplicated(S4Vectors::queryHits(h))) stop("Matched background feature masks must be disjoint.")
  labels <- rep("intergenic", length(gr))
  labels[S4Vectors::queryHits(h)] <- annotation$masks$feature[S4Vectors::subjectHits(h)]
  fractions <- matrix(0, nrow = length(gr), ncol = length(types), dimnames = list(NULL, types))
  h <- GenomicRanges::findOverlaps(gr, annotation$index, ignore.strand = FALSE)
  qi <- S4Vectors::queryHits(h); si <- S4Vectors::subjectHits(h)
  if (length(qi)) {
    spans <- pmin(IRanges::end(gr)[qi], IRanges::end(annotation$masks)[si]) -
      pmax(IRanges::start(gr)[qi], IRanges::start(annotation$masks)[si]) + 1L
    tab <- data.table::data.table(query = qi, type = annotation$masks$feature[si], span = spans)
    tab <- tab[, .(span = sum(span)), by = .(query, type)]
    fractions[cbind(tab$query, match(tab$type, types))] <- tab$span / IRanges::width(gr)[tab$query]
  }
  fractions[, "intergenic"] <- pmax(0, 1 - rowSums(fractions))
  list(label = labels, fractions = fractions)
}

motif_matching_genes <- function(gr, annotation) {
  genes <- annotation$genes
  h <- GenomicRanges::findOverlaps(motif_matching_centers(gr), genes, ignore.strand = FALSE)
  qi <- S4Vectors::queryHits(h); si <- S4Vectors::subjectHits(h)
  tab <- data.table::data.table(query = qi, gene = si,
    overlap = pmin(IRanges::end(gr)[qi], IRanges::end(genes)[si]) -
      pmax(IRanges::start(gr)[qi], IRanges::start(genes)[si]) + 1L,
    gene_key = genes$gene_key[si])
  data.table::setorder(tab, query, -overlap, gene_key)
  tab <- tab[!duplicated(query)]
  answer <- rep(NA_integer_, length(gr)); answer[tab$query] <- tab$gene
  answer
}

motif_matching_domains <- function(annotation, used_gene_indices) {
  genes <- annotation$genes[used_gene_indices]
  h <- GenomicRanges::findOverlaps(genes, annotation$index, ignore.strand = FALSE)
  qi <- S4Vectors::queryHits(h); si <- S4Vectors::subjectHits(h)
  out <- data.table::data.table(gene = used_gene_indices[qi], feature = annotation$masks$feature[si],
    start = pmax(IRanges::start(genes)[qi], IRanges::start(annotation$masks)[si]),
    end = pmin(IRanges::end(genes)[qi], IRanges::end(annotation$masks)[si]))
  split(out, paste(out$gene, out$feature, sep = ":"))
}

motif_matching_composition <- function(sequences) {
  counts <- Biostrings::alphabetFrequency(sequences, baseOnly = TRUE)
  data.frame(gc = (counts[, "G"] + counts[, "C"]) / Biostrings::width(sequences),
    unambiguous = rowSums(counts[, c("A", "C", "G", "T"), drop = FALSE]) == Biostrings::width(sequences))
}
