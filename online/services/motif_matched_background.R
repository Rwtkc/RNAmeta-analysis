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
