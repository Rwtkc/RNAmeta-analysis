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
