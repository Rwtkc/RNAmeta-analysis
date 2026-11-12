splicesite_default_controls <- function() {
  list(
    flank_size = 100
  )
}

site_profile_prepare_splicesite_windows <- function(species, exons_by_tx, flank_size = 100) {
  site_profile_cached_boundary_windows("splicesite", species, flank_size, function() {
    exons_flat <- BiocGenerics::unlist(exons_by_tx)
    exons_flat <- exons_flat[IRanges::width(exons_flat) > flank_size]
    transcript_ranges <- BiocGenerics::unlist(range(exons_by_tx))

    ss5p <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(exons_flat),
      ranges = IRanges::IRanges(
        ifelse(as.character(GenomicRanges::strand(exons_flat)) == "+", BiocGenerics::end(exons_flat), BiocGenerics::start(exons_flat)),
        ifelse(as.character(GenomicRanges::strand(exons_flat)) == "+", BiocGenerics::end(exons_flat), BiocGenerics::start(exons_flat))
      ),
      strand = GenomicRanges::strand(exons_flat)
    )
    names(ss5p) <- names(exons_flat)

    ss3p <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(exons_flat),
      ranges = IRanges::IRanges(
        ifelse(as.character(GenomicRanges::strand(exons_flat)) == "+", BiocGenerics::start(exons_flat), BiocGenerics::end(exons_flat)),
        ifelse(as.character(GenomicRanges::strand(exons_flat)) == "+", BiocGenerics::start(exons_flat), BiocGenerics::end(exons_flat))
      ),
      strand = GenomicRanges::strand(exons_flat)
    )
    names(ss3p) <- names(exons_flat)

    transcript5p <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(transcript_ranges),
      ranges = IRanges::IRanges(
        ifelse(as.character(GenomicRanges::strand(transcript_ranges)) == "+", BiocGenerics::end(transcript_ranges), BiocGenerics::start(transcript_ranges)),
        ifelse(as.character(GenomicRanges::strand(transcript_ranges)) == "+", BiocGenerics::end(transcript_ranges), BiocGenerics::start(transcript_ranges))
      ),
      strand = GenomicRanges::strand(transcript_ranges)
    )

    transcript3p <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(transcript_ranges),
      ranges = IRanges::IRanges(
        ifelse(as.character(GenomicRanges::strand(transcript_ranges)) == "+", BiocGenerics::start(transcript_ranges), BiocGenerics::end(transcript_ranges)),
        ifelse(as.character(GenomicRanges::strand(transcript_ranges)) == "+", BiocGenerics::start(transcript_ranges), BiocGenerics::end(transcript_ranges))
      ),
      strand = GenomicRanges::strand(transcript_ranges)
    )

    ss5p <- ss5p[setdiff(seq_along(ss5p), unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(ss5p, transcript5p))))]
    ss3p <- ss3p[setdiff(seq_along(ss3p), unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(ss3p, transcript3p))))]

    windows_5p <- GenomicRanges::resize(
      GenomicRanges::resize(ss5p, width = flank_size * 2 + 1L, fix = "center"),
      width = flank_size * 3 + 1L,
      fix = "start"
    )
    windows_3p <- GenomicRanges::resize(
      GenomicRanges::resize(ss3p, width = flank_size * 4 + 1L, fix = "center"),
      width = flank_size * 3 + 1L,
      fix = "start"
    )

    list(
      windows_5p = windows_5p,
      windows_3p = windows_3p,
      offset_5p = flank_size + 1L,
      offset_3p = flank_size * 2 + 1L
    )
  })
}
