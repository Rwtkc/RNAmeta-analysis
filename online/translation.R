translation_default_controls <- function() {
  list(
    flank_size = 100
  )
}

site_profile_prepare_translation_windows <- function(species, txlens, exons, flank_size = 100) {
  site_profile_cached_boundary_windows("translation", species, flank_size, function() {
    txlens_dt <- data.table::as.data.table(txlens)

    use_cds_only <- identical(species, "sce_R64") && !any(
      txlens_dt$utr5_len > flank_size & txlens_dt$utr3_len > flank_size &
        txlens_dt$cds_len > flank_size, na.rm = TRUE
    )
    if (use_cds_only) {
      usable <- txlens_dt[cds_len > flank_size, .(tx_name, utr5_len, cds_len, utr3_len)]
      tss_ranges <- usable[, .(
        seqnames = tx_name,
        start = utr5_len + 1,
        end = utr5_len + flank_size,
        strand = "*"
      )]
      tes_ranges <- usable[, .(
        seqnames = tx_name,
        start = utr5_len + cds_len - flank_size + 1,
        end = utr5_len + cds_len,
        strand = "*"
      )]
      x_domain <- c(1, flank_size)
      x_positions <- seq_len(flank_size)
    } else {
      usable <- txlens_dt[utr5_len > flank_size & utr3_len > flank_size & cds_len > flank_size, .(tx_name, utr5_len, cds_len, utr3_len)]
      tss_ranges <- usable[, .(
        seqnames = tx_name,
        start = utr5_len + 1 - flank_size,
        end = utr5_len + 1 + flank_size,
        strand = "*"
      )]
      tes_ranges <- usable[, .(
        seqnames = tx_name,
        start = utr5_len + cds_len - flank_size,
        end = utr5_len + cds_len + flank_size,
        strand = "*"
      )]
      x_domain <- c(-flank_size, flank_size)
      x_positions <- seq(-flank_size, flank_size)
    }

    tss_windows <- GenomicRanges::makeGRangesFromDataFrame(tss_ranges, keep.extra.columns = TRUE)
    tes_windows <- GenomicRanges::makeGRangesFromDataFrame(tes_ranges, keep.extra.columns = TRUE)
    names(tss_windows) <- as.character(tss_ranges$seqnames)
    names(tes_windows) <- as.character(tes_ranges$seqnames)
    GenomeInfoDb::seqlevels(tss_windows) <- unique(c(GenomeInfoDb::seqlevels(tss_windows), names(exons)))
    GenomeInfoDb::seqlevels(tes_windows) <- unique(c(GenomeInfoDb::seqlevels(tes_windows), names(exons)))

    list(
      tss_windows = tss_windows,
      tes_windows = tes_windows,
      x_positions = x_positions,
      x_domain = x_domain,
      rel_offset = if (use_cds_only) 0L else flank_size + 1L,
      flank_size = flank_size
    )
  })
}
