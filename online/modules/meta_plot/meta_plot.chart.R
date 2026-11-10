meta_plot_build_component_payload <- function(component_width, controls, pos) {
  component_structure <- data.frame(
    width = component_width,
    comp = names(component_width),
    stringsAsFactors = FALSE
  )
  component_structure <- component_structure[component_structure$width > 0, , drop = FALSE]

  if (!nrow(component_structure)) {
    return(list(trackBase = pos$rna_lgd_bl, labelY = pos$rna_comp_text - 0.035, segments = list(), separators = list()))
  }

  component_structure$end <- cumsum(component_structure$width)
  component_structure$start <- c(0, head(component_structure$end, -1)) + 0.001
  component_structure$mid <- (component_structure$start + component_structure$end) / 2
  component_structure$label <- component_structure$comp
  component_structure$label[component_structure$comp == "promoter"] <- if (identical(controls$headOrtail, "TRUE")) sprintf("Promoter(%skb)", controls$txpromoterLength / 1000) else ""
  component_structure$label[component_structure$comp == "tail"] <- if (identical(controls$headOrtail, "TRUE")) sprintf("Tail(%skb)", controls$txtailLength / 1000) else ""
  component_structure$label[component_structure$comp == "utr5"] <- "5'UTR"
  component_structure$label[component_structure$comp == "utr3"] <- "3'UTR"
  component_structure$label[component_structure$comp == "cds"] <- "CDS"
  component_structure$label[component_structure$comp == "ncrna"] <- "ncRNA"
  component_structure$label[component_structure$comp == "rna"] <- "RNA"
  component_structure$alpha <- ifelse(
    component_structure$comp %in% c("cds", "ncrna", "rna"),
    c(cds = 0.272, ncrna = 0.2, rna = 0.2)[component_structure$comp],
    0.99
  )
  component_structure$height <- ifelse(
    component_structure$comp == "cds",
    pos$rna_lgd_h_cds,
    ifelse(
      component_structure$comp %in% c("ncrna", "rna"),
      pos$rna_lgd_h_ncrna,
      ifelse(component_structure$comp %in% c("utr5", "utr3"), pos$rna_lgd_h_utr, pos$rna_lgd_h_flank)
    )
  )

  list(
    trackBase = pos$rna_lgd_bl,
    labelY = pos$rna_comp_text - 0.035,
    segments = lapply(seq_len(nrow(component_structure)), function(i) {
      list(
        start = component_structure$start[[i]],
        end = component_structure$end[[i]],
        mid = component_structure$mid[[i]],
        label = component_structure$label[[i]],
        alpha = component_structure$alpha[[i]],
        height = component_structure$height[[i]],
        component = component_structure$comp[[i]]
      )
    }),
    separators = lapply(head(component_structure$end, -1), function(value) {
      list(x = unname(value), y1 = pos$fig_top, y2 = pos$rna_lgd_bl)
    })
  )
}

meta_plot_build_chart_payload <- function(density_ci, component_width, controls, series_metadata = NULL) {
  peak <- if (identical(controls$enableCI, "TRUE") && "confidenceUp" %in% names(density_ci)) {
    max(density_ci$confidenceUp, na.rm = TRUE)
  } else {
    max(density_ci$density, na.rm = TRUE)
  }
  pos <- .generate_pos_para(peak)
  title <- c(tx = "Distribution on exon", mrna = "Distribution on mRNA", ncrna = "Distribution on lncRNA")[[controls$pltTxType]]
  palette <- c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF")
  groups <- split(density_ci, density_ci$group)

  list(
    title = title,
    yLabel = "Frequency",
    xDomain = c(0, 1),
    yDomain = c(pos$fig_bottom, pos$fig_top),
    showCI = identical(controls$enableCI, "TRUE") && "confidenceUp" %in% names(density_ci),
    components = meta_plot_build_component_payload(component_width, controls, pos),
    series = lapply(seq_along(groups), function(index) {
      group_data <- groups[[index]]
      group_name <- names(groups)[[index]]
      metadata_index <- NULL
      if (!is.null(series_metadata) && length(series_metadata)) {
        metadata_index <- which(vapply(series_metadata, function(entry) identical(entry$display_name, group_name), logical(1)))[1]
      }
      original_name <- if (!is.na(metadata_index) && length(metadata_index) == 1L) series_metadata[[metadata_index]]$original_name else group_name
      list(
        name = group_name,
        originalName = original_name,
        color = palette[((index - 1) %% length(palette)) + 1],
        values = lapply(seq_len(nrow(group_data)), function(i) {
          list(
            x = group_data$x[[i]],
            density = group_data$density[[i]],
            confidenceDown = if ("confidenceDown" %in% names(group_data)) group_data$confidenceDown[[i]] else NULL,
            confidenceUp = if ("confidenceUp" %in% names(group_data)) group_data$confidenceUp[[i]] else NULL
          )
        })
      )
    })
  )
}
