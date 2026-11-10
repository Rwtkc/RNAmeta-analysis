gene_matrix_intersection_rows <- function(gene_sets, sample_labels, sample_metadata = NULL) {
  rows <- list()

  for (set_size in seq_along(sample_labels)) {
    combinations <- combn(sample_labels, set_size, simplify = FALSE)

    for (combo in combinations) {
      genes <- if (length(combo) == 1) {
        gene_sets[[combo[[1]]]]
      } else {
        Reduce(intersect, gene_sets[combo])
      }

      rows[[length(rows) + 1]] <- list(
        sets = unname(as.list(combo)),
        originalSets = unname(lapply(combo, function(label) rnameta_lookup_original_name(sample_metadata, label, fallback = label))),
        size = length(genes),
        label = if (length(combo) == 1) combo[[1]] else "",
        genes = unname(as.list(sort(genes)))
      )
    }
  }

  rows
}

gene_matrix_build_chart_payload <- function(gene_sets, sample_labels, sample_metadata = NULL) {
  if (!length(gene_sets)) {
    return(NULL)
  }

  colors <- gene_matrix_venn_palette()
  color_map <- setNames(colors[seq_along(sample_labels)], sample_labels)

  list(
    title = "Gene Matrix",
    sampleLabels = unname(sample_labels),
    sampleMetadata = lapply(sample_labels, function(label) {
      list(
        name = label,
        originalName = rnameta_lookup_original_name(sample_metadata, label, fallback = label)
      )
    }),
    colors = as.list(color_map),
    intersections = gene_matrix_intersection_rows(gene_sets, sample_labels, sample_metadata = sample_metadata)
  )
}
