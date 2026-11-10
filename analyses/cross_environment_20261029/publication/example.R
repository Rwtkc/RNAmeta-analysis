# Run from the RNAmeta manuscript working directory after installing RNAmeta.
library(RNAmeta)
files <- file.path("comparative_GSE55572_20261025/input",
  c("Control.bed", "METTL3 KD.bed", "WTAP KD.bed"))
labels <- c("Control", "METTL3 KD", "WTAP KD")
ref <- function(suffix) file.path(
  "cross_environment_20261029/input/annotation",
  paste0("hg38.", suffix))

annotation <- RNAmeta_get_preliminary_analysis_data(
  files, ref("txdb.sqlite"), ref("txlens.rda"),
  ref("gff.rda"))
RNAmeta_peak_distribution_plot(annotation, labels)
RNAmeta_get_translation_plot(
  get_translation_data(annotation), labels)

meta <- RNAmeta_meta_plot(files, ref("txdb.sqlite"),
  plot_tx_type = "mRNA", head_or_tail = FALSE,
  tx_promoter_length = 0, tx_tail_length = 0,
  set_group_name = labels)
RNAmeta_meta_plot_export(meta, "mRNA.pdf")
