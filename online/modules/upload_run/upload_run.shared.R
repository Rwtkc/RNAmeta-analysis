upload_run_annotation_data_root <- function() {
  app_path("data")
}

upload_run_database_file_metadata <- function(path) {
  file_name <- basename(path)
  name_parts <- strsplit(sub("\\.bed$", "", file_name), " - ", fixed = TRUE)[[1]]
  name_parts <- trimws(name_parts)
  caller <- if (length(name_parts) >= 1L) tail(name_parts, 1L) else "Unknown"
  accession <- if (length(name_parts) >= 3L) name_parts[[3L]] else "Unknown"
  context_parts <- if (length(name_parts) > 4L) name_parts[4L:(length(name_parts) - 1L)] else character()
  sample <- if (length(context_parts)) context_parts[[1L]] else "Unspecified"
  treatment <- if (length(context_parts) > 1L) paste(context_parts[-1L], collapse = " - ") else "Unspecified"

  list(
    key = file_name,
    label = file_name,
    file_path = path,
    accession = accession,
    sample = sample,
    treatment = treatment,
    caller = caller
  )
}

upload_run_database_source <- function(
  species,
  directory,
  assembly,
  source_id = "m6a_atlas_v2",
  source_label = "m6A-Atlas v2.0 (MeRIP-seq)",
  modification = "m6A",
  assay = "MeRIP-seq",
  source_url = "https://rnamd.org/m6a/merip.php",
  citation = "m6A-Atlas v2.0"
) {
  bed_paths <- sort(list.files(directory, pattern = "\\.bed$", full.names = TRUE))
  if (!length(bed_paths)) {
    return(NULL)
  }

  list(
    id = source_id,
    label = source_label,
    species = species,
    info = list(
      modification = modification,
      assay = assay,
      assembly = assembly,
      source_url = source_url,
      citation = citation
    ),
    files = lapply(bed_paths, upload_run_database_file_metadata)
  )
}

upload_run_database_catalog <- function() {
  sources <- list(
    upload_run_database_source(
      species = "ara_TAIR10",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "ara_TAIR10"),
      assembly = "TAIR10"
    ),
    upload_run_database_source(
      species = "cel_WBcel235",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "cel_WBcel235"),
      assembly = "WBcel235"
    ),
    upload_run_database_source(
      species = "dme_BDGP6",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "dme_BDGP6"),
      assembly = "BDGP6"
    ),
    upload_run_database_source(
      species = "mm10",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "mm10"),
      assembly = "GRCm38 / mm10"
    ),
    upload_run_database_source(
      species = "mml.Mmul_10",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "mml.Mmul_10"),
      assembly = "Mmul_10"
    ),
    upload_run_database_source(
      species = "osa_IRGSP_1",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "osa_IRGSP_1"),
      assembly = "IRGSP 1.0"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "hg38"),
      assembly = "GRCh38 / hg38"
    ),
    upload_run_database_source(
      species = "ptc.Pop_tri_v3",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "ptc.Pop_tri_v3"),
      assembly = "Pop_tri_v3"
    ),
    upload_run_database_source(
      species = "ptr.Pan_tro_3.0",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "ptr.Pan_tro_3.0"),
      assembly = "Pan_tro_3.0"
    ),
    upload_run_database_source(
      species = "rn6",
      directory = app_path("data", "database_sources", "m6a_atlas_v2", "bed", "rn6"),
      assembly = "Rnor_6.0 / rn6"
    ),
    upload_run_database_source(
      species = "ara_TAIR10",
      directory = app_path("data", "database_sources", "repic_m6a", "bed", "ara_TAIR10"),
      assembly = "TAIR10",
      source_id = "repic_m6a",
      source_label = "REPIC m6A (MeRIP-seq)",
      assay = "MeRIP-seq / m6A-seq",
      source_url = "https://repicmod.uchicago.edu/repic/download.php",
      citation = "REPIC"
    ),
    upload_run_database_source(
      species = "ptr.Pan_tro_3.0",
      directory = app_path("data", "database_sources", "repic_m6a", "bed", "ptr.Pan_tro_3.0"),
      assembly = "Pan_tro_3.0",
      source_id = "repic_m6a",
      source_label = "REPIC m6A (MeRIP-seq)",
      assay = "MeRIP-seq / m6A-seq",
      source_url = "https://repicmod.uchicago.edu/repic/download.php",
      citation = "REPIC"
    ),
    upload_run_database_source(
      species = "dme_BDGP6",
      directory = app_path("data", "database_sources", "repic_m6a", "bed", "dme_BDGP6"),
      assembly = "BDGP6",
      source_id = "repic_m6a",
      source_label = "REPIC m6A (MeRIP-seq)",
      assay = "MeRIP-seq / m6A-seq",
      source_url = "https://repicmod.uchicago.edu/repic/download.php",
      citation = "REPIC"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "repic_m6a", "bed", "hg38"),
      assembly = "GRCh38 / hg38",
      source_id = "repic_m6a",
      source_label = "REPIC m6A (MeRIP-seq)",
      assay = "MeRIP-seq / m6A-seq",
      source_url = "https://repicmod.uchicago.edu/repic/download.php",
      citation = "REPIC"
    ),
    upload_run_database_source(
      species = "mm10",
      directory = app_path("data", "database_sources", "repic_m6a", "bed", "mm10"),
      assembly = "GRCm38 / mm10",
      source_id = "repic_m6a",
      source_label = "REPIC m6A (MeRIP-seq)",
      assay = "MeRIP-seq / m6A-seq",
      source_url = "https://repicmod.uchicago.edu/repic/download.php",
      citation = "REPIC"
    ),
    upload_run_database_source(
      species = "rn6",
      directory = app_path("data", "database_sources", "repic_m6a", "bed", "rn6"),
      assembly = "Rnor_6.0 / rn6",
      source_id = "repic_m6a",
      source_label = "REPIC m6A (MeRIP-seq)",
      assay = "MeRIP-seq / m6A-seq",
      source_url = "https://repicmod.uchicago.edu/repic/download.php",
      citation = "REPIC"
    ),
    upload_run_database_source(
      species = "sce_R64",
      directory = app_path("data", "database_sources", "repic_m6a", "bed", "sce_R64"),
      assembly = "R64-1-1",
      source_id = "repic_m6a",
      source_label = "REPIC m6A (MeRIP-seq)",
      assay = "MeRIP-seq / m6A-seq",
      source_url = "https://repicmod.uchicago.edu/repic/download.php",
      citation = "REPIC"
    ),
    upload_run_database_source(
      species = "ara_TAIR10",
      directory = app_path("data", "database_sources", "rmbase_v3", "bed", "ara_TAIR10"),
      assembly = "TAIR10",
      source_id = "rmbase_v3",
      source_label = "RMBase v3.0 (RNA modifications)",
      modification = "RNA modifications",
      assay = "Integrated RNA modification sites",
      source_url = "https://bioinformaticsscience.cn/rmbase/download.php",
      citation = "RMBase v3.0"
    ),
    upload_run_database_source(
      species = "dme_BDGP6",
      directory = app_path("data", "database_sources", "rmbase_v3", "bed", "dme_BDGP6"),
      assembly = "BDGP6",
      source_id = "rmbase_v3",
      source_label = "RMBase v3.0 (RNA modifications)",
      modification = "RNA modifications",
      assay = "Integrated RNA modification sites",
      source_url = "https://bioinformaticsscience.cn/rmbase/download.php",
      citation = "RMBase v3.0"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "rmbase_v3", "bed", "hg38"),
      assembly = "GRCh38 / hg38",
      source_id = "rmbase_v3",
      source_label = "RMBase v3.0 (RNA modifications)",
      modification = "RNA modifications",
      assay = "Integrated RNA modification sites",
      source_url = "https://bioinformaticsscience.cn/rmbase/download.php",
      citation = "RMBase v3.0"
    ),
    upload_run_database_source(
      species = "mm10",
      directory = app_path("data", "database_sources", "rmbase_v3", "bed", "mm10"),
      assembly = "GRCm38 / mm10",
      source_id = "rmbase_v3",
      source_label = "RMBase v3.0 (RNA modifications)",
      modification = "RNA modifications",
      assay = "Integrated RNA modification sites",
      source_url = "https://bioinformaticsscience.cn/rmbase/download.php",
      citation = "RMBase v3.0"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "geo_rbp", "bed", "hg38"),
      assembly = "GRCh38 / hg38",
      source_id = "geo_rbp",
      source_label = "GEO / ENCORE (RBP binding sites)",
      modification = "RBP binding",
      assay = "eCLIP",
      source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315407",
      citation = "ENCORE; DOI: 10.1016/j.cels.2026.101588"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "encode_eclip", "bed", "hg38"),
      assembly = "GRCh38 / hg38",
      source_id = "encode_eclip",
      source_label = "ENCODE eCLIP (RBP binding peaks)",
      modification = "RBP binding",
      assay = "eCLIP",
      source_url = "https://www.encodeproject.org/",
      citation = "ENCODE"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "m6aconquer", "bed", "hg38"),
      assembly = "GRCh38 / hg38",
      source_id = "m6aconquer",
      source_label = "m6AConquer (orthogonally validated m6A sites)",
      modification = "m6A",
      assay = "Orthogonal m6A site integration",
      source_url = "https://rnamd.org/m6aconquer/download.html",
      citation = "m6AConquer"
    ),
    upload_run_database_source(
      species = "mm10",
      directory = app_path("data", "database_sources", "m6aconquer", "bed", "mm10"),
      assembly = "GRCm38 / mm10",
      source_id = "m6aconquer",
      source_label = "m6AConquer (orthogonally validated m6A sites)",
      modification = "m6A",
      assay = "Orthogonal m6A site integration",
      source_url = "https://rnamd.org/m6aconquer/download.html",
      citation = "m6AConquer"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "epicancer_m6a", "bed", "hg38"),
      assembly = "GRCh38 / hg38",
      source_id = "epicancer_m6a",
      source_label = "EpiCancer m6A sites",
      modification = "m6A",
      assay = "m6A-seq / MeRIP-seq",
      source_url = "https://ngdc.cncb.ac.cn/epicancer/download",
      citation = "EpiCancer"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "encori_rbp", "bed", "hg38"),
      assembly = "GRCh38 / hg38",
      source_id = "encori_rbp",
      source_label = "ENCORI (CLIP-seq RBP binding sites)",
      modification = "RBP binding",
      assay = "CLIP-seq",
      source_url = "https://rnasysu.com/encori/download.php",
      citation = "ENCORI / starBase"
    ),
    upload_run_database_source(
      species = "mm10",
      directory = app_path("data", "database_sources", "encori_rbp", "bed", "mm10"),
      assembly = "GRCm38 / mm10",
      source_id = "encori_rbp",
      source_label = "ENCORI (CLIP-seq RBP binding sites)",
      modification = "RBP binding",
      assay = "CLIP-seq",
      source_url = "https://rnasysu.com/encori/download.php",
      citation = "ENCORI / starBase"
    ),
    upload_run_database_source(
      species = "hg19",
      directory = app_path("data", "database_sources", "rbpspot", "bed", "hg19", "SetA"),
      assembly = "GRCh37 / hg19",
      source_id = "rbpspot_hg19_seta",
      source_label = "RBPSpot SetA (CLIP-seq RBP binding sites)",
      modification = "RBP binding",
      assay = "CLIP-seq",
      source_url = "https://scbb.ihbt.res.in/RBPSpot/data/",
      citation = "RBPSpot"
    ),
    upload_run_database_source(
      species = "hg19",
      directory = app_path("data", "database_sources", "rbpspot", "bed", "hg19", "SetB"),
      assembly = "GRCh37 / hg19",
      source_id = "rbpspot_hg19_setb",
      source_label = "RBPSpot SetB (CLIP-seq RBP binding sites)",
      modification = "RBP binding",
      assay = "CLIP-seq",
      source_url = "https://scbb.ihbt.res.in/RBPSpot/data/",
      citation = "RBPSpot"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "rbpspot", "bed", "hg38", "SetA"),
      assembly = "GRCh38 / hg38",
      source_id = "rbpspot_hg38_seta",
      source_label = "RBPSpot SetA (CLIP-seq RBP binding sites)",
      modification = "RBP binding",
      assay = "CLIP-seq",
      source_url = "https://scbb.ihbt.res.in/RBPSpot/data/",
      citation = "RBPSpot"
    ),
    upload_run_database_source(
      species = "hg38",
      directory = app_path("data", "database_sources", "rbpspot", "bed", "hg38", "SetB"),
      assembly = "GRCh38 / hg38",
      source_id = "rbpspot_hg38_setb",
      source_label = "RBPSpot SetB (CLIP-seq RBP binding sites)",
      modification = "RBP binding",
      assay = "CLIP-seq",
      source_url = "https://scbb.ihbt.res.in/RBPSpot/data/",
      citation = "RBPSpot"
    )
  )

  extension_path <- app_path("data", "database_sources", "catalog_extension.json")
  if (file.exists(extension_path)) {
    extension <- jsonlite::read_json(extension_path, simplifyVector = FALSE)$databases
    extension <- lapply(extension, function(source) {
      source$files <- lapply(source$files, function(file) {
        file$file_path <- app_path("data", "database_sources", source$id, "bed", source$species, file$key)
        file
      })
      source
    })
    sources <- c(sources, extension)
  }

  Filter(Negate(is.null), sources)
}

upload_run_database_control_catalog <- function(catalog = upload_run_database_catalog(), include_files = TRUE) {
  lapply(catalog, function(source) {
    entry <- list(
      id = source$id,
      label = source$label,
      species = source$species,
      species_label = upload_run_display_species(source$species),
      file_count = length(source$files),
      info = source$info
    )
    if (isTRUE(include_files)) {
      entry$files <- lapply(source$files, function(file) {
        list(
          key = file$key,
          label = file$label,
          accession = file$accession,
          sample = file$sample,
          treatment = file$treatment,
          caller = file$caller
        )
      })
    }
    entry
  })
}

upload_run_database_selected_files <- function(species, source_id, file_keys, catalog = upload_run_database_catalog()) {
  source_matches <- Filter(
    function(source) {
      identical(source$id, as.character(source_id)) &&
        identical(source$species, as.character(species))
    },
    catalog
  )

  if (!length(source_matches)) {
    stop("The selected database is not available on this server.")
  }

  source <- source_matches[[1]]
  file_keys <- unique(as.character(unlist(file_keys, use.names = FALSE)))
  file_keys <- file_keys[nzchar(file_keys)]
  if (!length(file_keys)) {
    stop("Select at least one database BED file.")
  }

  source_keys <- vapply(source$files, function(file) file$key, character(1))
  missing_keys <- setdiff(file_keys, source_keys)
  if (length(missing_keys)) {
    stop("One or more selected database BED files are no longer available.")
  }

  selected <- source$files[match(file_keys, source_keys)]
  if (!all(vapply(selected, function(file) file.exists(file$file_path), logical(1)))) {
    stop("One or more selected database BED files could not be read from disk.")
  }

  list(
    source_label = source$label,
    files = lapply(selected, function(file) {
      list(
        file_name = file$label,
        file_path = file$file_path,
        source = "Database",
        source_label = source$label
      )
    })
  )
}

upload_run_count_bed_rows <- function(file_path, chunk_size = 100000L) {
  connection <- file(file_path, open = "r")
  on.exit(close(connection), add = TRUE)

  total <- 0L
  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    total <- total + length(lines)
    if (length(lines) < chunk_size) {
      break
    }
  }

  total
}

upload_run_database_bed_preview_page <- function(species, source_id, file_key, page = 1L, page_size = 10L, row_count_cache = NULL, catalog = NULL) {
  selected <- upload_run_database_selected_files(
    species = species,
    source_id = source_id,
    file_keys = file_key,
    catalog = if (is.null(catalog)) upload_run_database_catalog() else catalog
  )
  file_path <- selected$files[[1]]$file_path
  page <- max(1L, suppressWarnings(as.integer(page %||% 1L)))
  page_size <- max(1L, suppressWarnings(as.integer(page_size %||% 10L)))
  cache_key <- normalizePath(file_path, winslash = "/", mustWork = TRUE)
  total_rows <- NULL
  if (!is.null(row_count_cache) && exists(cache_key, envir = row_count_cache, inherits = FALSE)) {
    total_rows <- get(cache_key, envir = row_count_cache, inherits = FALSE)
  }
  if (is.null(total_rows)) {
    total_rows <- upload_run_count_bed_rows(file_path)
    if (!is.null(row_count_cache)) {
      assign(cache_key, total_rows, envir = row_count_cache)
    }
  }

  total_pages <- max(1L, ceiling(total_rows / page_size))
  page <- min(page, total_pages)
  preview_tbl <- read_bed_preview(
    file_path,
    n = page_size,
    skip = (page - 1L) * page_size
  )

  list(
    species = species,
    source_id = source_id,
    file_key = file_key,
    file_name = basename(file_path),
    columns = if (is.null(preview_tbl)) character() else names(preview_tbl),
    rows = if (is.null(preview_tbl)) list() else lapply(seq_len(nrow(preview_tbl)), function(index) as.list(preview_tbl[index, , drop = FALSE])),
    page = page,
    page_size = page_size,
    total = total_rows
  )
}

upload_run_species_display_map <- function() {
  c(
    "aca.ASM1452938v2" = "Adiantum capillus-veneris (ASM1452938v2)",
    "aca.v2" = "Amphidinium carterae (v2)",
    "aga.AgamP4" = "Anopheles gambiae (AgamP4)",
    "ara_TAIR10" = "Arabidopsis thaliana (TAIR10)",
    "atr.v1" = "Acer truncatum (v1)",
    "awo_gca_000247605" = "Acetobacterium woodii (DSM 1030, ASM24760v1)",
    "Basmati1.IGDBv1" = "Oryza sativa (Basmati1 IGDBv1)",
    "bdi.v3.0" = "Brachypodium distachyon (bdi.v3.0)",
    "bna.AST_PRJEB5043_v1" = "Brassica napus (AST_PRJEB5043_v1)",
    "bra.Brapa_1.0" = "Brassica rapa (Brapa 1.0)",
    "bsu_168" = "Bacillus subtilis (str. 168)",
    "bta.ARS-UCD1.2" = "Bos taurus (ARS-UCD1.2)",
    "bth_gca_000011065" = "Bacteroides thetaiotaomicron (VPI-5482, ASM1106v1)",
    "cac_gca_001042715" = "Clostridium aceticum (DSM 1496, GCA_001042715)",
    "cdr_gca_003096175" = "Clostridium drakei (DSM 1496, GCA_003096175)",
    "cel_WBcel235" = "Caenorhabditis elegans (WBcel235)",
    "CG14.IGDBv1" = "Oryza sativa (CG14 IGDBv1)",
    "clj_gca_000143685" = "Clostridium ljungdahlii (DSM 13528, ASM14368v1)",
    "clu.ROS_Cfam_1.0" = "Canis lupus familiaris (ROS Cfam 1.0)",
    "cre.v5.5" = "Chlamydomonas reinhardtii (v5.5)",
    "csa_chineseLong_v3" = "Cucumis sativus chineseLong v3 (chineseLong_v3)",
    "csa.ASM407v2" = "Cucumis sativus (ASM407v2)",
    "csi.AHAU_CSS_1" = "Camellia sinensis (AHAU_CSS_1)",
    "cvi_gca_000022005" = "Caulobacter vibrioides (NA1000, ASM2200v1)",
    "dcug.v1" = "Dalbergia cultrata (v1)",
    "dme_BDGP6" = "Drosophila melanogaster (BDGP6)",
    "dre_GRCz11" = "Danio rerio (GRCz11)",
    "eco_gca_000931565" = "Escherichia coli (K-12 substr. RV308, ASM93156v1)",
    "ecoli_k12" = "Escherichia coli (K-12)",
    "eli_gca_001481725" = "Eubacterium limosum (SA11, ASM148172v1)",
    "fjo_gca_000016645" = "Flavobacterium johnsoniae (UW101, ASM1664v1)",
    "G46.IGDBv1" = "Oryza sativa (G46 IGDBv1)",
    "gga.GRCg7b" = "Gallus gallus (GRCg7b)",
    "ggo.gorGor4" = "Gorilla gorilla (gorGor4)",
    "gma_v2" = "Glycine max (v2.1)",
    "gra.Graimondii2_0_v6" = "Gossypium raimondii (Graimondii2_0_v6)",
    "hg19" = "Homo sapiens (hg19)",
    "hg38" = "Homo sapiens (hg38)",
    "hg38.gencode_v46" = "Homo sapiens (hg38, GENECODE V46)",
    "hsa_gca_000006805" = "Halobacterium salinarum (NRC-1, ASM680v1)",
    "hsa_NRC1" = "Halobacterium salinarum (NRC1)",
    "hvo_DS2" = "Haloferax volcanii (DS2, ASM2568v1)",
    "IR64.IGDBv1" = "Oryza sativa (IR64 IGDBv1)",
    "kmi_gca_000963575" = "Klebsiella michiganensis (RC10, ASM96357v1)",
    "lch.v1" = "Litchi chinensis (v1)",
    "Lemont.IGDBv1" = "Oryza sativa (Lemont IGDBv1)",
    "lin_gca_000195795" = "Listeria innocua (Clip11262, ASM19579v1)",
    "LJ.IGDBv1" = "Oryza sativa (LJ IGDBv1)",
    "lmi.240118" = "Locusta migratoria (NCBI V2)",
    "lmi_cau" = "Locusta migratoria (CAU V1)",
    "lmo_gca_000196035" = "Listeria monocytogenes (EGD-e, ASM19603v1)",
    "lrh_gca_002287945" = "Lacticaseibacillus rhamnosus GG (DSM 14870, GCA_002287945)",
    "mab_gca_000069185" = "Mycobacteroides abscessus (ATCC19977, ASM6918v1)",
    "mdo.ASM211411v1" = "Malus domestica golden (ASM211411v1)",
    "mma.ASM706v1" = "Methanosarcina mazei (ASM706v1)",
    "mm10" = "Mus musculus (mm10)",
    "mm39" = "Mus musculus (mm39)",
    "mml.Mmul_10" = "Macaca mulatta (Mmul_10)",
    "mpo.MpTak1_v7.1" = "Marchantia polymorpha (MpTak1 v7.1)",
    "MSU.IGDBv1" = "Oryza sativa (MSU IGDBv1)",
    "msm_gca_000015005" = "Mycolicibacterium smegmatis (MC2-155, ASM1500v1)",
    "mtr.MedtrA17_4.0" = "Medicago truncatula (MedtrA17 4.0)",
    "mtu_gca_000016145" = "Mycobacterium tuberculosis (H37Ra, ASM1614v1)",
    "N22.IGDBv1" = "Oryza sativa (N22 IGDBv1)",
    "NamRoo.IGDBv1" = "Oryza sativa (NamRoo IGDBv1)",
    "nbe.v1" = "Nicotiana benthamiana (v1)",
    "ncr.NC12" = "Neurospora crassa (NC12)",
    "ntab.TN90" = "Nicotiana tabacum (TN90)",
    "oru.OR_W1943" = "Oryza rufipogon (OR_W1943)",
    "osa_rufipogon" = "Oryza sativa (Rufipogon)",
    "osai.ASM465v1" = "Oryza indica (ASM465v1)",
    "osa_huazhan" = "Oryza sativa (Huazhan 1.0)",
    "osa_IRGSP_1" = "Oryza sativa (IRGSP 1.0, RAP-DB version)",
    "osa_msu7" = "Oryza sativa (IRGSP 1.0, MSU7 version)",
    "pae_gca_002968755" = "Pseudomonas aeruginosa (AR441, GCA_002968755)",
    "pfl_gca_001307155" = "Pseudomonas fluorescens (FW300-N2E3, ASM130715v1)",
    "pfu_dsm_3638" = "Pyrococcus furiosus (DSM 3638)",
    "ppa.Phypa_V3" = "Physcomitrium patens (Phypa_V3)",
    "ptc.Pop_tri_v3" = "Populus trichocarpa (Pop_tri_v3)",
    "ptr.Pan_tro_3.0" = "Pan troglodytes (Pan_tro_3.0)",
    "rn5" = "Rattus norvegicus (rn5)",
    "rn6" = "Rattus norvegicus (rn6)",
    "sau_gca_001018655" = "Staphylococcus aureus (GCA_001018655)",
    "sav_gca_000009765" = "Streptomyces avermitilis (MA4680, ASM976v2)",
    "sbi.NCBIv3" = "Sorghum bicolor (NCBIv3)",
    "sce_R64" = "Saccharomyces cerevisiae (R64)",
    "scl_GCF_005519465" = "Streptomyces clavuligerus (ATCC 27064, GCF_005519465)",
    "sco_GCF_000203835" = "Streptomyces coelicolor (A3(2), ASM20383v1)",
    "sen_gca_003325055" = "Salmonella enterica (GCA_003325055)",
    "sgr_gca_000010605" = "Streptomyces griseus (JCM 4626, GCA_000010605)",
    "ska.Fugka2468_1" = "Symbiodinium kawagutii (Fugka2468_1)",
    "sli_GCF_000739105" = "Streptomyces lividans (TK24, GCF_000739105)",
    "sly_SL3" = "Solanum lycopersicum (SL v3.0)",
    "sme_gca_002197125" = "Sinorhizobium meliloti (M162, ASM219712v1)",
    "smo.v1" = "Selaginella moellendorffii (v1.0)",
    "spn_gca_000018965" = "Streptococcus pneumoniae (70585, ASM1896v1)",
    "ssc.Sscrofa11.1" = "Sus scrofa (Sscrofa11.1)",
    "sts_GCF_003932715" = "Streptomyces tsukubensis (NRRL 18488, GCF_003932715)",
    "stu.SolTub_3.0" = "Solanum tuberosum (SolTub 3.0)",
    "sve_gca_000253235" = "Streptomyces venezuelae (ATCC 10712, ASM25323v1)",
    "syn_gca_000478825" = "Synechocystis sp. (PCC 6714, ASM47882v2)",
    "tae.IWGSC" = "Triticum aestivum (IWGSC RefSeq 1.0)",
    "TM.IGDBv1" = "Oryza sativa (TM IGDBv1)",
    "Tumba.IGDBv1" = "Oryza sativa (Tumba IGDBv1)",
    "vna_gca_001680085" = "Vibrio natriegens (CCUG 16374, ASM168008v1)",
    "vvi.PN40024.v4" = "Vitis vinifera (PN40024.v4)",
    "vvu_gca_000039765" = "Vibrio vulnificus (CMCP6, ASM3976v1)",
    "WSSM.IGDBv1" = "Oryza sativa (WSSM IGDBv1)",
    "xtr.UCB_Xtro_10.0" = "Xenopus tropicalis (UCB Xtro 10.0)",
    "zma.B73_NAM5.0" = "Zea mays (B73 NAM 5.0)",
    "zma_RefGen_v4" = "Zea mays (B73 RefGen v4)",
    "zmo_gca_000007105" = "Zymomonas mobilis (ZM4, ASM710v1)"
  )
}

upload_run_annotation_codes_for_extension <- function(dir_path, extension) {
  if (!dir.exists(dir_path)) {
    return(character())
  }

  files <- list.files(dir_path, pattern = paste0("\\Q", extension, "\\E$"))
  sub(paste0("\\Q", extension, "\\E$"), "", files)
}

upload_run_species_label_from_code <- function(code) {
  display_map <- upload_run_species_display_map()
  label <- unname(display_map[code])

  if (!is.null(label) && nzchar(label)) {
    return(label)
  }

  code
}

upload_run_species_choice <- function(code) {
  list(
    value = code,
    label = upload_run_species_label_from_code(code)
  )
}

upload_run_choice_value <- function(choice) {
  if (is.list(choice) && !is.null(choice$value)) {
    return(as.character(choice$value))
  }

  as.character(choice)
}

upload_run_display_species <- function(species) {
  if (is.null(species) || !nzchar(species)) {
    return(species)
  }

  upload_run_species_label_from_code(species)
}

read_species_choices <- function() {
  data_root <- upload_run_annotation_data_root()
  txdb_codes <- upload_run_annotation_codes_for_extension(file.path(data_root, "txdb"), ".txdb.sqlite")
  txlens_codes <- upload_run_annotation_codes_for_extension(file.path(data_root, "txlens"), ".txlens.rda")
  gff_codes <- upload_run_annotation_codes_for_extension(file.path(data_root, "gff"), ".gff.rda")
  available_codes <- sort(Reduce(intersect, list(txdb_codes, txlens_codes, gff_codes)))

  if (!length(available_codes)) {
    return(character())
  }

  labels <- vapply(available_codes, upload_run_species_label_from_code, character(1))
  ordered_codes <- available_codes[order(tolower(labels))]

  lapply(ordered_codes, upload_run_species_choice)
}

upload_run_demo_species <- function(species_choices) {
  preferred_demo_species <- "ara_TAIR10"
  choice_values <- vapply(species_choices, upload_run_choice_value, character(1))

  if (preferred_demo_species %in% choice_values) {
    return(preferred_demo_species)
  }

  if (length(species_choices)) {
    return(upload_run_choice_value(species_choices[[1]]))
  }

  ""
}

upload_run_select_config <- function(id, database_catalog = upload_run_database_catalog()) {
  ns <- NS(id)
  species_choices <- read_species_choices()

  list(
    data_source_id = ns("data_source"),
    species_id = ns("species"),
    bed_file_id = ns("bed_file"),
    save_button_id = ns("save_upload"),
    database_save_input = ns("save_database_selection"),
    database_context_input = ns("database_context"),
    demo_button_id = ns("load_demo"),
    demo_species = upload_run_demo_species(species_choices),
    demo_file_name = "ac_diff_down.bed6.bed",
    data_source_choices = c("Upload File", "Select from Database", "Upload + Database"),
    species_choices = species_choices,
    database_catalog = upload_run_database_control_catalog(database_catalog, include_files = FALSE)
  )
}
