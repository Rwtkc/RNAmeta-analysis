# RNAmeta core analysis algorithms

Core analysis source code for RNAmeta Online, Meta Plot and RNAenrich.

## Links

- RNAmeta website and downloads: https://rnainformatics.cn/RNAmeta/
- R package: https://github.com/Rwtkc/RNAmeta_R
- Desktop application: https://github.com/Rwtkc/RNAmeta_desktop

## Source code

| Directory or file | Analysis |
| --- | --- |
| `online/peak_distribution.R` | Genomic feature construction, peak annotation and regional counts |
| `online/gene_matrix.R` | Peak-associated genes and gene-set comparison inputs |
| `online/gene_type.R` | RNA biotype annotation and counts |
| `online/exon_features.R` | Exon position, exon length and transcript exon counts |
| `online/transcript_length.R` | Lengths of peak-associated transcripts |
| `online/transcription.R` | Transcription start and end windows |
| `online/translation.R` | Translation start and end windows in transcript coordinates |
| `online/splice_sites.R` | Strand-aware 5′ and 3′ splice-boundary windows |
| `online/boundary_profiles.R` | Boundary mapping, density estimation and window-level coverage matrices |
| `online/motif/` | Sequence extraction, matched non-peak backgrounds, HOMER/MEME/STREME arguments and motif parsing |
| `online/browser.R` | Genomic-to-transcript coordinate mapping and annotation projection |
| `online/structure.R` | Transcript sequence extraction and RNAfold structure prediction |
| `online/shared.R` | Shared reference-path and value helpers |
| `meta_plot/` | Transcript components, coordinate normalization, weighted density estimation and bootstrap intervals |
| `rnaenrich/rnaenrich-core.ts` | BED normalization, exact and interval overlap, Jaccard indices, intersection membership, hypergeometric tests and BH adjustment |

The R functions are grouped by analysis, with shared computational and caching dependencies retained. The analysis functions require the corresponding R/Bioconductor packages, including data.table, GenomicRanges, GenomicFeatures, rtracklayer, Biostrings and Rsamtools. Motif discovery uses HOMER, MEME or STREME; structure prediction uses RNAfold. RNAenrich's computational core is written in TypeScript and has no React dependency.

Reference annotations and sequence files are available through the RNAmeta website. When using the R source directly, `options(rnameta.reference_root = "/path/to/reference-root")` selects the directory containing `data/txdb`, `data/gff`, `data/txlens` and `data/fa`. RNAfold is resolved through `RNAMETA_RNAFOLD` or the system `PATH`. The R package and Desktop repositories provide the packaged analysis interfaces.

## License

RNAmeta source code is distributed under GPL-3.0-or-later. External packages and tools retain their own licenses.
